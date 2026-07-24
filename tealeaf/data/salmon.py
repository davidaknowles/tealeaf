"""Salmon equivalence-class weights for fixed GLM observation designs."""

from __future__ import annotations

from dataclasses import dataclass
import gzip
import json
from pathlib import Path
import struct

import numpy as np
import scipy.sparse as sp


def _open_text(path):
    path = Path(path)
    return gzip.open(path, "rt") if path.suffix == ".gz" else path.open()


@dataclass
class WeightedEC:
    weighted_sum: np.ndarray
    count: float
    rows: int = 1

    @property
    def weights(self):
        return self.weighted_sum / self.count


def read_salmon_weighted_eqclasses(path, feature_to_index):
    """Read and count-collapse Salmon rich ECs keyed by global transcript IDs."""
    feature_to_index = dict(feature_to_index)
    by_key = {}
    with _open_text(path) as handle:
        num_targets = int(next(handle))
        num_rows = int(next(handle))
        names = [next(handle).rstrip("\n") for _ in range(num_targets)]
        missing = [name for name in names if name not in feature_to_index]
        if missing:
            raise ValueError(
                f"{len(missing)} Salmon targets are absent from alevin features; "
                f"first: {missing[0]}"
            )
        local_to_global = np.asarray(
            [feature_to_index[name] for name in names], dtype=np.int64
        )
        observed_rows = 0
        for line_number, line in enumerate(handle, start=num_targets + 3):
            fields = line.split()
            if not fields:
                continue
            size = int(fields[0])
            if len(fields) != 2 * size + 2:
                raise ValueError(
                    f"Salmon EC line {line_number} lacks rich weights"
                )
            local_ids = np.asarray(fields[1:1 + size], dtype=np.int64)
            weights = np.asarray(
                fields[1 + size:1 + 2 * size], dtype=np.float64
            )
            count = float(fields[-1])
            if count <= 0 or np.any(weights < 0) or not np.all(np.isfinite(weights)):
                raise ValueError(f"invalid Salmon EC values on line {line_number}")
            global_ids = local_to_global[local_ids]
            order = np.argsort(global_ids)
            key = tuple(global_ids[order])
            weights = weights[order]
            total = weights.sum()
            if total <= 0:
                raise ValueError(f"zero Salmon EC weight on line {line_number}")
            weights /= total
            contribution = count * weights
            existing = by_key.get(key)
            if existing is None:
                by_key[key] = WeightedEC(contribution, count)
            else:
                existing.weighted_sum += contribution
                existing.count += count
                existing.rows += 1
            observed_rows += 1
    if observed_rows != num_rows:
        raise ValueError(
            f"Salmon header reports {num_rows} EC rows, parsed {observed_rows}"
        )
    return by_key, {
        "salmon_targets": num_targets,
        "salmon_rows": observed_rows,
        "salmon_unique_transcript_sets": len(by_key),
        "salmon_duplicate_rows": observed_rows - len(by_key),
    }


def read_salmon_effective_lengths(path, feature_to_index):
    """Read primer-specific bias-corrected effective lengths from quant.sf."""
    lengths = np.full(len(feature_to_index), np.nan, dtype=np.float64)
    with open(path) as handle:
        header = next(handle).rstrip("\n").split("\t")
        name_col = header.index("Name")
        length_col = header.index("EffectiveLength")
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            index = feature_to_index.get(fields[name_col])
            if index is not None:
                lengths[index] = max(float(fields[length_col]), 1.0)
    return lengths


def summarize_positional_bias_models(path):
    """Summarize Salmon's normalized 20-bin positional models."""
    with gzip.open(path, "rb") as handle:
        num_models = struct.unpack("<I", handle.read(4))[0]
        bounds = list(struct.unpack(
            f"<{num_models}I", handle.read(4 * num_models)
        ))
        models = []
        for _ in range(num_models):
            size = struct.unpack("<I", handle.read(4))[0]
            values = np.frombuffer(handle.read(8 * size), dtype="<f8").copy()
            models.append(values)
        if handle.read(1):
            raise ValueError(f"unexpected trailing data in positional model {path}")
    centers = (np.arange(models[0].size) + 0.5) / models[0].size
    means = [float(np.dot(values, centers) / values.sum()) for values in models]
    return {
        "length_class_upper_bounds": bounds,
        "relative_position_means": means,
    }


def build_positional_ec_design(
    membership,
    features,
    salmon_eqclasses,
    salmon_quant,
    *,
    ec_counts=None,
    fallback_inverse_lengths=None,
):
    """Map Salmon positional-bias weights onto an alevin EC universe.

    Duplicate range-factorized Salmon rows are collapsed by raw fragment count.
    ECs absent from the bulk Salmon dump use inverse primer-specific effective
    lengths. The resulting conditional weights are normalized by transcript
    column to produce the fixed observation matrix used by the GLM.
    """
    membership = membership.tocsr()
    membership.sort_indices()
    features = list(features)
    feature_to_index = {name: index for index, name in enumerate(features)}
    if len(feature_to_index) != len(features):
        raise ValueError("alevin transcript features are not unique")
    eqclass_paths = (
        list(salmon_eqclasses)
        if isinstance(salmon_eqclasses, (list, tuple))
        else [salmon_eqclasses]
    )
    quant_paths = (
        [Path(path) for path in salmon_quant]
        if isinstance(salmon_quant, (list, tuple))
        else [Path(salmon_quant)]
    )
    if len(eqclass_paths) != len(quant_paths):
        raise ValueError("Salmon EC and quantification path counts differ")
    salmon_ecs = {}
    stats = {
        "salmon_targets": len(features),
        "salmon_rows": 0,
        "salmon_duplicate_rows": 0,
        "salmon_quantifications": len(quant_paths),
    }
    effective_sum = np.zeros(len(features), dtype=np.float64)
    effective_weight = np.zeros(len(features), dtype=np.float64)
    for eqclasses, quant in zip(eqclass_paths, quant_paths):
        run_ecs, run_stats = read_salmon_weighted_eqclasses(
            eqclasses, feature_to_index
        )
        stats["salmon_rows"] += run_stats["salmon_rows"]
        stats["salmon_duplicate_rows"] += run_stats["salmon_duplicate_rows"]
        for key, weighted in run_ecs.items():
            existing = salmon_ecs.get(key)
            if existing is None:
                salmon_ecs[key] = weighted
            else:
                existing.weighted_sum += weighted.weighted_sum
                existing.count += weighted.count
                existing.rows += weighted.rows
        lengths = read_salmon_effective_lengths(
            quant / "quant.sf", feature_to_index
        )
        meta_path = quant / "aux_info" / "meta_info.json"
        meta = json.loads(meta_path.read_text())
        run_weight = max(float(meta.get("num_mapped", 0)), 1.0)
        finite = np.isfinite(lengths)
        effective_sum[finite] += run_weight * lengths[finite]
        effective_weight[finite] += run_weight
    stats["salmon_unique_transcript_sets"] = len(salmon_ecs)
    stats["salmon_duplicate_rows"] = (
        sum(weighted.rows for weighted in salmon_ecs.values()) - len(salmon_ecs)
    )
    effective_lengths = np.full(len(features), np.nan, dtype=np.float64)
    known = effective_weight > 0
    effective_lengths[known] = effective_sum[known] / effective_weight[known]
    if fallback_inverse_lengths is None:
        fallback_inverse_lengths = np.ones(len(features), dtype=np.float64)
    fallback_inverse_lengths = np.asarray(fallback_inverse_lengths, dtype=np.float64)
    if fallback_inverse_lengths.shape != (len(features),):
        raise ValueError("fallback_inverse_lengths has the wrong shape")
    inverse_lengths = fallback_inverse_lengths.copy()
    known_lengths = np.isfinite(effective_lengths)
    inverse_lengths[known_lengths] = 1.0 / effective_lengths[known_lengths]

    data = np.empty(membership.nnz, dtype=np.float32)
    matched = np.zeros(membership.shape[0], dtype=bool)
    duplicate_matches = 0
    for ec_id in range(membership.shape[0]):
        start, stop = membership.indptr[ec_id:ec_id + 2]
        indices = membership.indices[start:stop]
        weighted = salmon_ecs.get(tuple(indices))
        if weighted is not None:
            values = weighted.weights
            matched[ec_id] = True
            duplicate_matches += weighted.rows > 1
        else:
            values = inverse_lengths[indices]
            if not np.all(np.isfinite(values)) or values.sum() <= 0:
                raise ValueError(f"no valid fallback weights for alevin EC {ec_id}")
            values = values / values.sum()
        data[start:stop] = values

    design = sp.csr_matrix(
        (data, membership.indices.copy(), membership.indptr.copy()),
        shape=membership.shape,
    )
    column_sums = np.asarray(design.sum(axis=0)).ravel()
    inverse = np.zeros_like(column_sums)
    inverse[column_sums > 0] = 1.0 / column_sums[column_sums > 0]
    design = (design @ sp.diags(inverse)).tocsr()

    stats.update({
        "alevin_ecs": int(membership.shape[0]),
        "matched_alevin_ecs": int(matched.sum()),
        "fallback_alevin_ecs": int((~matched).sum()),
        "duplicate_transcript_set_matches": int(duplicate_matches),
        "matched_ec_fraction": float(matched.mean()),
        "effective_lengths_available": int(known_lengths.sum()),
    })
    if ec_counts is not None:
        ec_counts = np.asarray(ec_counts).ravel()
        if ec_counts.shape != (membership.shape[0],):
            raise ValueError("ec_counts has the wrong shape")
        total = float(ec_counts.sum())
        matched_mass = float(ec_counts[matched].sum())
        stats.update({
            "alevin_molecule_count": total,
            "matched_molecule_count": matched_mass,
            "matched_molecule_fraction": matched_mass / total if total else 0.0,
        })
    return design, stats, effective_lengths
