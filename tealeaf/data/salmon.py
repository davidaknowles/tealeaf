"""Salmon equivalence-class weights for fixed GLM observation designs."""

from __future__ import annotations

from dataclasses import dataclass
import gzip
import json
from pathlib import Path
import struct

import numpy as np
import scipy.sparse as sp


def salmon_eqclass_paths(quant_paths):
    return [
        Path(path) / "aux_info" / "eq_classes.txt.gz"
        for path in quant_paths
    ]


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


def read_salmon_weighted_eqclasses(
    path,
    feature_to_index,
    *,
    allowed_keys=None,
):
    """Read and count-collapse Salmon rich ECs keyed by global transcript IDs."""
    feature_to_index = dict(feature_to_index)
    if allowed_keys is not None and not isinstance(allowed_keys, (set, frozenset)):
        allowed_keys = set(allowed_keys)
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
        retained_rows = 0
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
            observed_rows += 1
            if allowed_keys is not None and key not in allowed_keys:
                continue
            contribution = count * weights
            existing = by_key.get(key)
            if existing is None:
                by_key[key] = WeightedEC(contribution, count)
            else:
                existing.weighted_sum += contribution
                existing.count += count
                existing.rows += 1
            retained_rows += 1
    if observed_rows != num_rows:
        raise ValueError(
            f"Salmon header reports {num_rows} EC rows, parsed {observed_rows}"
        )
    return by_key, {
        "salmon_targets": num_targets,
        "salmon_rows": observed_rows,
        "salmon_retained_rows": retained_rows,
        "salmon_unique_transcript_sets": len(by_key),
        "salmon_duplicate_rows": retained_rows - len(by_key),
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


def validate_primer_positional_quantification(
    run_root,
    source_meta,
    *,
    expected_targets=None,
    expected_library_type=None,
):
    """Validate read conservation and Salmon products for one Parse run."""
    run_root = Path(run_root)
    source_meta = json.loads(Path(source_meta).read_text())
    stats = json.loads((run_root / "demultiplex_stats.json").read_text())
    count_fields = (
        "total",
        "polydT_exact",
        "polydT_corrected",
        "ranhex_exact",
        "ranhex_corrected",
        "unknown_or_ambiguous",
        "assigned",
    )
    if any(
        not isinstance(stats.get(field), int) or stats[field] < 0
        for field in count_fields
    ):
        raise ValueError("demultiplex statistics contain invalid counts")
    expected_total = int(source_meta["num_processed"])
    if stats["total"] != expected_total:
        raise ValueError(
            f"demultiplexed {stats['total']} reads; expected {expected_total}"
        )
    assigned = sum(
        stats[f"{primer}_{match}"]
        for primer in ("polydT", "ranhex")
        for match in ("exact", "corrected")
    )
    if stats["assigned"] != assigned:
        raise ValueError("assigned read count does not equal primer components")
    if assigned + stats["unknown_or_ambiguous"] != stats["total"]:
        raise ValueError("primer assignments do not conserve input reads")

    report = {
        "total_reads": expected_total,
        "assigned_reads": assigned,
        "assigned_fraction": assigned / expected_total if expected_total else 0.0,
        "unknown_or_ambiguous_reads": stats["unknown_or_ambiguous"],
        "primers": {},
    }
    for primer, stats_name in (("polydt", "polydT"), ("ranhex", "ranhex")):
        quant = run_root / f"salmon_{primer}"
        required = (
            quant / ".tealeaf_complete",
            quant / "quant.sf",
            quant / "aux_info" / "eq_classes.txt.gz",
            quant / "aux_info" / "meta_info.json",
            quant / "aux_info" / "obs5_pos.gz",
            quant / "aux_info" / "obs3_pos.gz",
            quant / "aux_info" / "exp5_pos.gz",
            quant / "aux_info" / "exp3_pos.gz",
        )
        missing = [str(path) for path in required if not path.is_file() or not path.stat().st_size]
        if missing:
            raise FileNotFoundError(
                f"missing positional Salmon products: {', '.join(missing)}"
            )
        meta = json.loads((quant / "aux_info" / "meta_info.json").read_text())
        library_types = meta.get("library_types", [])
        if (
            expected_library_type is not None
            and library_types != [expected_library_type]
        ):
            raise ValueError(
                f"{primer} library type is {library_types}; "
                f"expected {[expected_library_type]}"
            )
        primer_reads = stats[f"{stats_name}_exact"] + stats[f"{stats_name}_corrected"]
        if int(meta["num_processed"]) != primer_reads:
            raise ValueError(
                f"{primer} Salmon processed {meta['num_processed']} reads; "
                f"demultiplexer assigned {primer_reads}"
            )
        mapped = int(meta["num_mapped"])
        if mapped <= 0 or mapped > primer_reads:
            raise ValueError(f"{primer} Salmon mapped count is invalid")
        with gzip.open(quant / "aux_info" / "eq_classes.txt.gz", "rt") as handle:
            target_count = int(next(handle))
            eq_count = int(next(handle))
        if target_count <= 0 or eq_count <= 0:
            raise ValueError(f"{primer} rich EC header is empty")
        if expected_targets is not None and target_count != int(expected_targets):
            raise ValueError(
                f"{primer} has {target_count} targets; expected {expected_targets}"
            )
        with open(quant / "quant.sf") as handle:
            header = next(handle).rstrip("\n").split("\t")
            if header[:3] != ["Name", "Length", "EffectiveLength"]:
                raise ValueError(f"{primer} quant.sf has an unexpected header")
            quant_targets = sum(1 for _ in handle)
        if quant_targets != target_count:
            raise ValueError(
                f"{primer} quant.sf has {quant_targets} rows; expected {target_count}"
            )
        for model in ("obs5_pos.gz", "obs3_pos.gz", "exp5_pos.gz", "exp3_pos.gz"):
            summarize_positional_bias_models(quant / "aux_info" / model)
        report["primers"][primer] = {
            "processed_reads": primer_reads,
            "mapped_reads": mapped,
            "mapping_rate": mapped / primer_reads,
            "rich_equivalence_classes": eq_count,
            "targets": target_count,
            "library_types": library_types,
        }
    return report


def summarize_primer_positional_validations(paths):
    """Aggregate validated Parse primer runs into one audit summary."""
    paths = [Path(path) for path in paths]
    if not paths:
        raise ValueError("at least one positional validation report is required")
    reports = []
    run_ids = []
    for path in paths:
        report = json.loads(path.read_text())
        run_id = path.parent.name
        if run_id in run_ids:
            raise ValueError(f"duplicate positional validation run: {run_id}")
        run_ids.append(run_id)
        reports.append(report)
    summary = {
        "runs": len(reports),
        "run_ids": sorted(run_ids),
        "total_reads": sum(int(report["total_reads"]) for report in reports),
        "assigned_reads": sum(int(report["assigned_reads"]) for report in reports),
        "unknown_or_ambiguous_reads": sum(
            int(report["unknown_or_ambiguous_reads"]) for report in reports
        ),
        "primers": {},
    }
    if (
        summary["assigned_reads"] + summary["unknown_or_ambiguous_reads"]
        != summary["total_reads"]
    ):
        raise ValueError("aggregate primer validation does not conserve reads")
    summary["assigned_fraction"] = (
        summary["assigned_reads"] / summary["total_reads"]
        if summary["total_reads"] else 0.0
    )
    for primer in ("polydt", "ranhex"):
        primer_reports = [report["primers"][primer] for report in reports]
        processed = sum(int(report["processed_reads"]) for report in primer_reports)
        mapped = sum(int(report["mapped_reads"]) for report in primer_reports)
        rates = [float(report["mapping_rate"]) for report in primer_reports]
        summary["primers"][primer] = {
            "processed_reads": processed,
            "mapped_reads": mapped,
            "mapping_rate": mapped / processed if processed else 0.0,
            "minimum_run_mapping_rate": min(rates),
            "maximum_run_mapping_rate": max(rates),
            "rich_equivalence_classes_sum": sum(
                int(report["rich_equivalence_classes"])
                for report in primer_reports
            ),
            "targets": sorted({
                int(report["targets"]) for report in primer_reports
            }),
        }
    return summary


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
    allowed_keys = {
        tuple(membership.indices[start:stop])
        for start, stop in zip(membership.indptr[:-1], membership.indptr[1:])
    }
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
        "salmon_retained_rows": 0,
        "salmon_duplicate_rows": 0,
        "salmon_quantifications": len(quant_paths),
        "alevin_unique_transcript_sets": len(allowed_keys),
    }
    effective_sum = np.zeros(len(features), dtype=np.float64)
    effective_weight = np.zeros(len(features), dtype=np.float64)
    for eqclasses, quant in zip(eqclass_paths, quant_paths):
        run_ecs, run_stats = read_salmon_weighted_eqclasses(
            eqclasses,
            feature_to_index,
            allowed_keys=allowed_keys,
        )
        stats["salmon_rows"] += run_stats["salmon_rows"]
        stats["salmon_retained_rows"] += run_stats["salmon_retained_rows"]
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
