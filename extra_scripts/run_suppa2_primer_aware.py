#!/usr/bin/env python3
"""Run a primer-aware SUPPA2 event sensitivity analysis.

The event definitions are the standard SUPPA2 catalogue. Tealeaf's paired
primer preparation supplies primer-specific EC probabilities and counts; those
are aggregated to subject-by-cell-type samples and converted to expected
included/excluded event counts. A shared event logit with a primer offset is
estimated using beta-binomial-inspired iteratively reweighted logits, followed
by the same paired Wilcoxon approximation and within-contrast BH correction as
the pooled SUPPA2 comparator.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import sparse
from scipy.special import expit, logit

from tealeaf.sc import glm_cv

try:
    from extra_scripts.run_suppa2_full_data_comparison import bh, fast_paired_wilcoxon
except ModuleNotFoundError:
    from run_suppa2_full_data_comparison import bh, fast_paired_wilcoxon


EPS = 1e-6


def canonical(value):
    return str(value).split(".", 1)[0]


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--alevin-dir", required=True, type=Path)
    parser.add_argument("--salmon-ref", required=True, type=Path)
    parser.add_argument("--primer-pairs", required=True, type=Path)
    parser.add_argument("--probability-file", type=Path)
    parser.add_argument("--barcode-groups", required=True, type=Path)
    parser.add_argument("--event-catalog", required=True, type=Path)
    parser.add_argument("--contrasts", action="append", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--event-output", required=True, type=Path)
    parser.add_argument("--work-dir", required=True, type=Path)
    parser.add_argument("--fold", default="full_data")
    parser.add_argument("--minimum-pairs", type=int, default=8)
    parser.add_argument("--min-half-umis", type=float, default=500)
    parser.add_argument("--min-eq", type=int, default=5)
    parser.add_argument("--ec-design", default="weighted", choices=("binary", "weighted", "positional"))
    parser.add_argument("--primer-sampling-model", default="oligodt_tpm", choices=("effective_length", "oligodt_tpm", "all_tpm"))
    parser.add_argument("--event-design", default="theta", choices=("raw", "theta"), help="EC-to-transcript allocation used for event PSI; theta is the production abundance-fitting design and raw is an ablation using primer-specific EC-conditional probabilities.")
    parser.add_argument("--max-events", type=int, help="limit the catalogue for a smoke test")
    return parser.parse_args()


def load_catalog(path):
    return pd.read_csv(path, sep="\t", compression="gzip", dtype=str).fillna("")


def event_matrices(catalog, features):
    feature_index = {canonical(name): index for index, name in enumerate(features)}
    included_rows, included_cols = [], []
    excluded_rows, excluded_cols = [], []
    keep = []
    for event_index, row in enumerate(catalog.itertuples(index=False)):
        included = sorted({feature_index[canonical(x)] for x in row.included.split(",") if canonical(x) in feature_index})
        excluded = sorted({feature_index[canonical(x)] for x in row.excluded.split(",") if canonical(x) in feature_index})
        if not included or not excluded:
            continue
        keep.append(event_index)
        included_rows.extend([len(keep) - 1] * len(included))
        included_cols.extend(included)
        excluded_rows.extend([len(keep) - 1] * len(excluded))
        excluded_cols.extend(excluded)
    shape = (len(keep), len(features))
    included = sparse.csr_matrix((np.ones(len(included_rows)), (included_rows, included_cols)), shape=shape)
    excluded = sparse.csr_matrix((np.ones(len(excluded_rows)), (excluded_rows, excluded_cols)), shape=shape)
    return catalog.iloc[keep].reset_index(drop=True), included, excluded


def sample_group(value):
    parts = str(value).rsplit("__", 2)
    if len(parts) != 3:
        return None
    cell_type, _, subject = parts
    return f"{subject}__{cell_type}"


def prepare_event_counts(args, catalog):
    probability_file = args.probability_file or (args.alevin_dir / "gene_eqclass_probs.tsv.gz")
    prepared = glm_cv.prepare_paired_primer_glm_data(
        args.alevin_dir,
        args.salmon_ref,
        args.primer_pairs,
        ec_design=args.ec_design,
        regularization_target="theta",
        min_eq=args.min_eq,
        min_half_umis=args.min_half_umis,
        primer_sampling_model=args.primer_sampling_model,
        probability_file=probability_file,
        retain_event_probabilities=True,
    )
    pairs = pd.read_csv(args.primer_pairs, sep="\t", dtype=str)
    source_rows = np.asarray(prepared.metadata["source_rows"], dtype=int)
    if source_rows.shape[0] != len(prepared.barcodes):
        raise ValueError("prepared source rows do not match paired samples")
    pair_table = pairs.set_index("cell_id")
    poly_barcodes = pair_table.loc[prepared.barcodes, "polydt_barcode"].to_numpy()
    n_ec = prepared.cv_raw_counts.shape[1] // 2
    raw_poly = prepared.cv_raw_counts[:, :n_ec].tocsr()
    raw_hex = prepared.cv_raw_counts[:, n_ec:].tocsr()
    group_table = pd.read_csv(args.barcode_groups, header=None, names=["barcode", "group"], dtype=str)
    barcode_to_group = dict(zip(group_table.barcode, group_table.group))
    sample_names = np.asarray([sample_group(barcode_to_group.get(barcode, "")) for barcode in poly_barcodes], dtype=object)
    valid_samples = np.array([name is not None for name in sample_names])
    if not valid_samples.all():
        print(f"dropping {(~valid_samples).sum():,} paired samples without benchmark metadata")
        sample_names = sample_names[valid_samples]
        raw_poly = raw_poly[valid_samples]
        raw_hex = raw_hex[valid_samples]
    sample_levels = sorted(set(sample_names))
    sample_index = {name: index for index, name in enumerate(sample_levels)}
    sample_ids = np.asarray([sample_index[name] for name in sample_names], dtype=int)
    selector = sparse.coo_matrix((np.ones(len(sample_ids)), (sample_ids, np.arange(len(sample_ids)))), shape=(len(sample_levels), len(sample_ids))).tocsr()
    # Event PSI needs EC-conditional transcript probabilities. The legacy
    # route row-normalized the theta compatibility blocks, which are intended
    # for abundance fitting and include transcript exposure corrections. Use
    # the raw primer-specific probability blocks for event allocation.
    phi = []
    if args.event_design == "raw":
        event_blocks = prepared.metadata.get("event_probability_blocks")
        if event_blocks is None:
            raise ValueError("raw event probabilities were not retained")
    else:
        event_blocks = [prepared.compatibility[start : start + n_ec] for start in (0, n_ec)]
    for design in event_blocks:
        design = design.tocsr()
        row_totals = np.asarray(design.sum(axis=1)).ravel()
        inverse = np.divide(1.0, row_totals, out=np.zeros_like(row_totals), where=row_totals > 0)
        phi.append((sparse.diags(inverse) @ design).tocsr())
    catalog, included, excluded = event_matrices(catalog, prepared.features)
    included_by_transcript = [matrix @ included.T for matrix in phi]
    excluded_by_transcript = [matrix @ excluded.T for matrix in phi]
    counts_by_primer = []
    for raw, included_ec, excluded_ec in zip((raw_poly, raw_hex), included_by_transcript, excluded_by_transcript):
        sample_ec = (selector @ raw).tocsr()
        included_counts = (sample_ec @ included_ec).toarray().astype(float)
        excluded_counts = (sample_ec @ excluded_ec).toarray().astype(float)
        counts_by_primer.append((included_counts, excluded_counts))
    return sample_levels, catalog, counts_by_primer


def fit_primer_aware_psi(included, excluded, offset=None, dispersion=None, *, return_dispersion=False):
    """Estimate shared logits with a shrunk primer offset.

    The offset is a weighted logit difference between primer halves. The
    beta-binomial overdispersion enters the IRLS weights, which downweights
    low-count observations while retaining a shared biological logit.
    """
    totals = included + excluded
    proportions = np.divide(included, totals, out=np.full_like(included, np.nan), where=totals > 0)
    logits = logit(np.clip(proportions, EPS, 1 - EPS))
    if offset is None:
        weights = np.sqrt(np.maximum(totals, 1.0))
        difference = logits[:, 0] - logits[:, 1]
        valid = np.isfinite(difference)
        raw_offset = np.sum(weights[:, 0][valid] * difference[valid]) / np.sum(weights[:, 0][valid]) if valid.any() else 0.0
        # A weak beta-binomial shrinkage stabilizes events with little support.
        support = np.sum(totals[valid]) if valid.any() else 0.0
        offset = float(raw_offset * support / (support + 200.0))
    if dispersion is None:
        finite_proportion = np.isfinite(proportions)
        mean_proportion = np.divide(
            np.nansum(proportions, axis=1),
            finite_proportion.sum(axis=1),
            out=np.full(len(proportions), np.nan),
            where=finite_proportion.sum(axis=1) > 0,
        )
        valid_moment = np.isfinite(proportions).all(axis=1) & (totals.min(axis=1) > 1)
        numerator = np.sum(np.where(valid_moment[:, None], (proportions - mean_proportion[:, None]) ** 2 - mean_proportion[:, None] * (1.0 - mean_proportion[:, None]) / np.maximum(totals, 1.0), 0.0))
        denominator = np.sum(np.where(valid_moment[:, None], mean_proportion[:, None] * (1.0 - mean_proportion[:, None]) * (1.0 - 1.0 / np.maximum(totals, 1.0)), 0.0))
        raw_dispersion = max(0.0, float(numerator / denominator)) if denominator > 0 else 0.05
        support = float(np.sum(totals[valid_moment])) if valid_moment.any() else 0.0
        dispersion = float((raw_dispersion * support + 0.05 * 200.0) / (support + 200.0))
    offsets = np.array([0.5 * offset, -0.5 * offset])
    corrected = logits - offsets[None, :]
    corrected_weights = totals / (1.0 + (totals - 1.0) * dispersion)
    finite = np.isfinite(corrected)
    numerator = np.sum(np.where(finite, corrected * corrected_weights, 0.0), axis=1)
    denominator = np.sum(np.where(finite, corrected_weights, 0.0), axis=1)
    shared_logit = np.divide(numerator, denominator, out=np.full(len(included), np.nan), where=denominator > 0)
    result = expit(shared_logit), offset
    return (*result, dispersion) if return_dispersion else result


def main():
    args = parse_args()
    catalog = load_catalog(args.event_catalog)
    if args.max_events is not None:
        catalog = catalog.head(args.max_events).copy()
    sample_levels, catalog, counts_by_primer = prepare_event_counts(args, catalog)
    contrast_sets = []
    for path in args.contrasts:
        path_text = str(path)
        if "fold0" in path_text:
            fold = 0
        elif "fold1" in path_text:
            fold = 1
        else:
            fold = "full_data"
        contrast_sets.extend((fold, record) for record in json.loads(path.read_text()))
    row_lookup = {name: index for index, name in enumerate(sample_levels)}
    fold_sample_indices = {}
    for fold, contrast in contrast_sets:
        indices = fold_sample_indices.setdefault(fold, set())
        indices.update(row_lookup[name] for name in contrast["samples_a"] if name in row_lookup)
        indices.update(row_lookup[name] for name in contrast["samples_b"] if name in row_lookup)
    psi_by_fold = {}
    offsets_by_fold = {}
    dispersions_by_fold = {}
    for fold, sample_indices in fold_sample_indices.items():
        selected = np.asarray(sorted(sample_indices), dtype=int)
        fold_psi = []
        fold_offsets = []
        fold_dispersions = []
        for event_index in range(len(catalog)):
            included = np.column_stack([counts_by_primer[primer][0][selected, event_index] for primer in range(2)])
            excluded = np.column_stack([counts_by_primer[primer][1][selected, event_index] for primer in range(2)])
            _, offset, dispersion = fit_primer_aware_psi(included, excluded, return_dispersion=True)
            all_included = np.column_stack([counts_by_primer[primer][0][:, event_index] for primer in range(2)])
            all_excluded = np.column_stack([counts_by_primer[primer][1][:, event_index] for primer in range(2)])
            psi, _, _ = fit_primer_aware_psi(all_included, all_excluded, offset=offset, dispersion=dispersion, return_dispersion=True)
            fold_psi.append(psi)
            fold_offsets.append(offset)
            fold_dispersions.append(dispersion)
        psi_by_fold[fold] = np.asarray(fold_psi, dtype=float)
        offsets_by_fold[fold] = fold_offsets
        dispersions_by_fold[fold] = fold_dispersions
    catalog = catalog.copy()
    for fold, offsets in offsets_by_fold.items():
        catalog[f"primer_offset_logit_{fold}"] = offsets
        catalog[f"primer_dispersion_{fold}"] = dispersions_by_fold[fold]
    catalog["n_samples"] = np.isfinite(psi_by_fold[next(iter(psi_by_fold))]).sum(axis=1)
    args.event_output.parent.mkdir(parents=True, exist_ok=True)
    catalog.to_csv(args.event_output, sep="\t", index=False, compression="gzip")
    output_rows = []
    for fold, contrast in contrast_sets:
        pairs = [(row_lookup.get(a), row_lookup.get(b)) for a, b in zip(contrast["samples_a"], contrast["samples_b"])]
        pairs = [(a, b) for a, b in pairs if a is not None and b is not None]
        if len(pairs) < args.minimum_pairs:
            continue
        psi = psi_by_fold[fold]
        first = psi[:, [a for a, _ in pairs]]
        second = psi[:, [b for _, b in pairs]]
        valid = np.isfinite(first) & np.isfinite(second)
        enough = valid.sum(axis=1) >= args.minimum_pairs
        differences = second - first
        p_values = fast_paired_wilcoxon(differences, valid)
        p_values[~enough] = np.nan
        q_values = bh(p_values)
        effects = np.divide(np.nansum(np.where(valid, differences, np.nan), axis=1), valid.sum(axis=1), out=np.full(len(catalog), np.nan), where=valid.sum(axis=1) > 0)
        for event_index in np.flatnonzero(enough):
            event = catalog.iloc[event_index]
            output_rows.append({
                "method": "SUPPA2 (primer aware)",
                "contrast_id": contrast["contrast_id"],
                "effect": contrast.get("effect", "cell_type"),
                "stratum": contrast.get("stratum", "all"),
                "level_a": contrast["level_a"],
                "level_b": contrast["level_b"],
                "feature_id": event.feature_id,
                "event_type": event.event_type,
                "event_id": event.event_id,
                "p_value": float(p_values[event_index]),
                "q_value": float(q_values[event_index]),
                "effect_size": float(effects[event_index]),
                "n_subjects": int(valid[event_index].sum()),
                "fold": fold,
                "gene_id": event.gene_id,
                "gene_name": event.gene_name,
                "significant": bool(q_values[event_index] < 0.05),
                "criterion": "SUPPA2 primer-aware shared-logit paired Wilcoxon, BH q < 0.05",
            })
    result = pd.DataFrame(output_rows)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    result.to_csv(args.output, sep="\t", index=False, compression="gzip")
    print(f"wrote {len(result):,} primer-aware SUPPA2 tests and {len(catalog):,} events")


if __name__ == "__main__":
    main()
