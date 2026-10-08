#!/usr/bin/env python3
"""Frozen split-tail diagnostic of exonic/junction support in full STAR BAMs."""

import argparse
from dataclasses import asdict
import hashlib
import json
from pathlib import Path
import re

import numpy as np
import pandas as pd

from tealeaf.sc.differential import read_gtf_exons
from tealeaf.sc.local_read_support import LocalReadContrast, event_local_read_contrast, collect_indexed_local_read_support
from tealeaf.sc.sashimi import read_primer_cell_groups


def file_hash(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def event_envelope(event_id):
    """Published SUPPA coordinates, one-based endpoints to half-open span."""
    split_id = str(event_id).split(";", 1)
    if len(split_id) != 2:
        raise ValueError("SUPPA genomic event definition required")
    definition = split_id[1].split(":")
    if len(definition) < 5 or definition[0] not in ("SE", "RI", "MX", "A3", "A5", "AF", "AL") or definition[-1] not in ("+", "-"):
        raise ValueError("SUPPA genomic event definition required")
    coordinates = [int(value) for field in definition[2:-1] for value in field.split("-") if re.fullmatch(r"[0-9]+", value)]
    if len(coordinates) < 2 or min(coordinates) < 1 or min(coordinates) == max(coordinates):
        raise ValueError("nonempty positive SUPPA coordinate envelope required")
    return definition[1], definition[-1], (min(coordinates) - 1, max(coordinates))


def prepare(args):
    if args.output_dir.exists():
        raise ValueError("preserve earlier diagnostics, use a new directory")
    tables = []
    for fold in (0, 1):
        table = pd.read_csv(args.diagnostic_root / f"fold{fold}" / "diagnostics.tsv.gz", sep="\t")
        table["fold"] = fold
        tables.append(table)
    selected = pd.concat(tables, ignore_index=True)
    if selected.duplicated(["fold", "test_id"]).any():
        raise ValueError("unique frozen diagnostic requests required")
    selected["feature_id"] = selected.test_id.str.split("|", regex=False).str[0]
    catalog = pd.read_csv(args.event_catalog, sep="\t").set_index("feature_id", verify_integrity=True)
    if not set(selected.feature_id) <= set(catalog.index):
        raise ValueError("every frozen request must have its source event definition")
    source = catalog.loc[sorted(set(selected.feature_id))].copy()
    genes = read_gtf_exons(args.gtf)
    features, definitions = [], []
    for key, row in source.iterrows():
        feature = None
        error = ""
        try:
            chromosome, strand, span = event_envelope(row.event_id)
            record = genes[row.gene_id]
            if record["chromosome"] != chromosome or record["strand"] != strand:
                raise ValueError("event and gene annotation disagree")
            feature = event_local_read_contrast(key, row.gene_id, record, row.included.split(","), row.excluded.split(","), span)
            features.append(asdict(feature))
        except (KeyError, ValueError) as exc:
            error = str(exc)
        has_markers = feature is not None and any((feature.included_exons, feature.excluded_exons, feature.included_junctions, feature.excluded_junctions))
        definitions.append(dict(feature_id=key, gene_id=row.gene_id, event_type=row.event_type, status="ok" if has_markers else "no_class_unanimous_features" if feature else "unavailable", error=error, included_exon_segments=len(feature.included_exons) if feature else 0, excluded_exon_segments=len(feature.excluded_exons) if feature else 0, included_junctions=len(feature.included_junctions) if feature else 0, excluded_junctions=len(feature.excluded_junctions) if feature else 0))
    bams = sorted(args.starsolo_root.glob("*/Aligned.sortedByCoord.out.bam"))
    if len(bams) != args.expected_bams:
        raise ValueError("whole declared alignment-file family required")
    metadata = pd.read_csv(args.metadata, sep="\t", dtype=str)
    groups, _ = read_primer_cell_groups(args.metadata, args.primer_pairs, sorted(set(metadata.cell_type)))
    args.output_dir.mkdir(parents=True)
    selected.to_csv(args.output_dir / "selected_cases.tsv.gz", sep="\t", index=False)
    pd.DataFrame(definitions).to_csv(args.output_dir / "feature_definitions.tsv", sep="\t", index=False)
    diagnostic_inputs = [args.diagnostic_root / f"fold{fold}" / name for fold in (0, 1) for name in ("diagnostics.tsv.gz", "subject_influence.tsv.gz")]
    diagnostic_inputs.extend(args.bound_root / f"fold{fold}" / "diagnostics.tsv.gz" for fold in (0, 1))
    recipe = dict(events=features, barcode_groups=groups, bams=[dict(path=str(path.resolve()), size=path.stat().st_size, mtime_ns=path.stat().st_mtime_ns) for path in bams], inputs={str(path): file_hash(path) for path in (args.event_catalog, args.gtf, args.metadata, args.primer_pairs)}, diagnostic_inputs={str(path): file_hash(path) for path in diagnostic_inputs}, diagnostic_root=str(args.diagnostic_root), selected_sha256=file_hash(args.output_dir / "selected_cases.tsv.gz"), selection="identical prior strong raw-p tail and informative-subject-count-matched weak controls in both splits, no LR selection", scope="same-read local support diagnostic; unique STAR GX/NH plus per-file barcode/UB keys, not upstream EC-UMI counts, RNA PSI, a significance filter or a new inference model", production_changes=False)
    (args.output_dir / "recipe.json").write_text(json.dumps(recipe) + "\n")
    print(pd.DataFrame(definitions).groupby(["status", "event_type"]).size().to_string(), flush=True)


def collect(args):
    import pysam
    recipe_path = args.output_dir / "recipe.json"
    recipe = json.loads(recipe_path.read_text())
    if args.shard_index >= len(recipe["bams"]):
        raise ValueError("shard exceeds declared alignment family")
    packet = recipe["bams"][args.shard_index]
    bam = Path(packet["path"])
    if bam.stat().st_size != packet["size"] or bam.stat().st_mtime_ns != packet["mtime_ns"]:
        raise ValueError("alignment input changed since frozen preparation")
    shard = args.output_dir / f"shard_{args.shard_index}"
    if shard.exists():
        raise ValueError("preserve previous shard, do not overwrite")
    shard.mkdir()
    index = shard / "alignments.bai"
    pysam.index("-@", str(args.threads), "-o", str(index), str(bam))
    events = [decode_event(row) for row in recipe["events"]]
    counts, filters = collect_indexed_local_read_support(bam, index, recipe["barcode_groups"], events)
    rows = [dict(feature_id=key[0], subject=key[1], cell_type=key[2], primer=key[3], signature=key[4], count=count) for key, count in counts.items()]
    pd.DataFrame(rows, columns=["feature_id", "subject", "cell_type", "primer", "signature", "count"]).to_csv(shard / "support.tsv.gz", sep="\t", index=False)
    receipt = dict(input=packet, recipe_sha256=file_hash(recipe_path), shard=args.shard_index, requested_events=len(events), filters=dict(filters), complete=True)
    (shard / "manifest.json").write_text(json.dumps(receipt, indent=2) + "\n")
    print(json.dumps(receipt, indent=2), flush=True)


def decode_event(row):
    """Round-trip immutable feature intervals through a JSON recipe."""
    fields = {key: tuple(tuple(interval) for interval in value) if key.endswith("exons") or key.endswith("junctions") else value for key, value in row.items()}
    return LocalReadContrast(**fields)


def summarize_signature_counts(table):
    """Diagnostic barcode/UMI keys, conflicts do not count for either class."""
    result = table.copy()
    inc, exc = (result.signature.to_numpy(dtype=int) & bits != 0 for bits in (5, 10))
    local = result.signature.to_numpy(dtype=int) & 16 != 0
    for key, mask in dict(included_only=inc & ~exc, excluded_only=exc & ~inc, conflicting=inc & exc, local_nondiscriminating=local & ~inc & ~exc, distal_only=~local & ~inc & ~exc, gene_keys=np.ones(len(table), dtype=bool)).items():
        result[key] = result["count"] * mask
    return result.groupby(["feature_id", "subject", "cell_type", "primer"], sort=True)[["included_only", "excluded_only", "conflicting", "local_nondiscriminating", "distal_only", "gene_keys"]].sum().reset_index()


def collate(args):
    if args.public_dir is None or args.public_dir.exists():
        raise ValueError("new public output directory required")
    recipe_path = args.output_dir / "recipe.json"
    recipe = json.loads(recipe_path.read_text())
    selected_path = args.output_dir / "selected_cases.tsv.gz"
    if file_hash(selected_path) != recipe["selected_sha256"]:
        raise ValueError("frozen selection changed")
    if any(file_hash(path) != expected for path, expected in recipe["diagnostic_inputs"].items()):
        raise ValueError("frozen diagnostics or influence bounds changed")
    tables, receipts = [], []
    for index, packet in enumerate(recipe["bams"]):
        shard = args.output_dir / f"shard_{index}"
        receipt = json.loads((shard / "manifest.json").read_text())
        if receipt["complete"] is not True or receipt["shard"] != index or receipt["input"] != packet or receipt["recipe_sha256"] != file_hash(recipe_path) or receipt["requested_events"] != len(recipe["events"]):
            raise ValueError("every frozen alignment shard must be complete and compatible")
        table = pd.read_csv(shard / "support.tsv.gz", sep="\t")
        if table.duplicated(["feature_id", "subject", "cell_type", "primer", "signature"]).any() or not table.feature_id.isin([row["event_id"] for row in recipe["events"]]).all() or not table.signature.between(0, 31).all() or (table["count"] < 0).any():
            raise ValueError("unique nonnegative declared support signatures required")
        table["run"] = packet["library"] if "library" in packet else Path(packet["path"]).parent.name
        tables.append(table)
        receipts.append(receipt)
    signatures = pd.concat(tables, ignore_index=True)
    counts = summarize_signature_counts(signatures)
    count_lookup = counts.set_index(["feature_id", "subject", "cell_type", "primer"], verify_integrity=True).to_dict("index")
    selected = pd.read_csv(selected_path, sep="\t")
    features = pd.read_csv(args.output_dir / "feature_definitions.tsv", sep="\t").set_index("feature_id")
    groups, subjects = [], []
    for fold in (0, 1):
        evidence = pd.read_csv(Path(recipe["diagnostic_root"]) / f"fold{fold}" / "subject_influence.tsv.gz", sep="\t")
        bounds = pd.read_csv(args.bound_root / f"fold{fold}" / "diagnostics.tsv.gz", sep="\t").set_index("test_id", verify_integrity=True)
        for row in selected.loc[selected.fold.eq(fold)].itertuples(index=False):
            local = evidence.loc[evidence.test_id.eq(row.test_id)].copy()
            if local.empty or not set(local.panel) == {row.panel} or row.test_id not in bounds.index:
                raise ValueError("same original subject diagnostic and bound family required")
            levels = row.test_id.split("|")[-2:]
            count_columns = ["included_only", "excluded_only", "conflicting", "local_nondiscriminating", "distal_only", "gene_keys"]
            case = []
            for subject in local.subject:
                record = dict(fold=fold, test_id=row.test_id, feature_id=row.feature_id, subject=subject)
                for level_index, level in enumerate(levels):
                    for primer in ("poly(dT)", "random hexamer"):
                        match = count_lookup.get((row.feature_id, subject, level, primer), {})
                        for column in count_columns:
                            record[f"{column}_{level_index}_{'DT' if primer == 'poly(dT)' else 'RH'}"] = int(match.get(column, 0))
                    inc = record[f"included_only_{level_index}_DT"] + record[f"included_only_{level_index}_RH"]
                    exc = record[f"excluded_only_{level_index}_DT"] + record[f"excluded_only_{level_index}_RH"]
                    record[f"discriminating_keys_{level_index}"] = inc + exc
                    record[f"inclusion_fraction_{level_index}"] = inc / (inc + exc) if inc + exc else np.nan
                record["local_inclusion_difference"] = record["inclusion_fraction_1"] - record["inclusion_fraction_0"]
                case.append(record)
            case = pd.DataFrame(case).merge(local[["subject", "retained", "precision_share", "pseudo_effect"]], on="subject", how="left", validate="one_to_one")
            usable = case.retained.astype(str).str.lower().eq("true")
            measured = case.loc[usable]
            good = measured.discriminating_keys_0.ge(5) & measured.discriminating_keys_1.ge(5)
            pair = measured.local_inclusion_difference.dropna()
            dominant = measured.loc[measured.subject.eq(row.dominant_subject)]
            if len(dominant) != 1 or len(measured) != row.n_informative_subjects:
                raise ValueError("retained subject and dominant identity changed")
            groups.append(dict(fold=fold, test_id=row.test_id, panel=row.panel, feature_id=row.feature_id, event_type=features.at[row.feature_id, "event_type"], feature_status=features.at[row.feature_id, "status"], n_informative_subjects=len(measured), n_subjects_with_five_discriminating_keys_each_type=int(good.sum()), n_subjects_with_any_discriminating_keys_each_type=len(pair), n_retained_with_no_local_discriminating_keys=int((measured.discriminating_keys_0 + measured.discriminating_keys_1).eq(0).sum()), unavoidable_share_gt_90=bool(bounds.at[row.test_id, "largest_uniform_lower_bound"] > .9), original_maximum_precision_share=row.maximum_precision_share, dominant_minimum_discriminating_keys=int(dominant[["discriminating_keys_0", "discriminating_keys_1"]].min(axis=1).iloc[0]), mean_local_inclusion_difference=float(pair.mean()) if len(pair) else np.nan, raw_score_direction_agrees_with_local_mean=bool(pair.mean() * row.fitted_mean > 0) if len(pair) else np.nan))
            subjects.append(case)
    diagnostics = pd.DataFrame(groups)
    summaries = []
    for (fold, panel, unavoidable), local in diagnostics.groupby(["fold", "panel", "unavoidable_share_gt_90"], dropna=False):
        summaries.append(dict(fold=fold, panel=panel, unavoidable_share_gt_90=bool(unavoidable), requested=len(local), available_features=int(local.feature_status.eq("ok").sum()), at_least_four_subjects_with_five_discriminating_keys_each_type=int(local.n_subjects_with_five_discriminating_keys_each_type.ge(4).sum()), median_subjects_with_five_discriminating_keys_each_type=local.n_subjects_with_five_discriminating_keys_each_type.median(), median_dominant_minimum_discriminating_keys=local.dominant_minimum_discriminating_keys.median(), local_direction_available=int(local.raw_score_direction_agrees_with_local_mean.notna().sum()), raw_score_local_direction_agreement=local.raw_score_direction_agrees_with_local_mean.dropna().astype(float).mean()))
    args.public_dir.mkdir(parents=True)
    diagnostics.to_csv(args.public_dir / "diagnostics.tsv.gz", sep="\t", index=False)
    pd.concat(subjects, ignore_index=True).to_csv(args.public_dir / "subject_support.tsv.gz", sep="\t", index=False)
    signatures.to_csv(args.public_dir / "run_support_signatures.tsv.gz", sep="\t", index=False)
    features.to_csv(args.public_dir / "feature_definitions.tsv", sep="\t")
    pd.DataFrame(summaries).to_csv(args.public_dir / "summary.tsv", sep="\t", index=False)
    deduplication = "Exact UB keys are unioned within each physical library, never across libraries; this does not reproduce joint upstream UMI error correction." if recipe.get("library_union") else "Keys may duplicate across alignment files and differ from original EC UMI deduplication."
    manifest = dict(recipe_sha256=file_hash(recipe_path), shards=receipts, requested_tests=len(selected), diagnostics=len(diagnostics), selection=recipe["selection"], scope=recipe["scope"], caveats="Published coordinate envelopes and class-unanimous annotation markers do not define an exact read-path likelihood. One-base exon overlaps count. Unmodeled isoforms may contain markers. " + deduplication + " Pooled descriptive fractions combine primer-specific opportunities and are NOT RNA PSI or independent validation. Threshold five summarizes read support only, no events or subjects are removed from testing.", production_changes=False)
    (args.public_dir / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    print(pd.DataFrame(summaries).to_string(index=False), flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("mode", choices=("prepare", "collect", "collate"))
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--diagnostic-root", type=Path)
    parser.add_argument("--event-catalog", type=Path)
    parser.add_argument("--gtf", type=Path)
    parser.add_argument("--metadata", type=Path)
    parser.add_argument("--primer-pairs", type=Path)
    parser.add_argument("--starsolo-root", type=Path)
    parser.add_argument("--expected-bams", type=int, default=16)
    parser.add_argument("--shard-index", type=int)
    parser.add_argument("--threads", type=int, default=4)
    parser.add_argument("--public-dir", type=Path)
    parser.add_argument("--bound-root", type=Path, default=Path("analyses/event_sequence_variance_limits"))
    args = parser.parse_args()
    if args.mode == "prepare":
        if any(value is None for value in (args.diagnostic_root, args.event_catalog, args.gtf, args.metadata, args.primer_pairs, args.starsolo_root)):
            parser.error("all source inputs required for preparation")
        prepare(args)
    elif args.mode == "collect":
        if args.shard_index is None or args.shard_index < 0 or args.threads < 1:
            parser.error("nonnegative shard and positive threads required")
        collect(args)
    else:
        collate(args)


if __name__ == "__main__":
    main()
