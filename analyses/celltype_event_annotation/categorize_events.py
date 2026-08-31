"""Categorize splice structures in a robust calibrated discovery set."""

from __future__ import annotations

import argparse
import gzip
import json
import pickle
import re
from pathlib import Path

import pandas as pd

from tealeaf.sc.differential import classify_splice_block, splice_block_event_components


ATTRIBUTE = re.compile(r'(\S+)\s+"([^"]+)"')


def gene_names_from_gtf(path):
    names = {}
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt") as handle:
        for line in handle:
            if not line or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 9:
                continue
            attributes = dict(ATTRIBUTE.findall(fields[8]))
            gene_id = attributes.get("gene_id")
            if gene_id is not None:
                names[gene_id] = attributes.get("gene_name", gene_id)
    return names


def load_json(path):
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt") as handle:
        return json.load(handle)


def robust_catalog(significant_paths, primary_index, candidate_cache, blocks_path, gtf_path):
    tables = [pd.read_csv(path, sep="\t") for path in significant_paths]
    robust_ids = set(tables[0]["test_id"])
    for table in tables[1:]:
        robust_ids.intersection_update(table["test_id"])
    primary = tables[primary_index].set_index("test_id")
    candidates = pickle.loads(Path(candidate_cache).read_bytes())["candidates"]
    candidates = {candidate[0]: candidate for candidate in candidates}
    blocks = load_json(blocks_path)
    blocks = {(block["gene_id"], block["block_id"]): block for block in blocks}
    gene_names = gene_names_from_gtf(gtf_path)
    records = []
    for test_id in sorted(robust_ids):
        result = primary.loc[test_id]
        candidate = candidates[test_id]
        gene_id = candidate[2]
        block = blocks[(gene_id, candidate[1])]
        signatures = candidate[6]
        components = splice_block_event_components(block, signatures)
        records.append(
            {
                "test_id": test_id,
                "block_id": candidate[1],
                "gene_id": gene_id,
                "gene_name": gene_names.get(gene_id, gene_id),
                "event_type": classify_splice_block(block, signatures),
                "event_components": "; ".join(components),
                "chromosome": block["chromosome"],
                "strand": block["strand"],
                "left_anchor": json.dumps(block["left_anchor"]),
                "right_anchor": json.dumps(block["right_anchor"]),
                "n_tested_paths": len(signatures),
                "tested_path_signatures": json.dumps(signatures),
                "n_test_levels": int(
                    result.get("n_test_levels", len(candidate[9]))
                ),
                "statistic": result["statistic"],
                "empirical_p_value": result["p_value"],
                "fdr": result["fdr"],
            }
        )
    return pd.DataFrame.from_records(records)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--significant", type=Path, nargs="+", required=True)
    parser.add_argument("--primary-index", type=int, default=2)
    parser.add_argument("--candidate-cache", type=Path, required=True)
    parser.add_argument("--blocks", type=Path, required=True)
    parser.add_argument("--gtf", type=Path, required=True)
    parser.add_argument(
        "--literature-audit",
        type=Path,
        help="Optional literature-audit table to intersect with the discovery set.",
    )
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    if not 0 <= args.primary_index < len(args.significant):
        parser.error("--primary-index must select one --significant table")
    catalog = robust_catalog(args.significant, args.primary_index, args.candidate_cache, args.blocks, args.gtf)
    args.output_dir.mkdir(parents=True, exist_ok=True)
    catalog.to_csv(args.output_dir / "event_catalog.tsv", sep="\t", index=False)
    summary = catalog.groupby("event_type", as_index=False).size().rename(columns={"size": "n_events"})
    summary["fraction"] = summary["n_events"] / len(catalog)
    summary.sort_values(["n_events", "event_type"], ascending=[False, True]).to_csv(args.output_dir / "event_type_summary.tsv", sep="\t", index=False)
    unique_blocks = catalog.drop_duplicates("block_id")
    block_summary = unique_blocks.groupby("event_type", as_index=False).size().rename(columns={"size": "n_blocks"})
    block_summary["fraction"] = block_summary["n_blocks"] / len(unique_blocks)
    block_summary.sort_values(["n_blocks", "event_type"], ascending=[False, True]).to_csv(args.output_dir / "event_type_block_summary.tsv", sep="\t", index=False)
    if args.literature_audit is not None:
        literature = pd.read_csv(args.literature_audit, sep="\t")
        evidence = catalog[
            [
                "test_id",
                "n_tested_paths",
                "n_test_levels",
                "statistic",
                "empirical_p_value",
                "fdr",
            ]
        ].rename(columns={"n_test_levels": "n_cell_types"})
        literature.merge(evidence, on="test_id", how="inner").to_csv(
            args.output_dir / "literature_matches.tsv", sep="\t", index=False
        )


if __name__ == "__main__":
    main()
