"""Long-read direction agreement of full-data tests against barcode-tagged long-read molecules (e.g. Kinnex Parse PBMC).

Each molecule's exon chain is assigned to the tested path it fits (chain mapping, as with --chain-mapping in
assess_block_structure_filter.py), and molecules are pooled by long-read cell type. Cells take their type from a label
table (e.g. CellTypist) in the short-read vocabulary; replicates are two halves of the cells by RT well (the assessment
pools them). Writes lr_mapping.tsv.gz, lr_rank_summary.tsv and lr_rank.tsv.gz to --output-dir.
"""
import argparse
import gzip
import json
from pathlib import Path

import pandas as pd

from extra_scripts.assess_paired_inference_audit import long_read_assessment
from extra_scripts.assess_tilgner_long_read_replication import load_path_usage
from tealeaf.sc.junction_benchmark import benjamini_hochberg
from tealeaf.sc.longread_molecules import cell_labels, molecule_source


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--cache", type=Path, required=True, help="fit cache (<model>_full/merged/paired_path.tsv)")
    parser.add_argument("--models", required=True, help="comma-separated models whose full-data tests are pooled (e.g. internal and terminal fits)")
    parser.add_argument("--blocks", type=Path, required=True, help="blocks.json.gz of the fitted collection")
    parser.add_argument("--molecules", type=Path, required=True, help="molecules.parquet from longread_molecules.py collapse")
    parser.add_argument("--labels", type=Path, required=True, help="cell -> label table (tsv)")
    parser.add_argument("--label-column", default="Immune_All_Low_majority")
    parser.add_argument("--merge", nargs="*", default=[], help="label renames old=new after spaces become underscores")
    parser.add_argument("--suffix", default="", help="appended to labels to match level names (e.g. |Baseline)")
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    models = args.models.split(",")
    full = pd.concat([pd.read_csv(args.cache / f"{model}_full/merged/paired_path.tsv", sep="\t") for model in models], ignore_index=True)
    full["fdr"] = benjamini_hochberg(full.p_value.to_numpy())
    usage = pd.concat([load_path_usage(args.cache / f"{model}_full", set(full.test_id)) for model in models], ignore_index=True).set_index(["test_id", "cell_type", "path_number"]).proportion
    blocks = json.load(gzip.open(args.blocks, "rt"))
    native_blocks = {block["block_id"]: block for block in blocks.values()}
    labels = cell_labels(args.labels, args.label_column, args.merge, args.suffix)
    molecules = pd.read_parquet(args.molecules, columns=["cell", "gene", "chain"])
    matrix, features, columns, chains = molecule_source(molecules, labels)
    cells = labels.loc[labels.index.isin(set(molecules.cell))].groupby("cell_type").size()
    print(f"{len(molecules):,} molecules, {len(features):,} distinct chains; cells per type:\n{cells.to_string()}", flush=True)
    long_read_assessment(full, args.cache / f"{models[0]}_full", None, None, None, args.output_dir, ",".join(models), usage=usage, native_blocks=native_blocks, transcript_chains=chains, source=(matrix, features, columns))
    (args.output_dir / "manifest.json").write_text(json.dumps({"models": models, "molecules": str(args.molecules), "labels": str(args.labels), "label_column": args.label_column, "merge": args.merge, "suffix": args.suffix, "cells_per_type": cells.to_dict()}, indent=2) + "\n")


if __name__ == "__main__":
    main()
