#!/usr/bin/env python3
"""Split and long-read assessments of run_local_path_tests.py fits restricted by block structure.

Blocks are classified from the annotation alone (tealeaf.sc.local_path_reads.
block_structure on all annotated mature paths of the block): "internal" blocks
have both constitutive anchors, "terminal" blocks have alternative first or
last exons. The fits, p-values and reported usages are unchanged; tests outside
the retained structures are removed from the family before the split (gene
pairs without a retained test have p = 1) and long-read assessments. With
--terminal-model, terminal-block tests come from that model instead (fitted
with run_local_path_tests.py --structure terminal, e.g. on junction-only
classes) and are added to the internal-block tests of --model; the merged
family keeps each model's own p-values.
"""

import argparse
import gzip
import json
from pathlib import Path

import pandas as pd

from extra_scripts.assess_paired_inference_audit import fit_table, split_assessment, long_read_assessment
from extra_scripts.assess_tilgner_long_read_replication import load_path_usage
from tealeaf.sc.junction_benchmark import benjamini_hochberg
from tealeaf.sc.local_path_reads import block_structure

FILTERS = {"internal": lambda s: s["structure"] == "internal", "distinct_junctions": lambda s: s["structure"] == "internal" or not s["same_junctions"]}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cache", type=Path, required=True)
    parser.add_argument("--model", required=True)
    parser.add_argument("--local-reads", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--matrix-dir", type=Path, required=True)
    parser.add_argument("--gtf", type=Path, required=True)
    parser.add_argument("--block-cache", type=Path, required=True)
    parser.add_argument("--filter", choices=sorted(FILTERS), required=True)
    parser.add_argument("--terminal-model", help="model whose tests replace the terminal blocks (use with --filter internal)")
    parser.add_argument("--long-read-gtf", type=Path, help="annotation for long-read transcript names and exon chains (e.g. comprehensive GENCODE); default --gtf")
    parser.add_argument("--chain-mapping", action="store_true", help="assign long-read transcripts to tested paths by their exon chains within the block (--local-reads blocks), instead of through the annotated splice-block cache")
    parser.add_argument("--skip-split", action="store_true", help="long-read assessment only")
    parser.add_argument("--split-universe", type=Path, help="comparator-defined split universes (extra_scripts/build_split_universe.py); default: the archived matched universes")
    args = parser.parse_args()
    repo = Path(__file__).resolve().parents[1]
    blocks = json.load(gzip.open(args.local_reads / "blocks.json.gz", "rt"))
    structure = {block["block_id"]: block_structure(block["paths"][:block.get("n_annotated", len(block["signatures"]))]) for block in blocks.values()}
    keep = {block for block, value in structure.items() if FILTERS[args.filter](value)}
    output = args.output_dir
    output.mkdir(parents=True, exist_ok=True)
    name = f"{args.model}, {args.filter} blocks" + (f" + {args.terminal_model} terminal blocks" if args.terminal_model else "")
    sources = [(args.model, keep)] + ([(args.terminal_model, set(structure) - keep)] if args.terminal_model else [])
    folds = [pd.concat([fit_table(args.cache / f"{model}_{fold}/merged/paired_path.tsv").pipe(lambda table: table.loc[table.block_id.isin(selected)]) for model, selected in sources], ignore_index=True) for fold in (0, 1)]
    if not args.skip_split:
        split_assessment(folds, repo, output, name, args.split_universe)
    full = pd.concat([pd.read_csv(args.cache / f"{model}_full/merged/paired_path.tsv", sep="\t").pipe(lambda table: table.loc[table.block_id.isin(selected)]) for model, selected in sources], ignore_index=True)
    full["fdr"] = benjamini_hochberg(full.p_value.to_numpy())
    usage = pd.concat([load_path_usage(args.cache / f"{model}_full", set(full.test_id[full.block_id.isin(selected)])) for model, selected in sources], ignore_index=True)
    usage = usage.set_index(["test_id", "cell_type", "path_number"]).proportion
    native_blocks = transcript_chains = None
    if args.chain_mapping:
        from tealeaf.sc.differential import read_gtf_exons
        from extra_scripts.assess_tilgner_long_read_replication import stable_identifier
        native_blocks = {block["block_id"]: block for block in blocks.values()}
        genes = {gene.split(".")[0] for gene in full.gene_id}
        transcript_chains = {stable_identifier(transcript): sorted(tuple(exon) for exon in chain) for gene, record in read_gtf_exons(args.long_read_gtf or args.gtf).items() if gene.split(".")[0] in genes for transcript, chain in record["transcripts"].items()}
    long_read_assessment(full, args.cache / f"{args.model}_full", args.matrix_dir, args.long_read_gtf or args.gtf, args.block_cache, output, name, usage=usage, native_blocks=native_blocks, transcript_chains=transcript_chains)
    counts = pd.Series([value["structure"] + ("" if value["structure"] == "internal" else (", same junctions" if value["same_junctions"] else ", distinct junctions")) for value in structure.values()]).value_counts().to_dict()
    (output / "manifest.json").write_text(json.dumps({"model": args.model, "terminal_model": args.terminal_model, "filter": args.filter, "blocks_by_structure": counts, "retained_blocks": len(keep), "selection": "annotation-only block structure; fits and p-values unchanged"}, indent=2) + "\n")


if __name__ == "__main__":
    main()
