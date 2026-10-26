#!/usr/bin/env python3
"""Full-data discoveries per method, by event type.

Junction-cluster methods (LeafCutter, scQuint, MAJIQ Heterogen) are read from
--junction-table, either all-tested long-read source mappings or a comparison
all_tests table (cell-type effects are kept), and MAJIQ also from
--majiq-tests (edges collapsed to LSVs by their smallest p-value); LeafCutter clusters without a
gene take the unique gene of their junctions (--leafcutter-counts with
--junction-annotation). With --contrasts common only the cell-type pairs
tested by every method are kept; with all, each method keeps every pair it
tests, and methods with fewer pairs than the most are labelled with their
number of pairs. Every method's finite p-values form one BH family per method,
and discoveries are tests with BH-adjusted p at most --fdr. Two counts are
drawn: tests (each method's own unit, a block test, intron cluster, intron
group, LSV or event) and gene--cell-type-pair keys with at least one
discovery, the key taking the event type of its smallest-p discovery. Tealeaf
event classes come from categorize_block_tests.py; SUPPA2 and rMATS event
types map onto the same classes; junction-cluster methods have no event type.
Needs plotnine.
"""

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
import plotnine as p9

from tealeaf.sc.ds_benchmark import leafcutter_cluster_gene_map
from tealeaf.sc.junction_benchmark import benjamini_hochberg

JUNCTION = {"LeafCutter": "LeafCutter", "scQuint": "scQuint", "MAJIQ Heterogen": "MAJIQ"}
COMPARATOR_EVENTS = {"SE": "cassette exon", "A5SS": "alternative donor", "A3SS": "alternative acceptor", "MXE": "mutually exclusive exons", "RI": "retained intron"}
# categorical slots in fixed order; the classes absent from comparators fold into "other"
EVENTS = {"cassette exon": "#2a78d6", "alternative donor": "#eb6834", "alternative acceptor": "#1baf7a", "mutually exclusive exons": "#eda100", "retained intron": "#e87ba4", "alternative first exon": "#008300", "alternative last exon": "#6250d6", "other": "#e34948", "no event type": "#b5b3ad"}


def pair_key(a, b):
    return ["||".join(sorted(pair)) for pair in zip(a, b)]


def standard(table, method, event):
    table = table.assign(method=method, pair=pair_key(table.level_a, table.level_b), event=event, gene=table.gene_id.astype("string").str.split(".").str[0])
    table["p_value"] = pd.to_numeric(table.p_value, errors="coerce")
    return table.loc[np.isfinite(table.p_value), ["method", "pair", "gene", "event", "p_value"]]


def read_junction_table(path):
    table = pd.read_csv(path, sep="\t", usecols=lambda column: column in ("method", "effect", "level_a", "level_b", "feature_id", "gene_id", "p_value"), low_memory=False)
    return table.loc[table.effect.eq("cell_type")] if "effect" in table else table


def read_majiq_tests(directory):
    """MAJIQ Heterogen normalized edge tests (one file per contrast) collapsed to LSVs, each with its smallest edge p-value."""
    edges = pd.concat([pd.read_csv(path, sep="\t", usecols=["method", "effect", "level_a", "level_b", "feature_id", "gene_id", "p_value"]) for path in sorted(Path(directory).glob("*.tsv"))], ignore_index=True)
    edges = edges.loc[edges.effect.eq("cell_type")]
    edges["feature_id"] = edges.feature_id.str.split(":").str[:6].str.join(":")
    return edges.groupby(["method", "effect", "level_a", "level_b", "feature_id", "gene_id"], as_index=False).p_value.min()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--tealeaf-events", type=Path, required=True, help="test_events.tsv.gz of categorize_block_tests.py")
    parser.add_argument("--junction-table", type=Path, action="append", default=[])
    parser.add_argument("--majiq-tests", type=Path, help="directory of MAJIQ Heterogen normalized per-contrast edge tests, collapsed to LSVs; replaces MAJIQ rows of --junction-table")
    parser.add_argument("--suppa2", type=Path, required=True, help="SUPPA2 full-data tests")
    parser.add_argument("--suppa2-primer-aware", type=Path, required=True)
    parser.add_argument("--rmats", type=Path, required=True, help="paired rMATS full-data tests")
    parser.add_argument("--leafcutter-counts", type=Path, help="LeafCutter perind numerator counts, to map clusters without a gene")
    parser.add_argument("--junction-annotation", type=Path, help="junction table with chromosome, start, end and gene_id")
    parser.add_argument("--contrasts", choices=("common", "all"), default="common")
    parser.add_argument("--exclude-method", action="append", default=[])
    parser.add_argument("--fdr", type=float, default=0.05)
    parser.add_argument("--output", type=Path, required=True, help="figure; the counts are written next to it as .tsv")
    args = parser.parse_args()
    junctions = [read_junction_table(path) for path in args.junction_table]
    if args.majiq_tests:
        junctions = [table.loc[table.method.ne("MAJIQ Heterogen")] for table in junctions] + [read_majiq_tests(args.majiq_tests)]
    junctions = pd.concat(junctions, ignore_index=True).astype({"gene_id": "object"})
    leaf = junctions.method.eq("LeafCutter") & junctions.gene_id.isna()
    if leaf.any() and args.leafcutter_counts:
        genes = leafcutter_cluster_gene_map(args.junction_annotation, args.leafcutter_counts).set_index("feature_id").gene_id
        junctions.loc[leaf, "gene_id"] = junctions.loc[leaf, "feature_id"].map(genes)
    # scQuint intron groups are named gene_chromosome_position_strand
    named = junctions.method.eq("scQuint") & junctions.gene_id.isna()
    junctions.loc[named, "gene_id"] = junctions.loc[named, "feature_id"].str.rsplit("_", n=3).str[0]
    tables = [standard(local, JUNCTION[method], "no event type") for method, local in junctions.groupby("method")]
    tealeaf = pd.read_csv(args.tealeaf_events, sep="\t", usecols=["level_a", "level_b", "gene_id", "event", "p_value"])
    tealeaf = tealeaf.loc[tealeaf.event.ne("untested")]
    tables.append(standard(tealeaf, "Tealeaf", tealeaf.event.where(tealeaf.event.isin(EVENTS), "other")))
    columns = ["level_a", "level_b", "gene_id", "event_type", "p_value"]
    for method, path in (("SUPPA2", args.suppa2), ("SUPPA2 primer-aware", args.suppa2_primer_aware), ("rMATS", args.rmats)):
        local = pd.concat([chunk.loc[chunk.event_type.notna()] for chunk in pd.read_csv(path, sep="\t", usecols=columns, chunksize=500_000)], ignore_index=True)
        tables.append(standard(local, method, local.event_type.map(COMPARATOR_EVENTS).fillna("other")))
    tests = pd.concat(tables, ignore_index=True)
    tests = tests.loc[~tests.method.isin(args.exclude_method)]
    if args.contrasts == "common":
        tests = tests.loc[tests.pair.isin(set.intersection(*(set(local.pair) for _, local in tests.groupby("method"))))]
    tests["fdr"] = tests.groupby("method").p_value.transform(lambda p: benjamini_hochberg(p.to_numpy()))
    found = tests.loc[tests.fdr.le(args.fdr)]
    keys = found.loc[found.gene.notna()].sort_values("p_value").drop_duplicates(["method", "pair", "gene"])
    print("discoveries without a unique gene:", found.loc[found.gene.isna()].groupby("method").size().to_dict(), flush=True)
    counts = pd.concat([found.groupby(["method", "event"]).size().rename("count").reset_index().assign(unit="Tests"), keys.groupby(["method", "event"]).size().rename("count").reset_index().assign(unit="Gene–cell-type pairs")], ignore_index=True)
    totals = counts.groupby(["unit", "method"])["count"].sum().rename("total").reset_index()
    tested = tests.groupby("method").agg(tests=("p_value", "size"), contrasts=("pair", "nunique")).reset_index()
    print(f"{tests.pair.nunique()} contrasts", flush=True)
    print(tested.merge(totals.pivot(index="method", columns="unit", values="total").reset_index(), on="method").to_string(index=False), flush=True)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    counts.merge(tested, on="method").to_csv(args.output.with_suffix(".tsv"), sep="\t", index=False)
    pairs = tested.set_index("method").contrasts
    labels = {method: method if n == pairs.max() else f"{method} ({n} pairs)" for method, n in pairs.items()}
    order = [labels[method] for method in totals.loc[totals.unit.eq("Gene–cell-type pairs")].sort_values("total").method]
    for frame in (counts, totals):
        frame["method"] = pd.Categorical(frame.method.map(labels), categories=order)
        frame["unit"] = pd.Categorical(frame.unit, categories=["Tests", "Gene–cell-type pairs"])
    counts["event"] = pd.Categorical(counts.event, categories=list(EVENTS)[::-1])
    plot = (p9.ggplot(counts, p9.aes("method", "count", fill="event"))
            + p9.geom_col(width=0.7, color="white", size=0.3)
            + p9.geom_text(p9.aes("method", "total", label="total"), data=totals, inherit_aes=False, ha="left", nudge_y=0, size=8, format_string=" {:,}")
            + p9.facet_wrap("~unit", scales="free_x")
            + p9.coord_flip()
            + p9.scale_fill_manual(values=EVENTS, breaks=list(EVENTS), name="")
            + p9.scale_y_continuous(expand=(0, 0, 0.18, 0))
            + p9.labs(x="", y=f"Discoveries (BH {args.fdr:g})", fill="")
            + p9.theme_bw() + p9.theme(figure_size=(9, 3.6), legend_position="bottom", panel_grid_major_y=p9.element_blank(), panel_grid_minor=p9.element_blank(), panel_spacing_x=0.04, legend_key_size=10, legend_text=p9.element_text(size=8)))
    plot.save(args.output, verbose=False)


if __name__ == "__main__":
    main()
