#!/usr/bin/env python3
"""Data panels of the Tealeaf paper figures (docs/tealeaf2_paper.tex), one PDF per panel in --output-dir.

Panels are assembled with their letters in LaTeX. Inputs are the evaluation outputs named in PANELS below (paths
relative to the repository); every method except SUPPA2 is coverage-null filtered. Needs plotnine.
"""

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
import plotnine as p9
from mizani.breaks import breaks_extended

METHODS = {"Tealeaf": "#17324D", "MAJIQ": "#168C8C", "LeafCutter": "#7FB3D5", "scQuint": "#E6813B", "rMATS": "#9C6ADE", "SUPPA2": "#8C8C8C", "SUPPA2-PA": "#C8C8C8"}
NAMES = {"MAJIQ Heterogen": "MAJIQ", "rMATS (paired)": "rMATS", "rMATS paired JCEC": "rMATS", "SUPPA2 primer-aware": "SUPPA2-PA", "SUPPA2 (primer aware)": "SUPPA2-PA", "SUPPA2 (full data)": "SUPPA2",
         "SUPPA2 (full data); exact paired Wilcoxon": "SUPPA2", "SUPPA2 (primer aware); exact paired Wilcoxon": "SUPPA2-PA"}
PRIMERS = {"poly(dT)": "#2C7FB8", "random hexamer": "#D95F0E"}
EVENTS = {"cassette exon": "#2a78d6", "alternative donor": "#eb6834", "alternative acceptor": "#1baf7a", "mutually exclusive exons": "#eda100", "retained intron": "#e87ba4", "alternative first exon": "#008300", "alternative last exon": "#6250d6", "other": "#e34948", "no event type": "#b5b3ad"}
COMPARATORS = ["LeafCutter", "MAJIQ", "scQuint", "rMATS", "SUPPA2", "SUPPA2-PA"]


def theme(width, height, legend="right"):
    return (p9.theme_classic(base_size=7, base_family="DejaVu Sans")
            + p9.theme(figure_size=(width, height), legend_position=legend, legend_key_size=7, legend_text=p9.element_text(size=6), legend_title=p9.element_blank(),
                       strip_background=p9.element_blank(), strip_text=p9.element_text(size=7), axis_line=p9.element_line(size=0.4), axis_ticks=p9.element_line(size=0.4),
                       panel_grid_major_y=p9.element_line(color="#EBEBEB", size=0.3), plot_margin=0.01))


def save(plot, directory, name):
    plot.save(directory / f"{name}.pdf", verbose=False)
    print(name, flush=True)


def read_end_profiles(path):
    profiles = json.loads(Path(path).read_text())
    width = profiles["bin_width"]
    rows = [{"primer": primer, "end": {"five_prime": "Distance from 5' end", "three_prime": "Distance from 3' end"}[end], "distance": (index + 0.5) * width, "density": value}
            for primer, ends in profiles["profiles"].items() for end, values in ends.items() for index, value in enumerate(values)]
    return pd.DataFrame(rows)


def read_end_metaplot(paths, width=25, gap=900):
    """Mean per-gene share of read ends per 100 bases near the 5' and 3' ends, laid out as one gene (5' left, 3' right).

    Shares are of all the gene's reads, including those more than 2 kb from both ends (the omitted middle)."""
    rows = []
    for path in paths:
        for primer, genes in json.loads(Path(path).read_text()).items():
            for gene in genes:
                five, three = np.asarray(gene["five"], float), np.asarray(gene["three"], float)
                total = five.sum() + three.sum() + float(gene["middle"])
                if total > 0:
                    rows.append((primer, five / total, three / total))
    out = []
    for primer in sorted({row[0] for row in rows}):
        five = np.mean([row[1] for row in rows if row[0] == primer], axis=0)
        three = np.mean([row[2] for row in rows if row[0] == primer], axis=0)
        n = len(five)
        x_five = (np.arange(n) + 0.5) * width
        x_three = n * width + gap + (n * width - (np.arange(n) + 0.5) * width)
        scale = 100 * 100 / width  # percent of the gene's reads per 100 bases
        out.append(pd.DataFrame({"primer": primer, "x": x_five, "share": five * scale, "end": "5'"}))
        out.append(pd.DataFrame({"primer": primer, "x": x_three, "share": three * scale, "end": "3'"}))
    return pd.concat(out, ignore_index=True), n * width, gap


def metaplot(table, span, gap, width, height):
    breaks = [0, 1000, 2000, span + gap, span + gap + 1000, 2 * span + gap]
    labels = ["5' end", "+1 kb", "+2 kb", "-2 kb", "-1 kb", "3' end"]
    return (p9.ggplot(table, p9.aes("x", "share", color="primer", group="primer + end"))
            + p9.annotate("rect", xmin=span, xmax=span + gap, ymin=-np.inf, ymax=np.inf, fill="#F0F0F0")
            + p9.annotate("text", x=span + gap / 2, y=table.share.max() * 0.95, label="//", size=8, color="#888888")
            + p9.geom_line(size=0.6)
            + p9.scale_color_manual(values=PRIMERS)
            + p9.scale_x_continuous(breaks=breaks, labels=labels, expand=(0.01, 0))
            + p9.labs(x="Position of read end in the transcript (spliced bases)", y="Reads per 100 bases\n(% of the gene's reads)")
            + theme(width, height, "bottom"))


def split_table(path, universe):
    table = pd.read_csv(path, sep="\t")
    table = table.loc[table.universe.eq(universe) & ~table.comparison.str.contains("hybrid")].copy()
    table["comparator"] = table.comparison.replace(NAMES)
    table["method"] = np.where(table.method.eq("Tealeaf"), "Tealeaf", "Comparator")
    return table[["comparator", "method", "shared_genes", "replicated_bh"]]


def split_plot(table, width, height):
    order = [name for name in COMPARATORS if name in set(table.comparator)]
    labels = {name: f"{name}\n{genes:,} genes" for name, genes in table.drop_duplicates("comparator")[["comparator", "shared_genes"]].itertuples(index=False)}
    table = table.assign(comparator=pd.Categorical(table.comparator.map(labels), categories=[labels[name] for name in order[::-1]]), method=pd.Categorical(table.method, categories=["Comparator", "Tealeaf"]))
    table["label_y"] = table.replicated_bh + table.replicated_bh.max() * 0.02
    fill = {"Tealeaf": METHODS["Tealeaf"], "Comparator": "#B9C7D3"}
    return (p9.ggplot(table, p9.aes("comparator", "replicated_bh", fill="method"))
            + p9.geom_col(position=p9.position_dodge(width=0.8), width=0.75)
            + p9.geom_text(p9.aes(y="label_y", label="replicated_bh"), position=p9.position_dodge(width=0.8), ha="left", size=5.5)
            + p9.coord_flip() + p9.scale_fill_manual(values=fill, breaks=["Tealeaf", "Comparator"])
            + p9.scale_y_continuous(expand=(0, 0, 0.12, 0))
            + p9.labs(x="", y="Replicated genes")
            + theme(width, height, "bottom") + p9.theme(panel_grid_major_y=p9.element_blank(), axis_text_y=p9.element_text(size=6.5)))


def split_scatter(table, width, height, offsets):
    """Tealeaf against comparator replicated genes, one point per comparator and universe.

    offsets maps a label (comparator names joined by ", " where points coincide) to its comparator and the label
    offset (dx, dy) from the matched-universe point, as fractions of the axis range.
    """
    wide = table.pivot_table(index=["comparator", "universe"], columns="method", values="replicated_bh").reset_index()
    wide["comparator"] = pd.Categorical(wide.comparator, categories=[c for c in COMPARATORS if c in set(wide.comparator)])
    top = float(max(wide.Tealeaf.max(), wide.Comparator.max())) * 1.12
    matched = wide.loc[wide.universe.eq("matched")].set_index("comparator", drop=False)
    labels = pd.DataFrame([{"label": label, "comparator": comparator, "x": matched.loc[comparator, "Comparator"], "y": matched.loc[comparator, "Tealeaf"],
                            "lx": matched.loc[comparator, "Comparator"] + dx * top, "ly": matched.loc[comparator, "Tealeaf"] + dy * top, "ha": "left" if dx >= 0 else "right"}
                           for label, (comparator, dx, dy) in offsets.items()])
    plot = (p9.ggplot(wide, p9.aes("Comparator", "Tealeaf", color="comparator", shape="universe"))
            + p9.geom_abline(intercept=0, slope=1, linetype="dashed", color="#888888", size=0.3)
            + p9.geom_segment(p9.aes(x="x", y="y", xend="lx", yend="ly"), data=labels, inherit_aes=False, color="#AAAAAA", size=0.3)
            + p9.geom_point(size=2.0, fill="white", stroke=0.7))
    for ha in ("left", "right"):
        local = labels.loc[labels.ha.eq(ha)]
        if len(local):
            plot += p9.geom_text(p9.aes("lx", "ly", label="label", color="comparator"), data=local, inherit_aes=False, ha=ha, va="center", size=6, nudge_x=(0.01 if ha == "left" else -0.01) * top, show_legend=False)
    return (plot + p9.scale_color_manual(values=METHODS, guide=None)
            + p9.scale_shape_manual(values={"matched": "o", "comparator": "^"}, labels={"matched": "Matched universe", "comparator": "Comparator universe"})
            + p9.scale_x_continuous(limits=(0, top)) + p9.scale_y_continuous(limits=(0, top))
            + p9.coord_fixed()
            + p9.labs(x="Comparator replicated genes", y="Tealeaf replicated genes")
            + theme(width, height, "bottom"))


def split_both(path):
    rows = [split_table(path, universe).assign(universe=universe) for universe in ("matched", "comparator")]
    return pd.concat(rows, ignore_index=True)


def lr_curves(tealeaf_path, comparator_path, max_rank=200):
    tealeaf = pd.read_csv(tealeaf_path, sep="\t", low_memory=False).sort_values("rank").head(max_rank)
    agree = tealeaf.pooled_replicated.astype(str).str.lower().eq("true").to_numpy(float)
    curves = [pd.DataFrame({"method": "Tealeaf", "rank": np.arange(1, len(agree) + 1), "agreement": np.cumsum(agree) / np.arange(1, len(agree) + 1)})]
    comparators = pd.read_csv(comparator_path, sep="\t", low_memory=False)
    comparators["method"] = comparators.method.replace(NAMES)
    comparators = comparators.loc[comparators.method.isin(METHODS) & comparators["rank"].le(max_rank)]
    curves.append(comparators[["method", "rank", "cumulative_agreement"]].rename(columns={"cumulative_agreement": "agreement"}))
    return pd.concat(curves, ignore_index=True)


def lr_plot(curves, width, height, legend="right"):
    curves = curves.assign(method=pd.Categorical(curves.method, categories=list(METHODS)))
    return (p9.ggplot(curves, p9.aes("rank", "agreement", color="method"))
            + p9.geom_hline(yintercept=0.5, linetype="dashed", color="#888888", size=0.3)
            + p9.geom_line(size=0.6)
            + p9.scale_color_manual(values=METHODS)
            + p9.scale_y_continuous(limits=(0.35, 1.0), breaks=[0.4, 0.6, 0.8, 1.0])
            + p9.labs(x="Rank (by p-value)", y="Cumulative agreement\nwith long reads")
            + theme(width, height, legend))


def evaluable_plot(summary, width, height):
    counts = summary.loc[summary.cutoff.eq(100), ["method", "n_evaluable"]].copy()
    counts["method"] = counts.method.replace(NAMES)
    counts = counts.loc[counts.method.isin(METHODS)].drop_duplicates("method")
    counts["method"] = pd.Categorical(counts.method, categories=list(METHODS)[::-1])
    return (p9.ggplot(counts, p9.aes("method", "n_evaluable", fill="method"))
            + p9.geom_col(width=0.7, show_legend=False)
            + p9.geom_text(p9.aes(label="n_evaluable"), format_string="{:,.0f}", ha="left", nudge_y=0.05, size=5.5)
            + p9.coord_flip() + p9.scale_fill_manual(values=METHODS)
            + p9.scale_y_log10(limits=(300, 2e5), breaks=[1e3, 1e4, 1e5], expand=(0, 0), labels=lambda values: [f"{v:,.0f}" for v in values])
            + p9.labs(x="", y="Long-read-evaluable tests")
            + theme(width, height) + p9.theme(panel_grid_major_y=p9.element_blank()))


def discoveries_plot(path, width, height):
    counts = pd.read_csv(path, sep="\t")
    counts = counts.loc[counts.unit.eq("Gene–cell-type pairs")].copy()
    counts["method"] = counts.method.replace(NAMES)
    totals = counts.groupby("method", as_index=False)["count"].sum()
    order = totals.sort_values("count").method.tolist()
    counts["method"] = pd.Categorical(counts.method, categories=order)
    totals["method"] = pd.Categorical(totals.method, categories=order)
    counts["event"] = pd.Categorical(counts.event, categories=list(EVENTS)[::-1])
    return (p9.ggplot(counts, p9.aes("method", "count", fill="event"))
            + p9.geom_col(width=0.7, color="white", size=0.15)
            + p9.geom_text(p9.aes("method", "count", label="count"), data=totals, inherit_aes=False, ha="left", size=5.5, format_string=" {:,}")
            + p9.coord_flip() + p9.scale_fill_manual(values=EVENTS, breaks=list(EVENTS))
            + p9.scale_y_continuous(expand=(0, 0, 0.18, 0))
            + p9.labs(x="", y="Gene–cell-type pairs with a discovery")
            + theme(width, height) + p9.theme(panel_grid_major_y=p9.element_blank(), legend_key_size=6))


def null_by_coverage(majiq_path, tealeaf_path):
    majiq = pd.read_csv(majiq_path, sep="\t")
    majiq.columns = ["bin", "tests", "observed", "null"]
    upper = majiq.bin.str.extract(r",\s*([0-9.]+)\]")[0].astype(float)
    majiq["bin"] = np.select([upper <= 1, upper <= 2, upper <= 3], ["≤ 1", "1–2", "2–3"], "> 3")
    majiq = majiq.groupby("bin").apply(lambda d: pd.Series({"rate": np.average(d["null"], weights=d.tests)}), include_groups=False).reset_index().assign(method="MAJIQ, BH 0.05")
    tealeaf = pd.read_csv(tealeaf_path, sep="\t")
    tealeaf = tealeaf.loc[tealeaf.scenario.eq("subject 20, no row variation")].copy()
    upper = tealeaf.abs_log2_coverage_ratio.str.extract(r",\s*([0-9.]+)\]")[0].astype(float)
    tealeaf["bin"] = np.select([upper <= 1, upper <= 2, upper <= 3], ["≤ 1", "1–2", "2–3"], "> 3")
    tealeaf = tealeaf.assign(rate=tealeaf.reject05_converged, method="Tealeaf " + tealeaf.structure + ", p < 0.05")[["bin", "rate", "method"]]
    table = pd.concat([majiq, tealeaf], ignore_index=True)
    table["bin"] = pd.Categorical(table.bin, categories=["≤ 1", "1–2", "2–3", "> 3"])
    return table


def null_plot(table, width, height):
    colors = {"MAJIQ, BH 0.05": METHODS["MAJIQ"], "Tealeaf internal, p < 0.05": METHODS["Tealeaf"], "Tealeaf terminal, p < 0.05": "#5B7A99"}
    return (p9.ggplot(table, p9.aes("bin", "rate", color="method", group="method"))
            + p9.geom_hline(yintercept=0.05, linetype="dashed", color="#888888", size=0.3)
            + p9.geom_line(size=0.6) + p9.geom_point(size=1.2)
            + p9.scale_color_manual(values=colors)
            + p9.labs(x="|log2 coverage ratio|", y="Called under the\ncoverage-only null")
            + theme(width, height, "bottom") + p9.guides(color=p9.guide_legend(ncol=1)))


def removed_table(paths):
    rows = []
    for dataset, path in paths.items():
        table = pd.read_csv(path, sep="\t")
        table["dataset"] = dataset
        rows.append(table)
    table = pd.concat(rows, ignore_index=True)
    table["method"] = table.method.replace(NAMES)
    table["fraction"] = table.removed / table.simulated
    table["run"] = table.run.astype(str).map({"full": "Full data", "0": "Half 1", "1": "Half 2"})
    return table


def removed_plot(table, width, height):
    order = ["Tealeaf", "LeafCutter", "scQuint", "rMATS", "MAJIQ"]
    table = table.assign(method=pd.Categorical(table.method, categories=order[::-1]), shown=table.fraction.clip(lower=1e-3))
    return (p9.ggplot(table, p9.aes("method", "shown", color="method", shape="dataset"))
            + p9.geom_hline(yintercept=0.0072, linetype="dashed", color="#888888", size=0.3)
            + p9.geom_point(size=1.6, position=p9.position_dodge(width=0.5), fill="white")
            + p9.coord_flip() + p9.scale_color_manual(values=METHODS, guide=None)
            + p9.scale_shape_manual(values={"Brain": "o", "GSE320507": "^"})
            + p9.scale_y_log10(breaks=[0.001, 0.01, 0.1, 1], labels=["≤0.1%", "1%", "10%", "100%"], limits=(0.001, 1))
            + p9.labs(x="", y="Recomputed tests removed")
            + theme(width, height, "bottom") + p9.theme(panel_grid_major_y=p9.element_blank()))


def rmats_null_plot(path, width, height):
    table = pd.read_csv(path, sep="\t", keep_default_na=False)
    table = table.loc[table.event_type.isin(["SE", "MXE"])].pivot_table(index=["contrast", "event_type"], columns="input", values="bh05").reset_index()
    limit = max(table.real.max(), table.null.max()) * 1.1
    return (p9.ggplot(table, p9.aes("real", "null", color="event_type"))
            + p9.geom_abline(intercept=0, slope=1, linetype="dashed", color="#888888", size=0.3)
            + p9.geom_point(size=1.6)
            + p9.scale_x_log10(limits=(10, limit)) + p9.scale_y_log10(limits=(10, limit))
            + p9.scale_color_manual(values={"SE": EVENTS["cassette exon"], "MXE": EVENTS["mutually exclusive exons"]})
            + p9.labs(x="rMATS BH calls, data", y="rMATS BH calls,\ncoverage-only null")
            + theme(width, height, "right"))


def dataset_curves(paths, max_rank=200):
    curves = []
    for dataset, path in paths.items():
        table = pd.read_csv(path, sep="\t", low_memory=False)
        if "scope" in table:
            table = table.loc[table.scope.eq("all tested") & table.ranking.eq("continuous raw")]
        table = table.sort_values("rank").head(max_rank)
        agree = table.pooled_replicated.astype(str).str.lower().eq("true").to_numpy(float)
        curves.append(pd.DataFrame({"dataset": dataset, "rank": np.arange(1, len(agree) + 1), "agreement": np.cumsum(agree) / np.arange(1, len(agree) + 1)}))
    return pd.concat(curves, ignore_index=True)


def dataset_plot(curves, width, height):
    colors = {"GSE320507": "#17324D", "GSE327907": "#3E7CB1", "GSE334840": "#9DBCD4"}
    return (p9.ggplot(curves, p9.aes("rank", "agreement", color="dataset"))
            + p9.geom_hline(yintercept=0.5, linetype="dashed", color="#888888", size=0.3)
            + p9.geom_line(size=0.6) + p9.scale_color_manual(values=colors)
            + p9.scale_y_continuous(limits=(0.35, 1.0), breaks=[0.4, 0.6, 0.8, 1.0])
            + p9.labs(x="Rank (by p-value)", y="Tealeaf cumulative\nagreement")
            + theme(width, height, "right"))


def disease_plot(path, width, height):
    table = pd.read_csv(path, sep="\t")
    table["contrast"] = table.kind.replace({"P26 vs P12": "Random 26 vs 12", "LTCS vs HC+CONV": "LTCS vs HC + CONV"})
    order = ["Random 26 vs 12", "CONV vs HC", "LTCS vs CONV", "LTCS vs HC + CONV", "LTCS vs HC"]
    long = table.melt(id_vars=["contrast", "bh_within_kind"], value_vars=["below_01", "below_001"], var_name="level", value_name="rate")
    long["level"] = long.level.map({"below_01": "p < 0.01", "below_001": "p < 0.001"})
    long["contrast"] = pd.Categorical(long.contrast, categories=order[::-1])
    nominal = pd.DataFrame({"level": ["p < 0.01", "p < 0.001"], "rate": [0.01, 0.001]})
    labels = long.loc[long.level.eq("p < 0.01")].assign(label=lambda d: d.bh_within_kind.map(lambda n: f"{n} BH"))
    return (p9.ggplot(long, p9.aes("contrast", "rate"))
            + p9.geom_col(fill=METHODS["Tealeaf"], width=0.7)
            + p9.geom_hline(p9.aes(yintercept="rate"), data=nominal, linetype="dashed", color="#D95F0E", size=0.4)
            + p9.geom_text(p9.aes(label="label"), data=labels, ha="left", nudge_y=0.0006, size=5.5)
            + p9.facet_wrap("~level", scales="free_x") + p9.coord_flip()
            + p9.scale_y_continuous(expand=(0, 0, 0.3, 0), breaks=breaks_extended(n=3), labels=lambda values: [f"{100 * v:g}%" for v in values])
            + p9.labs(x="", y="Fraction of tests")
            + theme(width, height) + p9.theme(panel_grid_major_y=p9.element_blank()))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, default=Path("docs/figures/paper"))
    args = parser.parse_args()
    out = args.output_dir
    out.mkdir(parents=True, exist_ok=True)
    brain, pbmc = Path("analyses/coverage_null_filter_leaf5/brain"), Path("analyses/coverage_null_filter_leaf5/pbmc")

    table, span, gap = read_end_metaplot(sorted(Path(".cache/read_end_profiles").glob("library_*.json")))
    save(metaplot(table, span, gap, 3.3, 1.9), out, "fig1_profiles")

    brain_labels = {"MAJIQ": ("MAJIQ", 0.08, 0.02), "LeafCutter": ("LeafCutter", 0.02, -0.12), "scQuint": ("scQuint", -0.08, 0.07), "rMATS": ("rMATS", 0.08, 0.0), "SUPPA2-PA": ("SUPPA2-PA", 0.08, 0.04), "SUPPA2": ("SUPPA2", 0.08, -0.05)}
    save(split_scatter(split_both(brain / "split_filtered/split_metrics.tsv"), 3.0, 3.1, brain_labels), out, "fig2_split")
    save(lr_plot(lr_curves(brain / "tealeaf_lr_rank_filtered.tsv.gz", brain / "sign_agreement_by_rank_filtered.tsv"), 3.3, 2.0), out, "fig2_lr")
    summary = pd.read_csv(brain / "lr_summary.tsv", sep="\t")
    suppa = pd.read_csv(brain / "suppa2_exact_lr_summary.tsv", sep="\t")
    save(evaluable_plot(pd.concat([summary.loc[summary.table.eq("filtered")], suppa.loc[suppa.method.str.startswith("SUPPA2")]], ignore_index=True), 2.6, 2.0), out, "fig2_evaluable")
    save(discoveries_plot("docs/figures/discoveries_by_event.tsv", 4.0, 2.0), out, "fig2_discoveries")

    save(null_plot(null_by_coverage("analyses/pbmc/GSE320507_comparators/majiq_coverage_null/by_coverage_ratio.tsv", "analyses/pbmc/GSE320507_comparators/coverage_null/tealeaf_count_null_by_coverage.tsv"), 2.4, 2.3), out, "fig4_null_coverage")
    save(removed_plot(removed_table({"Brain": brain / "removed_tests.tsv", "GSE320507": pbmc / "removed_tests.tsv"}), 2.3, 2.3), out, "fig4_removed")
    save(rmats_null_plot("analyses/pbmc/GSE320507_comparators/coverage_null/rmats_null_vs_real_4contrasts.tsv", 2.2, 2.0), out, "fig4_rmats")

    pbmc_labels = {"MAJIQ": ("MAJIQ", 0.06, 0.08), "LeafCutter": ("LeafCutter", -0.08, 0.12), "rMATS": ("rMATS", 0.04, -0.14), "scQuint": ("scQuint", 0.08, -0.08), "SUPPA2, SUPPA2-PA": ("SUPPA2", 0.06, 0.12)}
    save(split_scatter(split_both(pbmc / "split_filtered/split_metrics.tsv"), 3.0, 3.1, pbmc_labels), out, "fig5_split")
    save(lr_plot(lr_curves(pbmc / "tealeaf_lr_rank_filtered.tsv.gz", pbmc / "sign_agreement_by_rank_filtered.tsv"), 3.3, 2.0), out, "fig5_lr")
    save(dataset_plot(dataset_curves({"GSE320507": "analyses/pbmc/GSE320507_lr_groupfix/lr_rank.tsv.gz", "GSE327907": "analyses/pbmc/GSE327907_lr/lr_rank.tsv.gz", "GSE334840": "analyses/pbmc/GSE334840_lr/lr_rank.tsv.gz"}), 3.0, 1.8), out, "fig5_datasets")
    save(disease_plot("analyses/pbmc/GSE320507_disease/calibration.tsv", 3.6, 1.8), out, "fig5_disease")


if __name__ == "__main__":
    main()
