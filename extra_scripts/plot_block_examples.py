#!/usr/bin/env python3
"""Figures for example block tests (tables from export_block_examples.py).

Per example: the block's path exon chains (introns longer than 300 bases
shortened), path usage per subject in the Parse data (molecule-share
composition of each subject/cell-type row, with the mean over subjects), and
path proportions in each long-read replicate. Needs plotnine.
"""

import argparse
from pathlib import Path

import pandas as pd
import plotnine as p9

PALETTE = ["#168C8C", "#E6813B", "#17324D", "#9C6ADE"]
TYPES = ["#2C7FB8", "#D95F0E"]


def theme(base):
    return p9.theme_classic(base_size=base) + p9.theme(plot_title=p9.element_text(size=base), legend_title=p9.element_text(size=base), legend_text=p9.element_text(size=base - 1), legend_key_size=base, axis_line=p9.element_line(size=0.4), axis_ticks=p9.element_line(size=0.4))


def figure(example, parse, long_read, exons, base=11, names=None):
    names = names or {}
    levels = [example.level_a, example.level_b]
    shown = [names.get(level, level) for level in levels]
    parse = parse.assign(cell_type=parse.cell_type.map(lambda level: names.get(level, level)))
    long_read = long_read.assign(cell_type=long_read.cell_type.map(lambda level: names.get(level, level)))
    paths = list(exons.sort_values("path").drop_duplicates("path").path_label)
    exons = exons.assign(path_label=pd.Categorical(exons.path_label, categories=paths[::-1]))
    exons["y"] = exons.path_label.cat.codes
    structure = (p9.ggplot(exons)
                 + p9.geom_rect(p9.aes(xmin="x0", xmax="x1", ymin="y - 0.3", ymax="y + 0.3", fill="factor(path)"))
                 + p9.scale_fill_manual(values=PALETTE, guide=None)
                 + p9.scale_y_continuous(breaks=list(range(len(paths))), labels=paths[::-1])
                 + p9.labs(x=f"{example.chromosome} ({example.strand} strand; long introns shortened)", y="", title=f"{example.label}, {example.structure} block, p = {example.p_value:.1e}")
                 + theme(base) + p9.theme(axis_text_x=p9.element_blank(), axis_ticks_major_x=p9.element_blank(), axis_line_x=p9.element_blank(), axis_title_x=p9.element_text(size=base - 1)))
    if example.strand == "-":
        structure += p9.scale_x_reverse()
    parse = parse.assign(cell_type=pd.Categorical(parse.cell_type, categories=shown), path_label=pd.Categorical(parse.path_label, categories=paths))
    means = parse.groupby(["cell_type", "path_label"], observed=True).proportion.mean().reset_index()
    dodge = p9.position_dodge(width=0.6)
    parse_plot = (p9.ggplot(parse, p9.aes("path_label", "proportion", color="cell_type"))
                  + p9.geom_point(position=p9.position_jitterdodge(jitter_width=0.15, dodge_width=0.6, random_state=1), alpha=0.6, size=base / 11 * 1.5, stroke=0)
                  + p9.geom_point(data=means, position=dodge, shape="_", size=base / 11 * 12, stroke=base / 11 * 1.2)
                  + p9.scale_color_manual(values=TYPES, name="")
                  + p9.scale_y_continuous(limits=(0, 1))
                  + p9.labs(x="", y="Path usage", title="Parse, per subject (bar: mean)")
                  + theme(base) + p9.theme(legend_position="bottom", legend_direction="vertical"))
    long_read = long_read.assign(cell_type=pd.Categorical(long_read.cell_type, categories=shown), path_label=pd.Categorical(long_read.path_label, categories=paths), replicate=long_read.replicate.astype(str))
    long_plot = (p9.ggplot(long_read, p9.aes("path_label", "proportion", color="cell_type", shape="replicate"))
                 + p9.geom_point(position=dodge, size=base / 11 * 3)
                 + p9.scale_color_manual(values=TYPES, guide=None)
                 + p9.scale_shape_manual(values=["o", "^"], name="Long-read replicate")
                 + p9.scale_y_continuous(limits=(0, 1))
                 + p9.labs(x="", y="Transcript share", title=f"Long reads ({'agrees' if example.long_read_agrees else 'disagrees'})")
                 + theme(base) + p9.theme(legend_position="bottom", legend_direction="vertical"))
    return structure / (parse_plot | long_plot)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--tables", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--width", type=float, default=8.)
    parser.add_argument("--height", type=float, default=7.)
    parser.add_argument("--base-size", type=float, default=11.)
    parser.add_argument("--level-label", action="append", default=[], help="LEVEL=label for cell-type names")
    args = parser.parse_args()
    names = dict(item.split("=", 1) for item in args.level_label)
    tables = {name: pd.read_csv(args.tables / f"{name}.tsv", sep="\t") for name in ("examples", "parse_usage", "long_read_usage", "exons")}
    args.output_dir.mkdir(parents=True, exist_ok=True)
    for example in tables["examples"].itertuples():
        select = lambda table: table.loc[table.test_id.eq(example.test_id)]
        name = f"example_{example.gene_name}_{example.block_id.split(':')[-1]}"
        for suffix, options in ((".pdf", {}), (".png", {"dpi": 110})):
            # a composed plot can be saved once, so it is rebuilt per format
            plot = figure(example, select(tables["parse_usage"]), select(tables["long_read_usage"]), select(tables["exons"]), args.base_size, names)
            plot.save(args.output_dir / f"{name}{suffix}", width=args.width, height=args.height, verbose=False, **options)
        print(name, flush=True)


if __name__ == "__main__":
    main()
