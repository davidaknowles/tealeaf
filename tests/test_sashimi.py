from collections import Counter

import pandas as pd
import pysam

from tealeaf.sc.sashimi import SashimiEvent, collect_bam_event_support, select_strongest_significant_contrasts, summarize_path_usage, write_ggsashimi_inputs


def test_collect_and_write_sashimi_support(tmp_path):
    bam_path = tmp_path / "reads.bam"
    header = {"HD": {"VN": "1.6"}, "SQ": [{"SN": "chr1", "LN": 1000}]}
    with pysam.AlignmentFile(bam_path, "wb", header=header) as target:
        alignment = pysam.AlignedSegment()
        alignment.query_name = "read"
        alignment.query_sequence = "A" * 20
        alignment.flag = 0
        alignment.reference_id = 0
        alignment.reference_start = 100
        alignment.mapping_quality = 60
        alignment.cigar = ((0, 10), (3, 100), (0, 10))
        alignment.query_qualities = pysam.qualitystring_to_array("I" * 20)
        alignment.set_tag("CB", "AAAA_CCCC")
        target.write(alignment)

    event = SashimiEvent("event", "chr1", 90, 230, "+")
    groups = {"AAAACCCC": ("subject", "cell", "poly(dT)")}
    exon_blocks, junctions, event_umis = collect_bam_event_support(bam_path, groups, (event,))
    assert exon_blocks == Counter({("event", "subject", "cell", "poly(dT)", 100, 110): 1, ("event", "subject", "cell", "poly(dT)", 210, 220): 1})
    assert junctions == Counter({("event", "subject", "cell", "poly(dT)", 110, 210): 1})
    assert event_umis == Counter({("event", "subject", "cell", "poly(dT)"): 1})

    sizes = pd.DataFrame([
        {"subject": "subject", "cell_type": "cell", "primer": "poly(dT)", "cells": 10},
        {"subject": "subject", "cell_type": "cell", "primer": "random hexamer", "cells": 10},
    ])
    all_junctions = junctions + Counter({("event", "subject", "cell", "random hexamer", 110, 210): 2})
    intron_path, exon_path, _ = write_ggsashimi_inputs(tmp_path / "output", event, exon_blocks, all_junctions, sizes, ("cell",))
    introns = pd.read_csv(intron_path, sep=" ")
    assert introns.loc[0, "cell_polydT_0"] == 100
    assert introns.loc[0, "cell_randomhexamer_0"] == 200
    assert not pd.read_csv(exon_path, sep=" ").empty


def test_select_strongest_significant_contrasts(tmp_path):
    fits = pd.DataFrame([
        {"test_id": "t1", "block_id": "b1", "level_a": "A", "level_b": "B", "mean_difference_norm": 0.4},
        {"test_id": "t2", "block_id": "b1", "level_a": "A", "level_b": "C", "mean_difference_norm": 0.9},
        {"test_id": "t3", "block_id": "b2", "level_a": "D", "level_b": "E", "mean_difference_norm": 0.7},
    ])
    catalog = pd.DataFrame([
        {"test_id": "t1", "block_id": "b1", "statistic": 4.0, "empirical_p_value": 0.01, "fdr": 0.02},
        {"test_id": "t2", "block_id": "b1", "statistic": 9.0, "empirical_p_value": 0.001, "fdr": 0.10},
        {"test_id": "t3", "block_id": "b2", "statistic": 5.0, "empirical_p_value": 0.005, "fdr": 0.03},
    ])
    fit_path = tmp_path / "fits.tsv"
    catalog_path = tmp_path / "catalog.tsv"
    fits.to_csv(fit_path, sep="\t", index=False)
    catalog.to_csv(catalog_path, sep="\t", index=False)
    selected = select_strongest_significant_contrasts([fit_path], catalog_path, ["b1", "b2"])
    assert selected["test_id"].tolist() == ["t1", "t3"]
    assert selected["mean_difference_norm"].tolist() == [0.4, 0.7]


def test_summarize_path_usage():
    usage = pd.DataFrame([
        {"test_id": "t1", "cell_type": "A", "path": "Path 1", "path_number": 1, "proportion": 0.2},
        {"test_id": "t1", "cell_type": "A", "path": "Path 1", "path_number": 1, "proportion": 0.4},
        {"test_id": "t1", "cell_type": "A", "path": "Path 2", "path_number": 2, "proportion": 0.8},
        {"test_id": "t1", "cell_type": "A", "path": "Path 2", "path_number": 2, "proportion": 0.6},
    ])
    summary = summarize_path_usage(usage)
    path_one = summary[summary["path_number"] == 1].iloc[0]
    assert path_one.n_subjects == 2
    assert abs(path_one.mean_proportion - 0.3) < 1e-12
    assert abs(path_one.se_proportion - 0.1) < 1e-12
