from collections import Counter

import pandas as pd
import pysam

from tealeaf.sc.sashimi import SashimiEvent, collect_bam_event_support, combine_path_usage_and_coverage, select_strongest_significant_contrasts, summarize_path_ordered_coverage, summarize_path_ordered_psi, summarize_path_usage, summarize_subject_feature_evidence, variable_path_features, write_ggsashimi_inputs


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


def test_collects_junction_competitors_outside_event_interval(tmp_path):
    bam_path = tmp_path / "reads.bam"
    header = {"HD": {"VN": "1.6"}, "SQ": [{"SN": "chr1", "LN": 1000}]}
    with pysam.AlignmentFile(bam_path, "wb", header=header) as target:
        alignment = pysam.AlignedSegment()
        alignment.query_name = "read"
        alignment.query_sequence = "A" * 20
        alignment.flag = 0
        alignment.reference_id = 0
        alignment.reference_start = 50
        alignment.mapping_quality = 60
        alignment.cigar = ((0, 10), (3, 150), (0, 10))
        alignment.query_qualities = pysam.qualitystring_to_array("I" * 20)
        alignment.set_tag("CB", "AAAA_CCCC")
        target.write(alignment)
    event = SashimiEvent("event", "chr1", 100, 230, "+", frozenset({210}))
    groups = {"AAAACCCC": ("subject", "cell", "poly(dT)")}
    _, junctions, _ = collect_bam_event_support(bam_path, groups, (event,))
    assert junctions == Counter({("event", "subject", "cell", "poly(dT)", 60, 210): 1})


def test_collect_can_deduplicate_cell_umis(tmp_path):
    bam_path = tmp_path / "reads.bam"
    header = {"HD": {"VN": "1.6"}, "SQ": [{"SN": "chr1", "LN": 1000}]}
    with pysam.AlignmentFile(bam_path, "wb", header=header) as target:
        for name in ("read1", "read2"):
            alignment = pysam.AlignedSegment()
            alignment.query_name = name
            alignment.query_sequence = "A" * 20
            alignment.flag = 0
            alignment.reference_id = 0
            alignment.reference_start = 100
            alignment.mapping_quality = 60
            alignment.cigar = ((0, 20),)
            alignment.query_qualities = pysam.qualitystring_to_array("I" * 20)
            alignment.set_tag("CB", "AAAA_CCCC")
            alignment.set_tag("UB", "UMI")
            alignment.set_tag("NH", 1)
            target.write(alignment)
    event = SashimiEvent("event", "chr1", 90, 130, "+")
    groups = {"AAAACCCC": ("subject", "cell", "poly(dT)")}
    exon_blocks, _, event_umis = collect_bam_event_support(bam_path, groups, (event,), deduplicate=True)
    assert sum(exon_blocks.values()) == 1
    assert sum(event_umis.values()) == 1


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


def test_summarize_path_ordered_coverage():
    event = SashimiEvent("event", "chr1", 90, 160, "+")
    groups = (("A", "poly(dT)"), ("A", "random hexamer"), ("B", "poly(dT)"), ("B", "random hexamer"))
    sizes = pd.DataFrame([{"subject": "s1", "cell_type": cell_type, "primer": primer, "cells": 10} for cell_type, primer in groups])
    exon_blocks = Counter({("event", "s1", "A", "poly(dT)", 100, 110): 1})
    junctions = Counter({("event", "s1", "A", "poly(dT)", 110, 120): 2})
    paths = (((100, 110), (120, 130)), ((100, 110), (140, 150)))
    coverage = summarize_path_ordered_coverage(event, paths, exon_blocks, junctions, sizes, ("A", "B"))
    a_polydt = coverage[(coverage["cell_type"] == "A") & (coverage["primer"] == "poly(dT)")]
    assert a_polydt.loc[a_polydt["feature_id"] == "P1:E1", "raw_value"].iloc[0] == 100.0
    assert a_polydt.loc[a_polydt["feature_id"] == "P1:J1", "raw_value"].iloc[0] == 200.0
    assert set(coverage.loc[coverage["coordinate"] == "chr1:101-110", "feature_id"]) == {"P1:E1", "P2:E1"}
    negative = summarize_path_ordered_coverage(SashimiEvent("event", "chr1", 90, 160, "-"), paths, exon_blocks, junctions, sizes, ("A", "B"))
    negative_path_one = negative[(negative["path_number"] == 1) & (negative["cell_type"] == "A") & (negative["primer"] == "poly(dT)")].sort_values("column_order")
    assert negative_path_one["coordinate"].tolist() == ["chr1:121-130", "chr1:111-121", "chr1:101-110"]
    assert negative_path_one["feature"].tolist() == ["E2", "J2", "E3"]
    usage = pd.DataFrame([
        {"test_id": "event", "cell_type": cell_type, "path": f"Path {path}", "path_number": path, "n_subjects": 1, "mean_proportion": value, "sd_proportion": 0.0, "se_proportion": 0.0}
        for cell_type, values in (("A", (0.25, 0.75)), ("B", (0.5, 0.5)))
        for path, value in enumerate(values, start=1)
    ])
    matrix = combine_path_usage_and_coverage(usage, coverage, "event")
    columns = matrix[["feature_id", "column_order"]].drop_duplicates().sort_values("column_order")
    assert columns["feature_id"].tolist() == ["P1:path", "P1:E1", "P1:J1", "P1:E2", "P2:path", "P2:E1", "P2:J2", "P2:E3"]


def test_summarize_path_ordered_psi():
    event = SashimiEvent("event", "chr1", 90, 180, "+")
    paths = (((100, 110), (120, 130), (160, 170)), ((100, 110), (140, 150), (160, 170)))
    exon_blocks = Counter()
    junctions = Counter()
    for primer in ("poly(dT)", "random hexamer"):
        for exon, count in (((100, 110), 10), ((120, 130), 8), ((140, 150), 2), ((160, 170), 10)):
            exon_blocks[("event", "s1", "A", primer, *exon)] = count
        for junction, count in (((110, 120), 8), ((110, 140), 2), ((130, 160), 8), ((150, 160), 2)):
            junctions[("event", "s1", "A", primer, *junction)] = count
    psi = summarize_path_ordered_psi(event, paths, exon_blocks, junctions, ("A",))
    polydt = psi[psi["primer"] == "poly(dT)"]
    values = polydt.drop_duplicates(["feature_type", "feature"]).set_index(["feature_type", "feature"])["raw_value"]
    assert values[("exon", "E2")] == 0.8
    assert values[("exon", "E3")] == 0.2
    assert values[("junction", "J1")] == 0.8
    assert values[("junction", "J2")] == 0.2
    assert values[("junction", "J3")] == 0.8
    assert values[("junction", "J4")] == 0.2
    assert set(polydt.loc[polydt["feature_type"] == "exon", "feature"]) == {"E2", "E3"}


def test_psi_matrix_retains_path_with_no_nonconstitutive_components():
    event = SashimiEvent("event", "chr1", 90, 180, "+")
    paths = (((160, 170),), ((100, 110), (160, 170)))
    psi = summarize_path_ordered_psi(event, paths, Counter(), Counter(), ("A",))
    usage = pd.DataFrame([
        {"test_id": "event", "cell_type": "A", "path": "Path 1", "path_number": 1, "n_subjects": 1, "mean_proportion": 0.25, "sd_proportion": 0.0, "se_proportion": 0.0},
        {"test_id": "event", "cell_type": "A", "path": "Path 2", "path_number": 2, "n_subjects": 1, "mean_proportion": 0.75, "sd_proportion": 0.0, "se_proportion": 0.0},
    ])
    matrix = combine_path_usage_and_coverage(usage, psi, "event")
    columns = matrix[["feature_id", "column_order"]].drop_duplicates().sort_values("column_order")
    assert columns["feature_id"].tolist() == ["P1:path", "P2:path", "P2:E1", "P2:J1"]


def test_variable_path_features_split_alternative_acceptor_segment():
    event = SashimiEvent("event", "chr1", 100, 300, "+")
    paths = (((100, 150), (200, 250)), ((100, 150), (220, 250)))
    features = variable_path_features(event, paths)
    exon = features[features["feature_type"] == "exon"].iloc[0]
    assert (exon.start, exon.end, exon.path_numbers) == (200, 220, (1,))
    junctions = features[features["feature_type"] == "junction"]
    assert set(zip(junctions.start, junctions.end, junctions.path_numbers)) == {(150, 200, (1,)), (150, 220, (2,))}


def test_summarize_subject_feature_evidence_uses_endpoint_depth_and_splice_sites():
    event = SashimiEvent("event", "chr1", 100, 400, "+")
    paths = (((100, 150), (200, 250), (350, 400)), ((100, 150), (280, 300), (350, 400)))
    usage = pd.DataFrame([
        {"subject": "s1", "cell_type": "A", "path_number": 1, "proportion": 0.8},
        {"subject": "s1", "cell_type": "A", "path_number": 2, "proportion": 0.2},
    ])
    exon_blocks = Counter()
    junctions = Counter()
    for primer in ("poly(dT)", "random hexamer"):
        for interval, count in (((100, 150), 10), ((200, 250), 8), ((280, 300), 2), ((350, 400), 10)):
            exon_blocks[("event", "s1", "A", primer, *interval)] = count
        for interval, count in (((150, 200), 8), ((150, 280), 2), ((250, 350), 8), ((300, 350), 2)):
            junctions[("event", "s1", "A", primer, *interval)] = count
    evidence = summarize_subject_feature_evidence(event, paths, ((100, 150), (350, 400)), exon_blocks, junctions, usage)
    polydt = evidence[evidence["primer"] == "poly(dT)"].set_index("feature_id")
    assert polydt.loc["E:200-250", "fitted_inclusion"] == 0.8
    assert polydt.loc["E:200-250", "observed_inclusion"] == 0.8
    assert polydt.loc["J:150-200", "observed_inclusion"] == 0.8
    assert polydt.loc["J:150-200", "evidence_denominator"] == 10
    assert polydt.loc["J:150-200", "endpoint_normalized_signal"] == 0.8
