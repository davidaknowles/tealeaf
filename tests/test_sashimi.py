from collections import Counter

import pandas as pd
import pysam

from tealeaf.sc.sashimi import SashimiEvent, collect_bam_event_support, write_ggsashimi_inputs


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
