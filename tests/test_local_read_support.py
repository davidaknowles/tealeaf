from dataclasses import asdict

import pandas as pd
import pysam
import pytest

from extra_scripts.audit_event_local_read_support import decode_event, event_envelope, summarize_signature_counts
from tealeaf.sc.local_read_support import event_local_read_contrast, alignment_local_signature, collect_indexed_local_read_support


def cassette(strand="+"):
    annotation = dict(chromosome="chr1", strand=strand, transcripts=dict(inc=((10, 20), (30, 40), (50, 60)), exc=((10, 20), (50, 60))))
    return event_local_read_contrast("event", "gene.1", annotation, ("inc",), ("exc",), (20, 50))


@pytest.mark.parametrize("strand", ["+", "-"])
def test_cassette_features_are_local_and_genomic_on_either_strand(strand):
    event = cassette(strand)
    assert event.included_exons == ((30, 40),)
    assert event.excluded_exons == ()
    assert event.included_junctions == ((20, 30), (40, 50))
    assert event.excluded_junctions == ((20, 50),)
    assert (event.gene_start, event.gene_end) == (10, 60)
    assert alignment_local_signature(((30, 35),), (), event) == 17
    assert alignment_local_signature(((10, 20), (50, 60)), ((20, 50),), event) == 8
    assert alignment_local_signature(((10, 15),), (), event) == 0


def test_feature_requires_class_unanimity_not_global_isoform_linkage():
    annotation = dict(chromosome="chr1", strand="+", transcripts=dict(i1=((10, 20), (30, 40), (50, 60)), i2=((10, 20), (50, 60)), exc=((10, 20), (50, 60))))
    event = event_local_read_contrast("e", "g", annotation, ("i1", "i2"), ("exc",), (20, 50))
    assert not event.included_exons and not event.excluded_exons
    assert not event.included_junctions and not event.excluded_junctions
    assert alignment_local_signature(((30, 35),), (), event) == 16


def test_exon_features_clip_to_definition_not_entire_transcript():
    annotation = dict(chromosome="chr1", strand="+", transcripts=dict(inc=((1, 10), (30, 40), (80, 100)), exc=((1, 20), (50, 60), (80, 100))))
    event = event_local_read_contrast("e", "g", annotation, ("inc",), ("exc",), (32, 55))
    assert event.included_exons == ((32, 40),)
    assert event.excluded_exons == ((50, 55),)
    assert alignment_local_signature(((10, 15),), (), event) == 0


@pytest.mark.parametrize("included,excluded,span", [((), ("exc",), (20, 50)), (("inc",), ("inc",), (20, 50)), (("missing",), ("exc",), (20, 50)), (("inc",), ("exc",), (50, 20))])
def test_invalid_event_classes_fail(included, excluded, span):
    with pytest.raises(ValueError):
        event_local_read_contrast("e", "g", dict(chromosome="chr1", strand="+", transcripts=dict(inc=((10, 20),), exc=((10, 20),))), included, excluded, span)


def test_json_feature_round_trip_supports_junction_set_operations():
    import json
    event = decode_event(json.loads(json.dumps(asdict(cassette()))))
    assert event == cassette()
    assert alignment_local_signature(((10, 20), (30, 40)), ((20, 30),), event) == 21


@pytest.mark.parametrize("kind,coordinates", [("SE", "20-30:40-50"), ("A3", "20-30:20-40"), ("AF", "10:20-50:30:40-50"), ("MX", "20-30:40-60:20-45:50-60")])
def test_suppa_envelope_does_not_use_gene_span(kind, coordinates):
    chromosome, strand, span = event_envelope(f"g;{kind}:chr1:{coordinates}:-")
    integers = [int(value) for value in coordinates.replace(":", "-").split("-")]
    assert chromosome == "chr1" and strand == "-"
    assert span == (min(integers) - 1, max(integers))


def write_bam(path, records):
    with pysam.AlignmentFile(path, "wb", header={"HD": {"VN": "1.6", "SO": "coordinate"}, "SQ": [{"SN": "chr1", "LN": 1000}]}) as target:
        for index, (start, cigar, umi, tags, flag) in enumerate(sorted(records, key=lambda row: row[0])):
            alignment = pysam.AlignedSegment()
            alignment.query_name = f"read{index}"
            length = sum(size for operation, size in cigar if operation in (0, 1, 4, 7, 8))
            alignment.query_sequence = "A" * length
            alignment.flag = flag
            alignment.reference_id = 0
            alignment.reference_start = start
            alignment.mapping_quality = 60
            alignment.cigar = cigar
            for tag, value in (dict(CB="AAAA", UB=umi, NH=1, GX="gene.1") | tags).items():
                alignment.set_tag(tag, value)
            target.write(alignment)
    pysam.index(str(path))


def test_collector_combines_read_loci_and_retains_umi_conflicts(tmp_path):
    path = tmp_path / "reads.bam"
    write_bam(path, [(30, ((0, 5),), "inc", {}, 0), (31, ((0, 5),), "inc", {}, 0), (30, ((0, 5),), "conflict", {}, 0), (10, ((0, 10), (3, 30), (0, 10)), "conflict", {}, 0), (10, ((0, 5),), "distal", {}, 0)])
    counts, filters = collect_indexed_local_read_support(path, str(path) + ".bai", {"AAAA": ("s", "ct", "poly(dT)")}, (cassette(),))
    assert counts[("event", "s", "ct", "poly(dT)", 17)] == 1
    assert counts[("event", "s", "ct", "poly(dT)", 25)] == 1
    assert counts[("event", "s", "ct", "poly(dT)", 0)] == 1
    assert sum(counts.values()) == 3
    assert filters["fetched_alignments"] == 5


def test_collector_does_not_count_intron_only_overlap(tmp_path):
    path = tmp_path / "reads.bam"
    write_bam(path, [(22, ((0, 4),), "intronic", {}, 0)])
    counts, _ = collect_indexed_local_read_support(path, str(path) + ".bai", {"AAAA": ("s", "ct", "random hexamer")}, (cassette(),))
    assert not counts


def test_collector_requires_unique_primary_gene_and_umi_tags(tmp_path):
    path = tmp_path / "reads.bam"
    write_bam(path, [(30, ((0, 5),), "good", {}, 0), (30, ((0, 5),), "multi", {"NH": 2}, 0), (30, ((0, 5),), "wrong", {"GX": "other"}, 0), (30, ((0, 5),), "ambiguous", {"GX": "gene.1;other"}, 0), (30, ((0, 5),), "secondary", {}, 256), (30, ((0, 5),), "noUB", {"UB": None}, 0), (30, ((0, 5),), "unknownCB", {"CB": "TTTT"}, 0)])
    counts, filters = collect_indexed_local_read_support(path, str(path) + ".bai", {"AAAA": ("s", "ct", "poly(dT)")}, (cassette(),))
    assert sum(counts.values()) == 1
    assert filters["multimapped"] == 1
    assert filters["nonunique_or_other_gene"] == 2
    assert filters["missing_tags"] == 1
    assert filters["nonprimary"] == 1
    assert filters["unrequested_barcode"] == 1


def test_signature_summary_does_not_turn_conflict_into_support():
    table = pd.DataFrame(dict(feature_id=["e"] * 5, subject=["s"] * 5, cell_type=["ct"] * 5, primer=["RH"] * 5, signature=[17, 8, 25, 16, 0], count=[2, 3, 7, 11, 13]))
    row = summarize_signature_counts(table).iloc[0]
    assert (row.included_only, row.excluded_only, row.conflicting, row.local_nondiscriminating, row.distal_only, row.gene_keys) == (2, 3, 7, 11, 13, 36)
    original = table.copy(deep=True)
    summarize_signature_counts(table)
    pd.testing.assert_frame_equal(table, original)
