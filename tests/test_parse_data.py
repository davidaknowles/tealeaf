from pathlib import Path

import pytest

from tealeaf.data.parse import (
    build_rt_primer_lookup,
    demultiplex_parse_transcript_reads,
    parse_primer_pairs,
)


def test_parse_primer_pairs_preserves_batch_and_ligation_prefix():
    rt = ["AAAA", "CCCC", "GGGG", "TTTT"]
    cells = ["b1:ACGTAAAA", "b1:ACGTGGGG", "b2:TGCAAAAA"]
    assert parse_primer_pairs(cells, rt) == [
        ("b1:ACGTAAAA", "b1:ACGTAAAA", "b1:ACGTGGGG"),
        ("b2:TGCAAAAA", "b2:TGCAAAAA", "b2:TGCAGGGG"),
    ]


def test_parse_primer_pairs_rejects_ambiguous_variable_suffixes():
    with pytest.raises(ValueError, match="ambiguous"):
        parse_primer_pairs(["prefixAAAA"], ["AAAA", "CCCC", "AA", "GG"])


def _write_fastq(path, records):
    path.write_text("".join(
        f"@{name}/{read}\n{sequence}\n+\n{'I' * len(sequence)}\n"
        for name, read, sequence in records
    ))


def test_rt_primer_lookup_only_corrects_unambiguous_primer_neighbors():
    exact, corrected = build_rt_primer_lookup(["AAAA", "TTTT", "CCCC", "GGGG"])
    assert exact["AAAA"] == "polydT"
    assert exact["CCCC"] == "ranhex"
    assert corrected["AAAT"] == "polydT"

    _, ambiguous = build_rt_primer_lookup(["AAAA", "TTTT", "AAAC", "GGGG"])
    assert "AAAG" not in ambiguous


def test_demultiplex_parse_transcript_reads(tmp_path):
    read1 = tmp_path / "r1.fastq"
    read2 = tmp_path / "r2.fastq"
    _write_fastq(read1, [
        ("exact_poly", 1, "ACGT"),
        ("corrected_hex", 1, "TGCA"),
        ("unknown", 1, "CCCC"),
    ])
    _write_fastq(read2, [
        ("exact_poly", 2, "NNAAAA"),
        ("corrected_hex", 2, "NNGGGA"),
        ("unknown", 2, "NNTATA"),
    ])
    polydt = tmp_path / "polydt.fastq"
    ranhex = tmp_path / "ranhex.fastq"
    stats = demultiplex_parse_transcript_reads(
        [read1],
        [read2],
        polydt,
        ranhex,
        ["AAAA", "TTTT", "CCCC", "GGGG"],
        rt_start=2,
    )
    assert stats == {
        "total": 3,
        "polydT_exact": 1,
        "polydT_corrected": 0,
        "ranhex_exact": 0,
        "ranhex_corrected": 1,
        "unknown_or_ambiguous": 1,
        "assigned": 2,
    }
    assert "@exact_poly/1" in polydt.read_text()
    assert "@corrected_hex/1" in ranhex.read_text()
    assert "unknown" not in polydt.read_text() + ranhex.read_text()


def test_demultiplex_rejects_mismatched_fastq_names(tmp_path):
    read1 = tmp_path / "r1.fastq"
    read2 = tmp_path / "r2.fastq"
    _write_fastq(read1, [("one", 1, "ACGT")])
    _write_fastq(read2, [("two", 2, "NNAAAA")])
    with pytest.raises(ValueError, match="names differ"):
        demultiplex_parse_transcript_reads(
            [read1],
            [read2],
            tmp_path / "poly.fastq",
            tmp_path / "hex.fastq",
            ["AAAA", "TTTT", "CCCC", "GGGG"],
            rt_start=2,
        )


def test_demultiplex_balances_prefix_without_duplicating_reads(tmp_path):
    read1_paths = []
    read2_paths = []
    for file_index in range(2):
        read1 = tmp_path / f"r1_{file_index}.fastq"
        read2 = tmp_path / f"r2_{file_index}.fastq"
        _write_fastq(read1, [
            (f"first_{file_index}", 1, "ACGT"),
            (f"second_{file_index}", 1, "TGCA"),
        ])
        _write_fastq(read2, [
            (f"first_{file_index}", 2, "NNAAAA"),
            (f"second_{file_index}", 2, "NNAAAA"),
        ])
        read1_paths.append(read1)
        read2_paths.append(read2)
    output = tmp_path / "poly.fastq"
    stats = demultiplex_parse_transcript_reads(
        read1_paths,
        read2_paths,
        output,
        tmp_path / "hex.fastq",
        ["AAAA", "TTTT", "CCCC", "GGGG"],
        rt_start=2,
        balanced_prefix_reads=1,
    )
    headers = [
        line.strip() for line in output.read_text().splitlines()
        if line.startswith("@")
    ]
    assert headers == [
        "@first_0/1", "@first_1/1", "@second_0/1", "@second_1/1"
    ]
    assert stats["total"] == 4
