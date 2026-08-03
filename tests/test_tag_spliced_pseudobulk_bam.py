import subprocess
import sys

import pandas as pd
import pysam


def alignment(header, name, barcode, umi, cigar=((0, 10), (3, 90), (0, 10))):
    value = pysam.AlignedSegment(header)
    value.query_name = name
    value.query_sequence = "A" * 20
    value.flag = 0
    value.reference_id = 0
    value.reference_start = 100
    value.mapping_quality = 255
    value.cigartuples = cigar
    value.query_qualities = pysam.qualitystring_to_array("I" * 20)
    value.set_tag("NH", 1)
    value.set_tag("CB", barcode)
    value.set_tag("UB", umi)
    return value


def test_tag_spliced_bam_deduplicates_and_sets_read_group(tmp_path):
    input_path = tmp_path / "input.bam"
    header = pysam.AlignmentHeader.from_dict(
        {"HD": {"VN": "1.6"}, "SQ": [{"SN": "chr1", "LN": 1000}]}
    )
    with pysam.AlignmentFile(input_path, "wb", header=header) as handle:
        handle.write(alignment(header, "a", "A_A_A", "TTT"))
        handle.write(alignment(header, "duplicate", "A_A_A", "TTT"))
        handle.write(alignment(header, "unspliced", "AAA", "GGG", ((0, 20),)))
        handle.write(alignment(header, "unlabeled", "CCC", "TTT"))
    metadata = tmp_path / "cells.tsv"
    pd.DataFrame(
        [["run1", "AAA", "type A", "case", "mouse1"]],
        columns=["run", "barcode", "cell_type", "condition", "subject"],
    ).to_csv(metadata, sep="\t", index=False)
    output = tmp_path / "output.bam"
    subprocess.run(
        [
            sys.executable,
            "extra_scripts/tag_spliced_pseudobulk_bam.py",
            "--input",
            str(input_path),
            "--cell-metadata",
            str(metadata),
            "--run",
            "run1",
            "--output",
            str(output),
        ],
        check=True,
    )
    with pysam.AlignmentFile(output, "rb") as handle:
        records = list(handle.fetch(until_eof=True))
        assert [record.get_tag("RG") for record in records] == ["mouse1__type_A"]
        assert handle.header.to_dict()["RG"] == [
            {"ID": "mouse1__type_A", "SM": "mouse1__type_A"}
        ]
