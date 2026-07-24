"""Parse Biosciences barcode and primer-pair helpers."""

from __future__ import annotations

from contextlib import contextmanager
import gzip
import io
from pathlib import Path
import signal
import subprocess


PRIMER_NAMES = ("polydT", "ranhex")


def read_barcodes(path) -> list[str]:
    with open(path) as handle:
        barcodes = [line.strip() for line in handle if line.strip()]
    if len(barcodes) != len(set(barcodes)):
        raise ValueError(f"duplicate barcodes in {path}")
    return barcodes


def parse_primer_pairs(cell_barcodes, rt_barcodes):
    """Pair poly(dT) and random-hexamer half cells by Parse RT barcode.

    Parse RT lists are ordered with all poly(dT) barcodes followed by the
    corresponding random-hexamer barcodes. Any prefix before the final RT
    barcode, including a batch prefix, is preserved.
    """
    rt_barcodes = list(rt_barcodes)
    if not rt_barcodes or len(rt_barcodes) % 2:
        raise ValueError("the Parse RT list must contain two equal halves")
    if len(rt_barcodes) != len(set(rt_barcodes)):
        raise ValueError("the Parse RT list contains duplicate barcodes")
    half = len(rt_barcodes) // 2
    poly_to_hex = dict(zip(rt_barcodes[:half], rt_barcodes[half:]))
    hex_to_poly = {value: key for key, value in poly_to_hex.items()}
    rt_to_poly = {
        **{barcode: barcode for barcode in poly_to_hex},
        **hex_to_poly,
    }
    rt_lengths = sorted({len(barcode) for barcode in rt_barcodes}, reverse=True)
    known = set(cell_barcodes)
    pairs = {}
    for barcode in cell_barcodes:
        matches = [
            barcode[-length:]
            for length in rt_lengths
            if len(barcode) >= length and barcode[-length:] in rt_to_poly
        ]
        if not matches:
            raise ValueError(f"cell barcode has no known Parse RT suffix: {barcode}")
        if len(matches) > 1:
            raise ValueError(f"cell barcode has ambiguous Parse RT suffix: {barcode}")
        rt = matches[0]
        prefix = barcode[: -len(rt)]
        poly_rt = rt_to_poly[rt]
        polydt = prefix + poly_rt
        ranhex = prefix + poly_to_hex[poly_rt]
        pairs[polydt] = (polydt, ranhex)
    return [
        (cell_id, polydt, ranhex)
        for cell_id, (polydt, ranhex) in sorted(pairs.items())
        if polydt in known or ranhex in known
    ]


def write_primer_pairs(cell_rows, rt_barcode_file, output) -> int:
    pairs = parse_primer_pairs(
        read_barcodes(cell_rows), read_barcodes(rt_barcode_file)
    )
    output = Path(output)
    output.parent.mkdir(parents=True, exist_ok=True)
    with open(output, "w") as handle:
        handle.write("cell_id\tpolydt_barcode\tranhex_barcode\n")
        for row in pairs:
            handle.write("\t".join(row) + "\n")
    return len(pairs)


def build_rt_primer_lookup(rt_barcodes, *, correct_hamming1=True):
    """Map exact and unambiguous one-edit Parse RT barcodes to primer type."""
    rt_barcodes = list(rt_barcodes)
    if not rt_barcodes or len(rt_barcodes) % 2:
        raise ValueError("the Parse RT list must contain two equal halves")
    if len(rt_barcodes) != len(set(rt_barcodes)):
        raise ValueError("the Parse RT list contains duplicate barcodes")
    lengths = {len(barcode) for barcode in rt_barcodes}
    if len(lengths) != 1:
        raise ValueError("FASTQ demultiplexing requires fixed-length RT barcodes")

    half = len(rt_barcodes) // 2
    exact = {
        barcode: PRIMER_NAMES[index >= half]
        for index, barcode in enumerate(rt_barcodes)
    }
    if not correct_hamming1:
        return exact, {}

    candidates = {}
    bases = "ACGT"
    for barcode, primer in exact.items():
        for index, original in enumerate(barcode):
            for base in bases:
                if base == original:
                    continue
                neighbor = barcode[:index] + base + barcode[index + 1:]
                candidates.setdefault(neighbor, set()).add(primer)
    corrected = {
        sequence: next(iter(primers))
        for sequence, primers in candidates.items()
        if sequence not in exact and len(primers) == 1
    }
    return exact, corrected


@contextmanager
def _open_fastq(path, *, pigz=None):
    path = Path(path)
    process = None
    if path.suffix == ".gz" and pigz:
        command = [str(pigz), "-dc"]
        if isinstance(pigz, tuple):
            pigz, threads = pigz
            command = [str(pigz), "-p", str(threads), "-dc"]
        command.append(str(path))
        process = subprocess.Popen(
            command,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
        )
        assert process.stdout is not None
        handle = io.TextIOWrapper(process.stdout)
    elif path.suffix == ".gz":
        handle = gzip.open(path, "rt")
    else:
        handle = path.open()
    try:
        yield handle
    finally:
        handle.close()
        if process is not None:
            stderr = process.stderr.read().decode() if process.stderr else ""
            return_code = process.wait()
            if return_code not in (0, -signal.SIGPIPE):
                raise RuntimeError(
                    f"pigz failed for {path} with status {return_code}: {stderr}"
                )


def _read_fastq_record(handle, path):
    lines = tuple(handle.readline() for _ in range(4))
    if not lines[0]:
        if any(lines[1:]):
            raise ValueError(f"truncated FASTQ record in {path}")
        return None
    if any(not line for line in lines[1:]):
        raise ValueError(f"truncated FASTQ record in {path}")
    if not lines[0].startswith("@") or not lines[2].startswith("+"):
        raise ValueError(f"malformed FASTQ record in {path}")
    if len(lines[1].rstrip("\r\n")) != len(lines[3].rstrip("\r\n")):
        raise ValueError(f"sequence and quality lengths differ in {path}")
    return lines


def _read_name(header):
    return header[1:].split(None, 1)[0].removesuffix("/1").removesuffix("/2")


def demultiplex_parse_transcript_reads(
    read1_paths,
    read2_paths,
    polydt_output,
    ranhex_output,
    rt_barcodes,
    *,
    rt_start=78,
    correct_hamming1=True,
    pigz=None,
    pigz_threads=None,
    max_reads=None,
    balanced_prefix_reads=None,
    progress_callback=None,
):
    """Stream Parse transcript reads into primer-specific FASTQ files.

    Split-seq v2 stores the transcript sequence in read 1 and the RT barcode at
    zero-based ``rt_start`` in read 2. Only read 1 is emitted because the
    resulting streams are inputs to bulk, single-end Salmon quantification.
    """
    read1_paths = [Path(path) for path in read1_paths]
    read2_paths = [Path(path) for path in read2_paths]
    if not read1_paths or len(read1_paths) != len(read2_paths):
        raise ValueError("read1_paths and read2_paths must be nonempty and paired")
    if balanced_prefix_reads is not None and max_reads is not None:
        raise ValueError("balanced_prefix_reads cannot be combined with max_reads")
    if pigz is not None and pigz_threads is not None:
        pigz = (pigz, int(pigz_threads))
    exact, corrected = build_rt_primer_lookup(
        rt_barcodes, correct_hamming1=correct_hamming1
    )
    rt_length = len(next(iter(exact)))
    stats = {
        "total": 0,
        "polydT_exact": 0,
        "polydT_corrected": 0,
        "ranhex_exact": 0,
        "ranhex_corrected": 0,
        "unknown_or_ambiguous": 0,
    }

    with open(polydt_output, "w") as polydt_handle, open(
        ranhex_output, "w"
    ) as ranhex_handle:
        outputs = {"polydT": polydt_handle, "ranhex": ranhex_handle}

        def consume(read1_path, read2_path, *, skip=0, limit=None):
            with _open_fastq(read1_path, pigz=pigz) as read1, _open_fastq(
                read2_path, pigz=pigz
            ) as read2:
                seen = 0
                emitted = 0
                while (
                    (max_reads is None or stats["total"] < int(max_reads))
                    and (limit is None or emitted < int(limit))
                ):
                    record1 = _read_fastq_record(read1, read1_path)
                    record2 = _read_fastq_record(read2, read2_path)
                    if record1 is None or record2 is None:
                        if record1 is not None or record2 is not None:
                            raise ValueError(
                                f"FASTQ pair has unequal records: "
                                f"{read1_path}, {read2_path}"
                            )
                        break
                    if _read_name(record1[0]) != _read_name(record2[0]):
                        raise ValueError(
                            f"FASTQ read names differ: "
                            f"{record1[0].strip()}, {record2[0].strip()}"
                        )
                    seen += 1
                    if seen <= int(skip):
                        continue
                    stats["total"] += 1
                    emitted += 1
                    if (
                        progress_callback is not None
                        and stats["total"] % 1_000_000 == 0
                    ):
                        progress_callback(stats)
                    sequence = record2[1].strip()
                    barcode = sequence[rt_start:rt_start + rt_length]
                    primer = exact.get(barcode)
                    match = "exact"
                    if primer is None:
                        primer = corrected.get(barcode)
                        match = "corrected"
                    if primer is None:
                        stats["unknown_or_ambiguous"] += 1
                        continue
                    outputs[primer].writelines(record1)
                    stats[f"{primer}_{match}"] += 1

        prefix = int(balanced_prefix_reads or 0)
        if prefix:
            for read1_path, read2_path in zip(read1_paths, read2_paths):
                consume(read1_path, read2_path, limit=prefix)
        for read1_path, read2_path in zip(read1_paths, read2_paths):
            consume(read1_path, read2_path, skip=prefix)

    stats["assigned"] = sum(
        stats[f"{primer}_{match}"]
        for primer in PRIMER_NAMES
        for match in ("exact", "corrected")
    )
    return stats
