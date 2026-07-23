"""Parse Biosciences barcode and primer-pair helpers."""

from __future__ import annotations

from pathlib import Path


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
    half = len(rt_barcodes) // 2
    poly_to_hex = dict(zip(rt_barcodes[:half], rt_barcodes[half:]))
    hex_to_poly = {value: key for key, value in poly_to_hex.items()}
    known = set(cell_barcodes)
    pairs = {}
    for barcode in cell_barcodes:
        matched = False
        for rt in rt_barcodes:
            if not barcode.endswith(rt):
                continue
            prefix = barcode[: -len(rt)]
            poly_rt = rt if rt in poly_to_hex else hex_to_poly[rt]
            polydt = prefix + poly_rt
            ranhex = prefix + poly_to_hex[poly_rt]
            cell_id = polydt
            pairs[cell_id] = (polydt, ranhex)
            matched = True
            break
        if not matched:
            raise ValueError(f"cell barcode has no known Parse RT suffix: {barcode}")
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
