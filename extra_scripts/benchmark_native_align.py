#!/usr/bin/env python3
"""Read-level comparison of the native target-locus aligner with STARsolo alignments.

The first --max-reads reads of one run are aligned natively
(tealeaf.sc.native_align, reference from native_collect.py --stage reference).
STAR alignments of the same reads (names <run>.<n>) are read from the run's
coordinate-sorted BAM over the target loci. For reads STAR places uniquely
(NH = 1) in a target locus, the script reports how often the native aligner
also reports a unique alignment, with identical aligned blocks and identical
junctions; for native unique alignments, what STAR reported; and for every
read entering a block of --blocks, agreement of its path class (compatibility
mask) between the two aligners. Writes summary.json and per-read tables.
"""

import argparse
import gzip
import json
from pathlib import Path
import time

import numpy as np
import pandas as pd

from tealeaf.sc import native_align as na
from tealeaf.sc.local_path_reads import compatibility_mask, path_structures
from extra_scripts.native_collect import READ_LENGTH, fastq_records, load_reference


def star_alignments(bam, index, reference, run, max_reads):
    """read number -> (NH, genomic blocks, junctions) for primary alignments in the target loci."""
    import pysam

    out = {}
    prefix = len(run) + 1
    with pysam.AlignmentFile(str(bam), "rb", index_filename=str(index)) as handle:
        for chromosome, start, length in zip(reference["locus_chromosome"], reference["locus_genome_start"], reference["locus_ref_end"] - reference["locus_ref_start"]):
            for alignment in handle.fetch(str(chromosome), int(start), int(start + length)):
                if alignment.is_secondary or alignment.is_supplementary or alignment.is_unmapped:
                    continue
                number = int(alignment.query_name[prefix:])
                if number > max_reads:
                    continue
                junctions, position = [], alignment.reference_start
                for operation, size in alignment.cigartuples:
                    if operation == 3:
                        junctions.append((position, position + size))
                    if operation in (0, 2, 3, 7, 8):
                        position += size
                out[number] = (alignment.get_tag("NH"), str(chromosome), tuple(alignment.get_blocks()), tuple(junctions), str(alignment.get_tag("GX")) if alignment.has_tag("GX") else "-")
    return out


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--reference-dir", type=Path, required=True)
    parser.add_argument("--recipe", type=Path, required=True)
    parser.add_argument("--fastq-dir", type=Path, required=True)
    parser.add_argument("--run", required=True)
    parser.add_argument("--max-reads", type=int, default=2_000_000)
    parser.add_argument("--blocks", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--rust-library", type=Path, help="also time the Rust kernel (rust/native_align) on the same reads and check identical output")
    parser.add_argument("--known-junctions", type=Path, help="second pass: junctions.tsv.gz of a first pass (native_collect.py --stage junctions)")
    parser.add_argument("--known-min-reads", type=int, default=20)
    parser.add_argument("--known-min-libraries", type=int, default=4)
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    reference = load_reference(args.reference_dir / "reference.npz", known_junctions=args.known_junctions, min_reads=args.known_min_reads, min_libraries=args.known_min_libraries)
    sequences = []
    for chunk in fastq_records(args.fastq_dir / f"{args.run}_1.fastq.gz"):
        sequences.append(chunk)
        if sum(len(c) for c in sequences) >= args.max_reads:
            break
    reads = np.ascontiguousarray(na.CODE[np.concatenate(sequences)[:args.max_reads]])
    tables = (reference["ref"], reference["index_codes"], reference["index_positions"], reference["buckets"], reference["unique"], na.K, reference["locus_ref_start"], reference["locus_ref_end"], reference["tx_space_start"], reference["tx_ref_start"], reference["tx_length"], reference["exon_offsets"], reference["exon_tx_start"], reference["exon_ref_start"], reference["exon_length"], reference["junctions"])
    clock = time.time()
    na.align_reads(reads[:1000], *tables)
    compile_seconds = time.time() - clock
    clock = time.time()
    numba_out = na.align_reads(reads, *tables)
    seconds = time.time() - clock
    status, orientation, score, matched, mismatched, n_blocks, blocks = numba_out
    kernels = {"threads": int(__import__("os").environ.get("NUMBA_NUM_THREADS", __import__("os").cpu_count())), "numba_first_call_seconds": compile_seconds, "numba_seconds": seconds, "numba_reads_per_second": len(reads) / seconds}
    if args.rust_library:
        aligner = na.rust_aligner(args.rust_library, threads=kernels["threads"])
        clock = time.time()
        rust_out = aligner(reads, *tables)
        kernels["rust_seconds"] = time.time() - clock
        kernels["rust_reads_per_second"] = len(reads) / kernels["rust_seconds"]
        kernels["rust_speedup"] = seconds / kernels["rust_seconds"]
        differing = np.zeros(len(reads), dtype=bool)
        for a, b in zip(numba_out, rust_out):
            differing |= (a != b).reshape(len(reads), -1).any(axis=1)
        kernels["reads_with_different_output"] = int(differing.sum())
    print(kernels, flush=True)
    recipe = json.loads(args.recipe.read_text())
    run = next(run for library in recipe["bams"] for run in library["runs"] if Path(run["path"]).parent.name == args.run)
    clock = time.time()
    star = star_alignments(run["path"], run["index"], reference, args.run, args.max_reads)
    star_seconds = time.time() - clock
    locus = np.searchsorted(reference["locus_ref_start"], blocks[:, 0, 0], side="right") - 1
    native = {}
    for r in np.flatnonzero(status > 0):
        shift = int(reference["locus_genome_start"][locus[r]] - reference["locus_ref_start"][locus[r]])
        segments = tuple((int(blocks[r, i, 0]) + shift, int(blocks[r, i, 1]) + shift) for i in range(n_blocks[r]))
        junctions = tuple((segments[i][1], segments[i + 1][0]) for i in range(len(segments) - 1) if segments[i + 1][0] - segments[i][1] >= na.MIN_INTRON)
        native[int(r) + 1] = (int(status[r]), str(reference["locus_chromosome"][locus[r]]), segments, junctions)
    names = ("unmapped", "unique", "multi", "not_genome_unique", "too_short")
    rows = []
    for number in sorted(set(star) | {n for n, value in native.items() if value[0] == na.UNIQUE}):
        s, v = star.get(number), native.get(number)
        rows.append({"read": number, "star_nh": s[0] if s else 0, "star_gx": s[4] if s else "", "native_status": names[v[0]] if v else "unmapped",
                     "same_blocks": bool(s and v and s[1] == v[1] and tuple(map(tuple, s[2])) == v[2]), "same_junctions": bool(s and v and s[1] == v[1] and s[3] == v[3]),
                     "star_spliced": bool(s and s[3]), "native_spliced": bool(v and v[3]), "star_blocks": json.dumps(s[2]) if s else "", "native_blocks": json.dumps(v[2]) if v else ""})
    table = pd.DataFrame(rows)
    table.to_csv(args.output_dir / "reads.tsv.gz", sep="\t", index=False)
    star_unique = table[table.star_nh == 1]
    native_unique = table[table.native_status == "unique"]
    both = star_unique[star_unique.native_status == "unique"]
    summary = {
        "reads": int(len(reads)), "kernels": kernels, "star_scan_seconds": star_seconds,
        "native_status": dict(zip(names, np.bincount(status, minlength=5).tolist())),
        "star_unique_in_loci": int(len(star_unique)), "star_multi_in_loci": int((table.star_nh > 1).sum()),
        "recall": float(len(both) / max(len(star_unique), 1)),
        "native_status_of_star_unique": star_unique.native_status.value_counts().to_dict(),
        "same_blocks_given_both": float(both.same_blocks.mean()), "same_junctions_given_both": float(both.same_junctions.mean()),
        "spliced_star_unique": int(star_unique.star_spliced.sum()), "spliced_recall": float(((star_unique.star_spliced) & (star_unique.native_status == "unique") & star_unique.same_junctions).sum() / max(star_unique.star_spliced.sum(), 1)),
        "native_unique": int(len(native_unique)), "native_unique_star_nh": {str(k): int(v) for k, v in native_unique.star_nh.value_counts().items()},
    }
    # block classes
    block_set = json.load(gzip.open(args.blocks, "rt"))
    by_chromosome = {}
    for key, block in block_set.items():
        structures = path_structures(block["paths"])
        by_chromosome.setdefault(block["chromosome"], []).append(structures)
    windows = {chromosome: (np.array([x[2][0] for x in values]), np.array([x[2][1] for x in values])) for chromosome, values in by_chromosome.items()}
    agree = total = only_star = only_native = 0
    for number in table.read.tolist():
        s, v = star.get(number), native.get(number)
        s = s if s and s[0] == 1 else None
        v = v if v and v[0] == na.UNIQUE else None
        if s is None and v is None:
            continue
        chromosome = (s or v)[1]
        if chromosome not in windows:
            continue
        spans = [value[2] for value in (s, v) if value is not None and value[1] == chromosome]
        low, high = min(x[0][0] for x in spans), max(x[-1][1] for x in spans)
        starts, ends = windows[chromosome]
        for j in np.flatnonzero((starts < high) & (ends > low)):
            structures = by_chromosome[chromosome][j]
            start, end = structures[2]
            masks = [None if value is None or value[1] != chromosome or not any(a < end and b > start for a, b in value[2]) else compatibility_mask(value[2], value[3], structures) for value in (s, v)]
            if masks == [None, None]:
                continue
            if masks[0] is None:
                only_native += 1
            elif masks[1] is None:
                only_star += 1
            else:
                total += 1
                agree += masks[0] == masks[1]
    summary.update(block_class_pairs=total, block_class_agreement=float(agree / max(total, 1)), block_entries_star_only=only_star, block_entries_native_only=only_native)
    (args.output_dir / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2), flush=True)


if __name__ == "__main__":
    main()
