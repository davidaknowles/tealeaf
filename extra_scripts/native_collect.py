#!/usr/bin/env python3
"""Dependency-free block-local read collection from FASTQ with the native target-locus aligner.

Stages:
  reference  loci and annotated transcripts of the candidate genes, k-mer index
             and genome k-mer uniqueness (tealeaf.sc.native_align), saved as
             reference.npz.
  align      one sequencing run. Read 2 carries the UMI (bases 1-10) and three
             8-base barcode segments (bases 11-18, 49-56, 79-86), each corrected
             to its whitelist within two mismatches when the nearest entry is
             unique; read 1 is aligned. Unique alignments of reads with a
             barcode of the run's library and a UMI without N or homopolymer
             are stored per chunk; junction reads of all unique alignments are
             counted with their largest overhang, as STAR counts SJ.out.tab.
  junctions  junction support per library with STAR's SJ.out.tab filters for
             unannotated junctions (overhang >= 12; intron <= 50/100/200 kb for
             1/2/3 reads), and the block set with novel paths from it
             (tealeaf.sc.local_blocks.annotated_blocks), saved as blocks.json.gz.
  collect    one library: 1MM_All UMI collapsing within cell and locus, block
             classes (AND over a molecule's reads), and counts of all molecules
             and of junction-informative molecules.
  collate    counts, blocks and read opportunities for the all-read and the
             junction-only collections, in the format of
             collect_local_path_reads.py; with --simulated-opportunities the
             opportunities come from error-free reads passed through the
             aligner instead of from path geometry.

With --known-junctions (junctions.tsv.gz of a first pass), canonical-motif
junctions of the first pass that pass --novel-min-reads and
--novel-min-libraries (the support required of novel paths) join the annotated
junctions for rescue of clipped read ends, overhang minimums and scoring, as
in STAR's two-pass mode; align, junctions and collect then use runs_pass2/,
junctions_pass2.tsv.gz and blocks_pass2.json.gz.

collect --sense-only drops reads antisense to the block's gene (read 1 is
sense to the RNA in this protocol). --rust-library runs the Rust alignment
kernel (rust/native_align) in place of numba; both give identical output.
"""

import argparse
from functools import partial
import gzip
import json
from multiprocessing import Pool
import os
from pathlib import Path
import pickle
import subprocess
import time

import numpy as np
import pandas as pd

from tealeaf.sc import native_align as na
from tealeaf.sc.differential import read_gtf_exons
from tealeaf.sc.local_blocks import annotated_blocks, block_genes, block_opportunities, junction_block_opportunities, junction_support

READ_LENGTH = 151
UMI = (0, 10)
SEGMENTS = ((10, 18), (48, 56), (78, 86))
REFERENCE_KEYS = ("ref", "locus_ref_start", "locus_ref_end", "locus_genome_start", "locus_chromosome", "locus_genes", "tx_space_start", "tx_ref_start", "tx_length", "tx_gene", "exon_offsets", "exon_tx_start", "exon_ref_start", "exon_length", "junctions")


def candidate_genes(caches):
    genes = set()
    for cache in caches:
        genes.update(candidate[2] for candidate in pickle.load(open(cache, "rb"))["candidates"])
    return genes


def load_reference(path, index=True, known_junctions=None, min_reads=1, min_libraries=1):
    """Reference arrays; known_junctions (a first-pass junctions.tsv.gz) adds its canonical-motif junctions with at least min_reads reads in min_libraries libraries to reference["junctions"]."""
    data = np.load(path, allow_pickle=False)
    keys = REFERENCE_KEYS + (("index_codes", "index_positions", "buckets", "unique") if index else ())
    reference = {key: data[key] for key in keys}
    if known_junctions is not None:
        table = pd.read_csv(known_junctions, sep="\t")
        table = table[(table["strand"] != ".") & (table["reads"] >= min_reads) & (table["libraries"] >= min_libraries)]
        reference["junctions"] = na.add_junctions(reference, table["chromosome"], table["donor"], table["acceptor"])
    return reference


def suffix(args):
    return "_pass2" if args.known_junctions else ""


def fastq_records(path, chunk_bytes=1 << 27):
    """Yield arrays (n x read length, uint8 bytes) of the sequence lines of a gzipped FASTQ."""
    process = subprocess.Popen(["zcat", str(path)], stdout=subprocess.PIPE, bufsize=1 << 24)
    rest = b""
    try:
        while True:
            data = process.stdout.read(chunk_bytes)
            buffer = rest + data
            if not buffer:
                break
            array = np.frombuffer(buffer, dtype=np.uint8)
            newlines = np.flatnonzero(array == 10)
            n = len(newlines) // 4
            if n == 0:
                if not data:
                    break
                rest = buffer
                continue
            seq_start = newlines[0:4 * n:4] + 1
            seq_end = newlines[1:4 * n:4]
            lengths = seq_end - seq_start
            if (lengths != lengths[0]).any():
                raise ValueError(f"{path}: variable read lengths are not supported")
            yield array[seq_start[:, None] + np.arange(lengths[0])[None, :]]
            rest = buffer[newlines[4 * n - 1] + 1:]
            if not data:
                break
    finally:
        process.kill()
        process.wait()


def paired_batches(read1, read2, batch):
    """Yield (read 1, read 2) sequence arrays of up to batch records, in lockstep."""
    sources = (fastq_records(read1), fastq_records(read2))
    pending = [[], []]
    held = [0, 0]
    exhausted = [False, False]
    while True:
        for i, source in enumerate(sources):
            while not exhausted[i] and held[i] < batch:
                try:
                    piece = next(source)
                    pending[i].append(piece)
                    held[i] += len(piece)
                except StopIteration:
                    exhausted[i] = True
        n = min(batch, held[0], held[1])
        if n == 0:
            if held[0] != held[1]:
                raise ValueError("read 1 and read 2 have different record counts")
            return
        out = []
        for i in range(2):
            merged = np.concatenate(pending[i]) if len(pending[i]) > 1 else pending[i][0]
            out.append(merged[:n])
            pending[i] = [merged[n:]] if len(merged) > n else []
            held[i] -= n
        yield out[0], out[1]


def numpy_batches(read1, read2, batch, tables, lookup, size):
    """Yield (read-1 codes, cell, UMI code, UMI valid) per batch, as tealeaf.sc.native_align.rust_fastq_pairs."""
    homopolymers = [sum(base * 4 ** i for i in range(UMI[1] - UMI[0])) for base in range(4)]
    for first, second in paired_batches(read1, read2, batch):
        codes2 = na.CODE[second]
        segment = [tables[i][na.pack(codes2[:, a:b])[0]] for i, (a, b) in enumerate(SEGMENTS)]
        valid_cell = (segment[0] >= 0) & (segment[1] >= 0) & (segment[2] >= 0)
        cell = np.where(valid_cell, (segment[0] * size + segment[1]) * size + segment[2], 0)
        umi, umi_n = na.pack(codes2[:, UMI[0]:UMI[1]])
        yield na.CODE[first], np.where(valid_cell, lookup[cell], -1), umi, ~umi_n & ~np.isin(umi, homopolymers)


def barcode_tools(whitelist_dir, barcodes):
    """Segment lookup tables and a cell-code -> library barcode index lookup."""
    whitelists = [Path(whitelist_dir, f"barcode_segment_{i}.txt").read_text().split() for i in (1, 2, 3)]
    tables = [na.segment_table(words) for words in whitelists]
    positions = [{word: i for i, word in enumerate(words)} for words in whitelists]
    size = len(whitelists[0])
    lookup = np.full(size ** 3, -1, dtype=np.int64)
    for local, barcode in enumerate(barcodes):
        parts = [barcode[8 * i:8 * i + 8] for i in range(3)]
        lookup[(positions[0][parts[0]] * size + positions[1][parts[1]]) * size + positions[2][parts[2]]] = local
    return tables, lookup, size


def target_genes(args):
    """Target genes: with --all-multiexon-genes every gene of --gtf with a multi-exon transcript, with --all-block-genes every gene with an annotated block of 2..--max-paths paths, otherwise the candidate caches' genes."""
    if args.all_multiexon_genes:
        return {gene for gene, record in read_gtf_exons(args.gtf).items() if any(len(chain) > 1 for chain in record["transcripts"].values())}
    return block_genes(args.block_cache, args.max_paths) if args.all_block_genes else candidate_genes(args.candidate_cache)


def stage_reference(args):
    import pysam

    start = time.time()
    genome = pysam.FastaFile(str(args.genome))
    targets = target_genes(args)
    genes = {gene: record for gene, record in read_gtf_exons(args.gtf).items() if gene in targets}
    reference = na.build_reference(genome, genes, pad=args.pad)
    index = na.build_index(reference)
    built = time.time()
    unique = na.genome_uniqueness(reference, index, genome)
    args.output_dir.mkdir(parents=True, exist_ok=True)
    np.savez(args.output_dir / "reference.npz", **reference, **index, unique=unique)
    summary = {"genes": len(genes), "loci": len(reference["locus_ref_start"]), "transcripts": len(reference["tx_ref_start"]), "reference_bases": int(len(reference["ref"])), "locus_bases": int(reference["tx_space_start"]), "indexed_kmers": int(len(index["index_codes"])), "genome_unique_fraction": float(unique[index["index_positions"]].mean()), "index_seconds": built - start, "uniqueness_seconds": time.time() - built}
    (args.output_dir / "reference.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(summary, flush=True)


def stage_align(args):
    recipe = json.loads(args.recipe.read_text())
    library = next(library for library in recipe["bams"] if any(Path(run["path"]).parent.name == args.run for run in library["runs"]))
    barcodes = sorted(library["barcode_groups"])
    tables, lookup, size = barcode_tools(args.whitelist_dir, barcodes)
    reference = load_reference(args.output_dir / "reference.npz", known_junctions=args.known_junctions, min_reads=args.novel_min_reads, min_libraries=args.novel_min_libraries)
    kernel = na.rust_aligner(args.rust_library) if args.rust_library else na.align_reads
    out = args.output_dir / f"runs{suffix(args)}" / args.run
    out.mkdir(parents=True, exist_ok=True)
    status_counts = np.zeros(5, dtype=np.int64)
    junction_keys, junction_overhangs = [], []
    stats = {"reads": 0, "library_barcode": 0, "valid_umi": 0, "stored": 0}
    timing = {"read": 0., "align": 0., "total": time.time()}
    read1, read2 = (args.fastq_dir / f"{args.run}_{mate}.fastq.gz" for mate in (1, 2))
    if args.rust_library and not args.numpy_reader:
        batches = na.rust_fastq_pairs(args.rust_library, read1, read2, tables, lookup, size, UMI, SEGMENTS, args.chunk_reads)
    else:
        batches = numpy_batches(read1, read2, args.chunk_reads, tables, lookup, size)
    chunk = 0
    while True:
        clock = time.time()
        batch = next(batches, None)
        timing["read"] += time.time() - clock
        if batch is None or (args.max_reads and stats["reads"] >= args.max_reads):
            break
        codes1, local, umi, valid_umi = (part[:args.max_reads - stats["reads"]] for part in batch) if args.max_reads else batch
        clock = time.time()
        status, orientation, _, _, _, n_blocks, blocks = kernel(np.ascontiguousarray(codes1), reference["ref"], reference["index_codes"], reference["index_positions"], reference["buckets"], reference["unique"], na.K, reference["locus_ref_start"], reference["locus_ref_end"], reference["tx_space_start"], reference["tx_ref_start"], reference["tx_length"], reference["exon_offsets"], reference["exon_tx_start"], reference["exon_ref_start"], reference["exon_length"], reference["junctions"])
        timing["align"] += time.time() - clock
        status_counts += np.bincount(status, minlength=5)
        unique = status == na.UNIQUE
        for j in range(na.MAX_BLOCKS - 1):
            has = unique & (n_blocks > j + 1)
            gap = blocks[:, j + 1, 0] - blocks[:, j, 1]
            spliced = has & (gap >= na.MIN_INTRON)
            junction_keys.append((blocks[spliced, j, 1] << 32) | blocks[spliced, j + 1, 0])
            junction_overhangs.append(np.minimum(blocks[spliced, j, 1] - blocks[spliced, j, 0], blocks[spliced, j + 1, 1] - blocks[spliced, j + 1, 0]))
        keep = unique & (local >= 0) & valid_umi
        offsets, starts, ends = na.flatten_blocks(n_blocks[keep], blocks[keep])
        locus = np.searchsorted(reference["locus_ref_start"], blocks[keep, 0, 0], side="right") - 1
        np.savez(out / f"chunk_{chunk:04d}.npz", read=np.flatnonzero(keep) + stats["reads"], cell=local[keep].astype(np.int32), umi=umi[keep].astype(np.int32), locus=locus.astype(np.int32), orientation=orientation[keep], offsets=offsets, starts=starts.astype(np.int32), ends=ends.astype(np.int32))
        stats["reads"] += len(codes1)
        stats["library_barcode"] += int((local >= 0).sum())
        stats["valid_umi"] += int(((local >= 0) & valid_umi).sum())
        stats["stored"] += int(keep.sum())
        print(f"chunk {chunk}: {stats['reads']} reads, {dict(zip(('unmapped', 'unique', 'multi', 'not_genome_unique', 'too_short'), status_counts.tolist()))}, read {timing['read']:.0f} s, align {timing['align']:.0f} s, total {time.time() - timing['total']:.0f} s", flush=True)
        chunk += 1
    keys = np.concatenate(junction_keys) if junction_keys else np.zeros(0, dtype=np.int64)
    overhangs = np.concatenate(junction_overhangs) if junction_overhangs else np.zeros(0, dtype=np.int64)
    order = np.argsort(keys, kind="stable")
    unique_keys, first_index, reads = np.unique(keys[order], return_index=True, return_counts=True)
    max_overhang = np.maximum.reduceat(overhangs[order], first_index) if len(order) else np.zeros(0, dtype=np.int64)
    np.savez(out / "junctions.npz", keys=unique_keys, reads=reads, max_overhang=max_overhang)
    timing["total"] = time.time() - timing["total"]
    stats.update(status=dict(zip(("unmapped", "unique", "multi", "not_genome_unique", "too_short"), status_counts.tolist())), seconds=timing, library=library["library"], threads=int(os.environ.get("NUMBA_NUM_THREADS", os.cpu_count())))
    (out / "stats.json").write_text(json.dumps(stats, indent=2) + "\n")
    print(stats, flush=True)


def run_junctions(path, reference, annotated):
    """Genomic junction support of one run after STAR's SJ.out.tab filters (unannotated junctions only)."""
    data = np.load(path)
    keys, reads, overhang = data["keys"], data["reads"], data["max_overhang"]
    donor, acceptor = keys >> 32, keys & ((1 << 32) - 1)
    intron = acceptor - donor
    is_annotated = np.isin(keys, annotated)
    passes = (overhang >= 12) & ((reads >= 4) | ((reads == 3) & (intron <= 200000)) | ((reads == 2) & (intron <= 100000)) | ((reads == 1) & (intron <= 50000)))
    keep = is_annotated | passes
    donor, acceptor, reads = donor[keep], acceptor[keep], reads[keep]
    locus = np.searchsorted(reference["locus_ref_start"], donor, side="right") - 1
    shift = reference["locus_genome_start"][locus] - reference["locus_ref_start"][locus]
    motif = na.junction_motifs(reference["ref"], donor, acceptor)
    strand = np.where(motif == 0, ".", np.where(motif % 2 == 1, "+", "-"))
    chromosome = reference["locus_chromosome"][locus].astype(str)
    return {(c, d, a): (n, s) for c, d, a, n, s in zip(chromosome.tolist(), (donor + shift).tolist(), (acceptor + shift).tolist(), reads.tolist(), strand.tolist())}


def stage_junctions(args):
    recipe = json.loads(args.recipe.read_text())
    reference = load_reference(args.output_dir / "reference.npz", index=False)
    reads, libraries, strands = {}, {}, {}
    for index, library in enumerate(recipe["bams"]):
        for run in library["runs"]:
            for key, (n, strand) in run_junctions(args.output_dir / f"runs{suffix(args)}" / Path(run["path"]).parent.name / "junctions.npz", reference, reference["junctions"]).items():
                reads[key] = reads.get(key, 0) + n
                libraries.setdefault(key, set()).add(index)
                strands[key] = strand
    support = {key: (reads[key], len(libraries[key]), strands[key]) for key in reads}
    pd.DataFrame([(c, d, a, *value) for (c, d, a), value in sorted(support.items())], columns=["chromosome", "donor", "acceptor", "reads", "libraries", "strand"]).to_csv(args.output_dir / f"junctions{suffix(args)}.tsv.gz", sep="\t", index=False)
    genes = read_gtf_exons(args.gtf)
    options = {"min_reads": args.novel_min_reads, "min_libraries": args.novel_min_libraries, "max_novel": args.novel_max, "terminal_length": args.novel_terminal_length}
    blocks = annotated_blocks(args.candidate_cache, args.block_cache, args.max_paths, junction_support(support), options, genes, genes=target_genes(args))
    with gzip.open(args.blocks, "wt") as handle:
        json.dump(blocks, handle)
    novel = sum(len(block["signatures"]) - block["n_annotated"] for block in blocks.values())
    print(f"{len(support)} junctions; {len(blocks)} blocks, {sum(':N' in block['block_id'] for block in blocks.values())} new, {novel} novel paths", flush=True)


def stage_collect(args):
    recipe = json.loads(args.recipe.read_text())
    library = recipe["bams"][args.library_index]
    barcodes = sorted(library["barcode_groups"])
    groups = [tuple(library["barcode_groups"][barcode]) for barcode in barcodes]
    reference = load_reference(args.output_dir / "reference.npz", index=False)
    blocks = json.load(gzip.open(args.blocks, "rt"))
    keys, tables = na.block_tables(blocks, reference)
    parts = {name: [] for name in ("cell", "umi", "locus", "orientation", "offsets", "starts", "ends")}
    total = 0
    for run in library["runs"]:
        for path in sorted((args.output_dir / f"runs{suffix(args)}" / Path(run["path"]).parent.name).glob("chunk_*.npz")):
            data = np.load(path)
            for name in ("cell", "umi", "locus", "orientation", "starts", "ends"):
                parts[name].append(data[name].astype(np.int64))
            parts["offsets"].append(data["offsets"][:-1] + total)
            total += int(data["offsets"][-1])
    cell, umi, locus, orientation = (np.concatenate(parts[name]) for name in ("cell", "umi", "locus", "orientation"))
    offsets = np.r_[np.concatenate(parts["offsets"]), total].astype(np.int64)
    starts, ends = np.concatenate(parts["starts"]), np.concatenate(parts["ends"])
    if cell.max(initial=0) >= 1 << 19 or locus.max(initial=0) >= 1 << 24 or len(keys) >= 1 << 24:
        raise ValueError("cell, locus or block index exceeds the molecule key fields")
    triple = (cell << 44) | (locus << 20) | umi
    unique, inverse = np.unique(triple, return_inverse=True)
    representative = na.collapse_umis(unique >> 20, unique & ((1 << 20) - 1), UMI[1] - UMI[0])
    umi = representative[inverse]
    read, block, mask, informative = na.block_masks(locus, offsets, starts, ends, tables["locus"], tables["window_start"], tables["window_end"], tables["order"], tables["path_offsets"], tables["exon_offsets"], tables["exon_start"], tables["exon_end"], tables["junction_offsets"], tables["junction_keys"], tables["variable_offsets"], tables["variable_keys"])
    antisense = (orientation[read] == 0) != (tables["strand"][block] == 1)
    if args.sense_only:
        read, block, mask, informative = read[~antisense], block[~antisense], mask[~antisense], informative[~antisense]
    if args.drop_unfit_reads:
        fit = mask != 0
        read, block, mask, informative = read[fit], block[fit], mask[fit], informative[fit]
    molecule = (cell[read] << 44) | (umi[read] << 24) | block
    order = np.argsort(molecule, kind="stable")
    molecule, mask, informative = molecule[order], mask[order], informative[order]
    first = np.r_[0, np.flatnonzero(np.diff(molecule)) + 1]
    molecule_mask = np.bitwise_and.reduceat(mask, first)
    molecule_informative = np.logical_or.reduceat(informative, first)
    molecule_cell = molecule[first] >> 44
    molecule_block = molecule[first] & ((1 << 24) - 1)
    for name, selected in (("all", np.ones(len(first), dtype=bool)), ("junction", molecule_informative)):
        frame = pd.DataFrame({"block": molecule_block[selected], "cell": molecule_cell[selected], "mask": molecule_mask[selected]})
        filters = {"molecules": int(len(frame)), "incompatible_molecules": int((frame["mask"] == 0).sum()), "reads": int(len(cell)), "blocks": len(keys), "antisense_read_entries": int(antisense.sum()), "sense_only": bool(args.sense_only)}
        if name == "junction":
            filters["no_variable_junction"] = int((~molecule_informative).sum())
        frame = frame[frame["mask"] != 0]
        counts = frame.groupby(["block", "cell", "mask"]).size().reset_index(name="count")
        group = np.array(groups, dtype=object)[counts["cell"].to_numpy()]
        rows = pd.DataFrame({"path_key": np.array(keys, dtype=object)[counts["block"].to_numpy()], "subject": [g[0] for g in group], "cell_type": [g[1] for g in group], "primer": [g[2] for g in group], "mask": counts["mask"].to_numpy(), "count": counts["count"].to_numpy()})
        rows = rows.groupby(["path_key", "subject", "cell_type", "primer", "mask"], as_index=False)["count"].sum()
        out = args.collection_dir / name / f"part_{args.library_index}"
        out.mkdir(parents=True, exist_ok=True)
        rows.to_csv(out / "counts.tsv.gz", sep="\t", index=False)
        (out / "filters.json").write_text(json.dumps(filters, indent=2) + "\n")
        print(f"library {library['library']} {name}: {filters}", flush=True)


def stage_collate(args):
    recipe = json.loads(args.recipe.read_text())
    blocks = json.load(gzip.open(args.blocks, "rt"))
    simulated = None
    if args.simulated_opportunities:
        reference = load_reference(args.output_dir / "reference.npz", known_junctions=args.known_junctions, min_reads=args.novel_min_reads, min_libraries=args.novel_min_libraries)
        tables = [reference[name] for name in ("ref", "index_codes", "index_positions", "buckets", "unique")] + [na.K] + [reference[name] for name in ("locus_ref_start", "locus_ref_end", "tx_space_start", "tx_ref_start", "tx_length", "exon_offsets", "exon_tx_start", "exon_ref_start", "exon_length", "junctions")]
        kernel = na.rust_aligner(args.rust_library) if args.rust_library else na.align_reads
        clock = time.time()
        simulated = na.simulated_opportunities(blocks, reference, lambda reads: kernel(reads, *tables), args.read_length)
        print(f"simulated opportunities for {len(simulated)} blocks in {time.time() - clock:.0f} s", flush=True)
    for variant, (name, function) in enumerate((("all", block_opportunities), ("junction", junction_block_opportunities))):
        source = args.collection_dir / name
        out = (args.collate_dir or args.collection_dir) / name
        out.mkdir(parents=True, exist_ok=True)
        parts = sorted(source.glob("part_*/counts.tsv.gz"))
        if len(parts) != len(recipe["bams"]):
            raise ValueError(f"{name}: expected {len(recipe['bams'])} parts, found {len(parts)}")
        table = pd.concat([pd.read_csv(path, sep="\t") for path in parts], ignore_index=True)
        table = table.groupby(["path_key", "subject", "cell_type", "primer", "mask"], as_index=False)["count"].sum()
        table.to_csv(out / "counts.tsv.gz", sep="\t", index=False)
        with gzip.open(out / "blocks.json.gz", "wt") as handle:
            json.dump(blocks, handle)
        if simulated is not None:
            opportunities = {key: {str(mask): vector.tolist() for mask, vector in simulated[key][variant].items()} for key in blocks if key in simulated}
        else:
            with Pool(len(os.sched_getaffinity(0))) as pool:
                opportunities = dict(zip(blocks, pool.map(partial(function, read_length=args.read_length), list(blocks.values()), chunksize=8)))
        with gzip.open(out / "opportunities.json.gz", "wt") as handle:
            json.dump(opportunities, handle)
        filters = {}
        for path in sorted(source.glob("part_*/filters.json")):
            for key, value in json.loads(path.read_text()).items():
                filters[key] = filters.get(key, 0) + value
        (out / "filters.json").write_text(json.dumps(filters, indent=2) + "\n")
        print(f"{name}: {len(blocks)} blocks, {len(table)} count rows, {int(table['count'].sum())} molecules; {filters}", flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--stage", choices=("reference", "align", "junctions", "collect", "collate"), required=True)
    parser.add_argument("--output-dir", type=Path, required=True, help="reference, aligned runs and junction support")
    parser.add_argument("--recipe", type=Path)
    parser.add_argument("--candidate-cache", type=Path, action="append")
    parser.add_argument("--block-cache", type=Path)
    parser.add_argument("--gtf", type=Path)
    parser.add_argument("--genome", type=Path)
    parser.add_argument("--pad", type=int, default=0, help="reference: bases added to each side of a target gene span (1 kb changed under 0.2%% of block molecules)")
    parser.add_argument("--fastq-dir", type=Path)
    parser.add_argument("--whitelist-dir", type=Path)
    parser.add_argument("--run")
    parser.add_argument("--chunk-reads", type=int, default=4_000_000)
    parser.add_argument("--max-reads", type=int, default=0)
    parser.add_argument("--library-index", type=int, default=0)
    parser.add_argument("--blocks", type=Path, help="blocks.json.gz written by the junctions stage and classified into by collect and collate (default: <output-dir>/blocks.json.gz, blocks_pass2.json.gz with --known-junctions)")
    parser.add_argument("--collection-dir", type=Path, help="collection root; all/ and junction/ hold the two read-class collections")
    parser.add_argument("--sense-only", action="store_true", help="collect: drop reads antisense to the block's gene")
    parser.add_argument("--drop-unfit-reads", action="store_true", help="collect: leave reads compatible with no path out of the molecule's class (AND) and junction flag instead of discarding the molecule")
    parser.add_argument("--simulated-opportunities", action="store_true", help="collate: read opportunities from error-free reads passed through the aligner (tealeaf.sc.native_align.simulated_opportunities)")
    parser.add_argument("--rust-library", type=Path, help="use the Rust alignment kernel and FASTQ reader (rust/native_align) instead of numba and numpy")
    parser.add_argument("--numpy-reader", action="store_true", help="align: read FASTQ with zcat and numpy even with --rust-library")
    parser.add_argument("--known-junctions", type=Path, help="second pass: junctions.tsv.gz of the first pass; junctions passing --novel-min-reads and --novel-min-libraries become known junctions")
    parser.add_argument("--collate-dir", type=Path, help="collate: write the collated collections here instead of --collection-dir")
    parser.add_argument("--read-length", type=int, default=READ_LENGTH, help="collate: read 1 length for read opportunities (151 for the microglia data)")
    parser.add_argument("--max-paths", type=int, default=30)
    parser.add_argument("--novel-min-reads", type=int, default=20)
    parser.add_argument("--novel-min-libraries", type=int, default=4)
    parser.add_argument("--novel-max", type=int, default=8)
    parser.add_argument("--all-multiexon-genes", action="store_true", help="reference and junctions: target every gene of --gtf with a multi-exon transcript (genes without an annotated multi-path block can still get new blocks)")
    parser.add_argument("--all-block-genes", action="store_true", help="reference and junctions: target every gene with an annotated block of 2..--max-paths paths (--block-cache) instead of the candidate caches' genes")
    parser.add_argument("--novel-terminal-length", type=int, default=0, help="length of unannotated first or last exons from deep intronic novel splice sites (tealeaf.sc.novel_paths); 0 skips them")
    args = parser.parse_args()
    if args.blocks is None:
        args.blocks = args.output_dir / f"blocks{suffix(args)}.json.gz"
    {"reference": stage_reference, "align": stage_align, "junctions": stage_junctions, "collect": stage_collect, "collate": stage_collate}[args.stage](args)


if __name__ == "__main__":
    main()
