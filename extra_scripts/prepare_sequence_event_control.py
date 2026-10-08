"""Prepare isolated reference-opportunity inputs without changing count support."""

import argparse
import hashlib
import json
from pathlib import Path
import pickle

import numpy as np
import pandas as pd

from extra_scripts.run_suppa2_tealeaf_hybrid import canonical, supported_gene_transcripts
from tealeaf.sc.sequence_ec import replace_ec_kernel_rows_inplace


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-root", type=Path, required=True)
    parser.add_argument("--kernel-root", type=Path, required=True)
    parser.add_argument("--source", choices=("parsimony_binary", "original_binary"), required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--shard-count", type=int, default=16)
    args = parser.parse_args()
    if args.output_dir.exists():
        raise FileExistsError("use a new experimental directory, never overwrite prepared inputs")
    source = args.input_root / f"{args.source}_paired"
    shards = [args.kernel_root / args.source / f"shard_{index}" for index in range(args.shard_count)]
    settings = [json.loads((path / "settings.json").read_text()) for path in shards]
    recipe = settings[0]
    if any(value != recipe for value in settings) or recipe["source"] != args.source or recipe["shard_count"] != args.shard_count or recipe["read_lengths"] != [150] or recipe["terminal_start_windows"] != ["full"] or recipe["background_fraction"] != .01 or not recipe["selection"].startswith("whole screened gene catalog"):
        raise ValueError("requires complete fixed L150/full-window/background .01 catalog")
    frames = [pd.read_csv(path / "genes.tsv.gz", sep="\t") for path in shards]
    records = pd.concat(frames, ignore_index=True)
    expected = {(gene, primer) for gene in recipe["requested_genes"] for primer in (0, 1)}
    if len(records) != len(expected) or set(map(tuple, records[["gene_id", "primer"]].values)) != expected or not records.source.eq(args.source).all() or not records.status.eq("ok").all() or not records.full_positive_support.astype(str).str.lower().eq("true").all() or not records.retained_molecule_fraction.eq(1).all():
        raise ValueError("missing, duplicate, failed or support-losing kernel family")
    with (source / "prepared.pkl").open("rb") as handle:
        metadata, counts, genes, gene_tx, gene_ecs, original = pickle.load(handle)
    designs = tuple(value.astype(float).tocsr(copy=True) for value in original)
    for value in designs:
        value.sort_indices()
    features_bytes = (source / "features.txt").read_bytes()
    features = features_bytes.decode().splitlines()
    lookup = {canonical(value): index for index, value in enumerate(genes)}
    if len(lookup) != len(genes):
        raise ValueError("ambiguous canonical gene identifiers")
    touched = set()
    hashes = {}
    for shard_index, path in enumerate(shards):
        for gene_id in recipe["requested_genes"][shard_index::args.shard_count]:
            gene = lookup[canonical(gene_id)]
            transcripts = supported_gene_transcripts(gene, gene_tx, gene_ecs, designs)
            packet = path / f"{gene_id}_L150_Wfull_E0.01.npz"
            with np.load(packet, allow_pickle=False) as kernel:
                rows = np.asarray(gene_ecs[gene])
                if not np.array_equal(kernel["ec_indices"], rows) or not np.array_equal(kernel["transcripts"], [features[index] for index in transcripts]) or float(kernel["background_fraction"]) != .01:
                    raise ValueError("kernel EC/transcript alignment or recipe mismatch")
                if touched.intersection(map(int, rows)):
                    raise ValueError("screened gene EC rows overlap")
                replace_ec_kernel_rows_inplace(designs, rows, np.asarray(transcripts), (kernel["dt"], kernel["rh"]))
                touched.update(map(int, rows))
            hashes[gene_id] = hashlib.sha256(packet.read_bytes()).hexdigest()
    for before, after in zip(original, designs):
        before = before.tocsr()
        if not np.array_equal(before.indptr, after.indptr) or not np.array_equal(before.indices, after.indices) or not (after.data > 0).all():
            raise ValueError("global EC compatibility structure changed")
    args.output_dir.mkdir(parents=True)
    with (args.output_dir / "prepared.pkl").open("xb") as handle:
        pickle.dump((metadata, counts, genes, gene_tx, gene_ecs, designs), handle, protocol=pickle.HIGHEST_PROTOCOL)
    (args.output_dir / "features.txt").write_bytes(features_bytes)
    manifest = dict(source=args.source, source_prepared=str((source / "prepared.pkl").resolve()), kernel_recipe=recipe, kernel_sha256=hashes, patched_genes=len(hashes), patched_EC_rows=len(touched), features_sha256=hashlib.sha256(features_bytes).hexdigest(), count_objects="original unmodified count matrices, metadata, gene structures and transcript ordering", support="every global CSR index/indptr and every positive source compatibility preserved", primer_UMIs=[float(value.sum()) for value in counts], interpretation="experimental exact gene-local single-read opportunity maps with fixed background, not real aligner/UMI or validated sampling law", production_changes=False)
    (args.output_dir / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    print(f"{args.source}, {len(hashes)} genes, {len(touched)} EC rows, unchanged counts and compatibility support", flush=True)


if __name__ == "__main__":
    main()
