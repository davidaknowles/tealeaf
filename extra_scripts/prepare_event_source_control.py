"""Prepare paired-primer EC input controls with existing Tealeaf machinery.

This diagnostic changes the upstream count source and/or compatibility design,
not the subject split, coverage thresholds or event definitions. It writes a
new cache and its exact feature order, never overwriting production inputs.
"""

import argparse
import hashlib
import json
from pathlib import Path
import pickle

import numpy as np

from extra_scripts.run_differential_splicing import aggregate_inputs, gene_structures
from tealeaf.sc.glm_cv import prepare_paired_primer_glm_data


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--alevin-dir", type=Path, required=True)
    parser.add_argument("--salmon-ref", type=Path, required=True)
    parser.add_argument("--primer-pairs", type=Path, required=True)
    parser.add_argument("--transcript-to-gene", type=Path, required=True)
    parser.add_argument("--barcode-groups", type=Path, required=True)
    parser.add_argument("--ec-design", choices=("binary", "weighted"), default="binary")
    parser.add_argument("--min-half-umis", type=float, default=500)
    parser.add_argument("--min-cells", type=int, default=20)
    parser.add_argument("--min-pseudobulk-umis", type=float, default=100_000)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    if (args.output_dir / "prepared.pkl").exists():
        raise FileExistsError("use a new control directory rather than replacing a prepared cache")
    args.output_dir.mkdir(parents=True, exist_ok=True)
    prepared = prepare_paired_primer_glm_data(args.alevin_dir, args.salmon_ref, args.primer_pairs, ec_design=args.ec_design, regularization_target="theta", min_eq=5, min_half_umis=args.min_half_umis, primer_sampling_model="oligodt_tpm")
    groups, _, _, counts = aggregate_inputs(args, prepared)
    structures = gene_structures(prepared, args.transcript_to_gene)
    with (args.output_dir / "prepared.pkl").open("xb") as handle:
        pickle.dump((groups, counts, *structures), handle, protocol=pickle.HIGHEST_PROTOCOL)
    (args.output_dir / "features.txt").write_text("\n".join(prepared.features) + "\n")
    genes, transcripts, ecs, designs = structures
    global_support = np.logical_or.reduce([np.asarray(design.sum(axis=0)).ravel() > 0 for design in designs])
    manifest = {"settings": {key: str(value) if isinstance(value, Path) else value for key, value in vars(args).items()}, "source": str(args.alevin_dir.resolve()), "features_sha256": hashlib.sha256((args.output_dir / "features.txt").read_bytes()).hexdigest(), "annotation_inputs_sha256": {key: hashlib.sha256(getattr(args, key).read_bytes()).hexdigest() for key in ("primer_pairs", "barcode_groups", "transcript_to_gene")}, "paired_cells": len(prepared.barcodes), "pseudobulks": len(groups), "genes": len(genes), "transcripts": len(prepared.features), "globally_supported_transcripts": int(global_support.sum()), "ECs": [design.shape[0] for design in designs], "gene_EC_median": float(np.median([len(value) for value in ecs])), "genes_below_EC128": int(sum(0 < len(value) <= 128 for value in ecs)), "primer_UMIs": [float(values.sum()) for values in counts], "selection": "original barcode annotations, paired-half thresholds and pooled EC5; no LR or significance selection", "production_changes": False}
    (args.output_dir / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    print(json.dumps(manifest, indent=2), flush=True)


if __name__ == "__main__":
    main()
