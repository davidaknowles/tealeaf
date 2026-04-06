# tealeaf – Copilot Instructions

## Overview
tealeaf is a Python bioinformatics tool for alternative splicing analysis. It is a modified version of [LeafCutter](https://davidaknowles.github.io/leafcutter/) that quantifies excised introns using isoform abundance and transcriptome annotation rather than raw junction reads.

## Install & Run
```bash
# development install
pip install -e .

# conda environment
conda env create -f requirements.yml

# CLI entry points (after install)
tealeaf-map       # map_gen step
tealeaf-cluster   # bulk/pseudobulk clustering
tealeaf-sc        # single-cell clustering
tealeaf-ggsashimi # sashimi plotting
```

No automated test suite exists in this repository.

## Architecture

The pipeline has three sequential stages:

1. **`tealeaf/map_gen/tealeaf_map_gen.py`** – Parses a GTF annotation (Gencode or Stringtie) to generate isoform→intron and isoform→exon maps. Outputs space-separated TSVs and (for single-cell) sparse matrices in NPZ format. Only needs to run once per annotation file.

2. **`tealeaf/clustering/tealeaf_clustering.py`** – Bulk/pseudobulk mode. Reads Salmon `quant.sf` files, maps isoform TPM/normalized-counts onto introns, builds initial overlap-based clusters, then refines them via `shared_functions.py`. Outputs `{outprefix}refined_cluster` and `{outprefix}ratio_count`.

3. **`tealeaf/sc/tealeaf_sc.py`** – Single-cell mode. Reads alevin-fry output (eq-class matrix), runs an EM step via `sc_utils.py`, generates pseudobulk samples (metacell or bootstrapping), then calls the same clustering logic as step 2.

**`tealeaf/__main__.py`** dispatches CLI entry points using `runpy.run_path()` – each submodule is a self-contained script that uses `if __name__ == '__main__':` guards.

**Shared utilities:**
- `tealeaf/shared_functions.py` – `build_init_cluster`, `process_clusters`, `compute_ratio`; used by both `clustering` and `sc` submodules.
- `tealeaf/utils.py` – `timing_decorator` (wraps functions to print elapsed time) and `write_options_to_file`.

## Key Conventions

- **Argument parsing uses `optparse.OptionParser`**, not `argparse`. All CLI scripts follow this pattern.
- **Space-separated files** (not commas) are the standard intermediate format: `pd.read_csv(..., sep=' ')`.
- **Sparse matrices** (isoform × intron, isoform × exon) are stored as `.npz` files via `scipy.sparse.save_npz`/`load_npz`.
- **`@timing_decorator`** from `tealeaf.utils` should be applied to any major top-level function to log elapsed time.
- **Output file names** follow `{out_prefix}<name>` – the prefix is user-supplied and defaults to `'Leafcutter_'` or `'leafcutter_'` depending on the module.
- **Cluster definition modes**: `1` = overlap only, `2` = overlap + shared intron splice site, `3` = overlap + shared intron splice site + shared exon splice site (default).
- **Normalization scale modes**: `junction` (default), `local`, `global` — controls how Salmon TPM/counts are converted to intron-level values.
- `clustering/tealeaf_clustering.py` imports the shared functions from `tealeaf.shared_functions` — do not duplicate them locally.
- Outputs of the clustering steps are compatible with downstream LeafCutter tools (`leafcutter_ds`, `leafviz`).
