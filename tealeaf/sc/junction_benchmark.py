"""Shared junction-count representation for differential-splicing benchmarks."""

from __future__ import annotations

from dataclasses import dataclass
import json
from pathlib import Path
import re

import numpy as np
import pandas as pd
from scipy import sparse
from scipy.io import mmread


STAR_SJ_COLUMNS = (
    "chromosome",
    "start",
    "end",
    "strand_code",
    "motif",
    "annotated",
    "unique_reads",
    "multi_reads",
    "max_overhang",
)


def benjamini_hochberg(pvalues) -> np.ndarray:
    values = np.asarray(pvalues, dtype=float)
    order = np.argsort(values)
    ranked = values[order]
    adjusted = np.minimum.accumulate(
        (ranked * len(ranked) / np.arange(1, len(ranked) + 1))[::-1]
    )[::-1]
    result = np.empty_like(adjusted)
    result[order] = np.minimum(adjusted, 1.0)
    return result


def _read_lines(path: Path) -> list[str]:
    with open(path) as handle:
        return [line.rstrip("\n") for line in handle]


def normalize_starsolo_barcode(barcode: str) -> str:
    """Remove separators STAR adds between complex-barcode segments."""
    return barcode.replace("_", "")


def _solo_file(run_dir: Path, feature: str, filename: str) -> Path:
    candidates = (
        run_dir / "Solo.out" / feature / "raw" / filename,
        run_dir / feature / "raw" / filename,
    )
    for candidate in candidates:
        if candidate.exists():
            return candidate
    raise FileNotFoundError(f"missing STARsolo {feature} file {filename} in {run_dir}")


def read_starsolo_sj(run_dir: str | Path):
    """Read one STARsolo sparse cell-by-junction count matrix."""
    run_dir = Path(run_dir)
    matrix = sparse.csr_matrix(mmread(_solo_file(run_dir, "SJ", "matrix.mtx"))).T
    barcodes = [
        normalize_starsolo_barcode(value)
        for value in _read_lines(_solo_file(run_dir, "SJ", "barcodes.tsv"))
    ]
    feature_path = _solo_file(run_dir, "SJ", "features.tsv")
    features = pd.read_csv(feature_path, sep="\t", header=None)
    if features.shape[1] == 1:
        sj_path = run_dir / "SJ.out.tab"
        if not sj_path.exists():
            raise ValueError(f"{feature_path} has IDs only and {sj_path} is absent")
        all_junctions = pd.read_csv(sj_path, sep="\t", header=None, names=STAR_SJ_COLUMNS)
        features = all_junctions.iloc[features.iloc[:, 0].astype(int).to_numpy() - 1]
    else:
        features = features.iloc[:, : len(STAR_SJ_COLUMNS)].copy()
        features.columns = STAR_SJ_COLUMNS[: features.shape[1]]
    if matrix.shape != (len(barcodes), len(features)):
        raise ValueError(
            f"STARsolo shape {matrix.shape} disagrees with "
            f"{len(barcodes)} barcodes and {len(features)} junctions"
        )
    features = features.reset_index(drop=True)
    features["strand"] = features["strand_code"].map({0: ".", 1: "+", 2: "-"})
    return matrix, barcodes, features


def prepare_splitseq_whitelists(barcodes: list[str], output_dir: str | Path) -> list[Path]:
    """Write the three segment whitelists for a 24-base Split-seq barcode."""
    clean = sorted({barcode.strip().split(":")[-1] for barcode in barcodes if barcode.strip()})
    invalid = [barcode for barcode in clean if len(barcode) != 24 or re.search("[^ACGT]", barcode)]
    if invalid:
        raise ValueError(f"invalid 24-base Split-seq barcodes, including {invalid[0]!r}")
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    paths = []
    for segment in range(3):
        values = sorted({barcode[8 * segment : 8 * (segment + 1)] for barcode in clean})
        path = output_dir / f"barcode_segment_{segment + 1}.txt"
        path.write_text("".join(f"{value}\n" for value in values))
        paths.append(path)
    return paths


@dataclass
class JunctionBundle:
    counts: sparse.csr_matrix
    junctions: pd.DataFrame
    samples: pd.DataFrame

    def validate(self) -> None:
        if self.counts.shape != (len(self.samples), len(self.junctions)):
            raise ValueError("junction bundle dimensions disagree")
        required = {"sample_id", "cell_type", "condition", "subject", "n_cells"}
        if missing := required - set(self.samples):
            raise ValueError(f"sample metadata is missing {sorted(missing)}")
        if self.samples["sample_id"].duplicated().any():
            raise ValueError("sample IDs are not unique")
        if (self.counts.data < 0).any():
            raise ValueError("junction counts must be nonnegative")

    def save(self, prefix: str | Path) -> None:
        self.validate()
        prefix = Path(prefix)
        prefix.parent.mkdir(parents=True, exist_ok=True)
        sparse.save_npz(prefix.with_name(prefix.name + "_counts.npz"), self.counts)
        self.junctions.to_csv(
            prefix.with_name(prefix.name + "_junctions.tsv.gz"), sep="\t", index=False
        )
        self.samples.to_csv(
            prefix.with_name(prefix.name + "_samples.tsv"), sep="\t", index=False
        )
        manifest = {
            "format": "tealeaf-junction-bundle-v1",
            "samples": len(self.samples),
            "junctions": len(self.junctions),
            "nonzero_counts": int(self.counts.nnz),
            "total_umis": int(self.counts.sum()),
        }
        prefix.with_name(prefix.name + "_manifest.json").write_text(
            json.dumps(manifest, indent=2) + "\n"
        )

    @classmethod
    def load(cls, prefix: str | Path) -> "JunctionBundle":
        prefix = Path(prefix)
        result = cls(
            counts=sparse.load_npz(prefix.with_name(prefix.name + "_counts.npz")),
            junctions=pd.read_csv(
                prefix.with_name(prefix.name + "_junctions.tsv.gz"), sep="\t"
            ),
            samples=pd.read_csv(
                prefix.with_name(prefix.name + "_samples.tsv"), sep="\t", dtype=str
            ),
        )
        result.samples["n_cells"] = result.samples["n_cells"].astype(int)
        result.validate()
        return result


def aggregate_starsolo_runs(
    run_dirs: dict[str, str | Path], cell_metadata: pd.DataFrame
) -> JunctionBundle:
    """Aggregate STARsolo SJ UMI counts into subject-by-cell-type pseudobulks."""
    required = {"run", "barcode", "cell_type", "condition", "subject"}
    if missing := required - set(cell_metadata):
        raise ValueError(f"cell metadata is missing {sorted(missing)}")
    metadata = cell_metadata.copy().astype(str)
    metadata["sample_id"] = (
        metadata["subject"] + "__" + metadata["cell_type"].str.replace(r"[^A-Za-z0-9_.-]", "_", regex=True)
    )
    if metadata[["run", "barcode"]].duplicated().any():
        raise ValueError("run/barcode cell keys are not unique")
    sample_metadata = (
        metadata.groupby("sample_id", sort=True)
        .agg(
            cell_type=("cell_type", "first"),
            condition=("condition", "first"),
            subject=("subject", "first"),
            n_cells=("barcode", "size"),
        )
        .reset_index()
    )
    inconsistent = metadata.groupby("sample_id")[["cell_type", "condition", "subject"]].nunique()
    if (inconsistent > 1).any(axis=None):
        raise ValueError("a pseudobulk sample has inconsistent metadata")
    sample_index = {value: index for index, value in enumerate(sample_metadata.sample_id)}
    junction_index: dict[tuple, int] = {}
    junction_rows: list[dict] = []
    all_rows, all_cols, all_data = [], [], []

    for run, run_dir in run_dirs.items():
        local_metadata = metadata[(metadata["run"] == str(run)) | (metadata["run"] == "*")]
        if local_metadata.empty:
            continue
        counts, barcodes, junctions = read_starsolo_sj(run_dir)
        barcode_index = {barcode: index for index, barcode in enumerate(barcodes)}
        matched = local_metadata[local_metadata["barcode"].isin(barcode_index)].copy()
        if matched.empty:
            continue
        cell_rows = matched["barcode"].map(barcode_index).to_numpy()
        group_cols = matched["sample_id"].map(sample_index).to_numpy()
        assignment = sparse.csr_matrix(
            (np.ones(len(matched), dtype=np.int8), (np.arange(len(matched)), group_cols)),
            shape=(len(matched), len(sample_metadata)),
        )
        pseudobulk = (assignment.T @ counts[cell_rows]).tocoo()
        local_to_global = np.empty(len(junctions), dtype=np.int64)
        for index, row in junctions.iterrows():
            key = (str(row.chromosome), int(row.start), int(row.end), str(row.strand))
            global_index = junction_index.get(key)
            if global_index is None:
                global_index = len(junction_rows)
                junction_index[key] = global_index
                junction_rows.append(row.to_dict())
            local_to_global[index] = global_index
        all_rows.append(pseudobulk.row)
        all_cols.append(local_to_global[pseudobulk.col])
        all_data.append(pseudobulk.data)

    if not all_data:
        raise ValueError("no STARsolo barcodes matched the cell metadata")
    matrix = sparse.coo_matrix(
        (np.concatenate(all_data), (np.concatenate(all_rows), np.concatenate(all_cols))),
        shape=(len(sample_metadata), len(junction_rows)),
    ).tocsr()
    matrix.sum_duplicates()
    nonempty = np.asarray(matrix.sum(axis=1)).ravel() > 0
    return JunctionBundle(
        counts=matrix[nonempty],
        junctions=pd.DataFrame(junction_rows),
        samples=sample_metadata.loc[nonempty].reset_index(drop=True),
    )


def annotate_scquint_groups(bundle: JunctionBundle, gtf_path: str | Path) -> JunctionBundle:
    """Annotate observed junctions and define scQuint three-prime groups once."""
    junctions = bundle.junctions.copy()

    def contig(value):
        value = str(value)
        return value[3:] if value.startswith("chr") else value

    endpoints = set()
    for row in junctions.itertuples():
        chromosome = contig(row.chromosome)
        endpoints.add((chromosome, int(row.start)))
        endpoints.add((chromosome, int(row.end)))
    boundary_genes: dict[tuple[str, int], set[tuple[str, str]]] = {}
    gene_pattern = re.compile(r'gene_id "([^"]+)"')
    name_pattern = re.compile(r'gene_name "([^"]+)"')
    with open(gtf_path) as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 9 or fields[2] != "exon":
                continue
            gene_match = gene_pattern.search(fields[8])
            if gene_match is None:
                continue
            name_match = name_pattern.search(fields[8])
            gene = (gene_match.group(1), name_match.group(1) if name_match else gene_match.group(1))
            chromosome = contig(fields[0])
            for key in ((chromosome, int(fields[3]) - 1), (chromosome, int(fields[4]) + 1)):
                if key in endpoints:
                    boundary_genes.setdefault(key, set()).add(gene)

    gene_ids, gene_names, groups = [], [], []
    for row in junctions.itertuples():
        chromosome = contig(row.chromosome)
        genes = boundary_genes.get((chromosome, int(row.start)), set()) | boundary_genes.get(
            (chromosome, int(row.end)), set()
        )
        if len(genes) != 1 or row.strand not in ("+", "-"):
            gene_ids.append(None)
            gene_names.append(None)
            groups.append(None)
            continue
        gene_id, gene_name = next(iter(genes))
        anchor = int(row.end) if row.strand == "+" else int(row.start)
        gene_ids.append(gene_id)
        gene_names.append(gene_name)
        groups.append(f"{gene_name}_{row.chromosome}_{anchor}_{row.strand}")
    junctions["gene_id"] = gene_ids
    junctions["gene_name"] = gene_names
    junctions["intron_group"] = groups
    sizes = junctions.intron_group.value_counts(dropna=True)
    junctions["intron_group_size"] = junctions.intron_group.map(sizes)
    return JunctionBundle(bundle.counts, junctions, bundle.samples)


def write_leafcutter_junctions(bundle: JunctionBundle, output_dir: str | Path) -> Path:
    """Export pseudobulk UMI junction counts as regtools-compatible BED12 files."""
    bundle.validate()
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    paths = []
    for sample_index, sample in bundle.samples.iterrows():
        path = output_dir / f"{sample.sample_id}.junc"
        row = bundle.counts.getrow(sample_index)
        with open(path, "w") as handle:
            for feature_index, count in zip(row.indices, row.data):
                junction = bundle.junctions.iloc[feature_index]
                start, end = int(junction.start), int(junction.end)
                handle.write(
                    f"{junction.chromosome}\t{start - 1}\t{end}\t.\t{int(count)}\t"
                    f"{junction.strand}\t{start - 1}\t{end}\t0\t2\t1,1\t0,{end-start}\n"
                )
        paths.append(path)
    list_path = output_dir / "junction_files.txt"
    list_path.write_text("".join(f"{path}\n" for path in paths))
    return list_path


def plan_pairwise_contrasts(samples: pd.DataFrame, min_subjects: int = 4) -> list[dict]:
    """Define replicate-aware condition and paired cell-type contrasts."""
    contrasts = []
    for cell_type, local in samples.groupby("cell_type", sort=True):
        conditions = sorted(local.condition.unique())
        for first_index, first in enumerate(conditions):
            for second in conditions[first_index + 1 :]:
                selected = local[local.condition.isin([first, second])]
                counts = selected.groupby("condition").subject.nunique()
                if min(counts.get(first, 0), counts.get(second, 0)) < min_subjects:
                    continue
                contrasts.append(
                    {
                        "contrast_id": f"condition__{cell_type}__{first}__{second}",
                        "effect": "condition",
                        "stratum": str(cell_type),
                        "level_a": str(first),
                        "level_b": str(second),
                        "samples_a": selected.loc[selected.condition == first, "sample_id"].tolist(),
                        "samples_b": selected.loc[selected.condition == second, "sample_id"].tolist(),
                        "paired_subjects": [],
                    }
                )
    cell_types = sorted(samples.cell_type.unique())
    for first_index, first in enumerate(cell_types):
        for second in cell_types[first_index + 1 :]:
            first_subjects = set(samples.loc[samples.cell_type == first, "subject"])
            second_subjects = set(samples.loc[samples.cell_type == second, "subject"])
            shared = sorted(first_subjects & second_subjects)
            if len(shared) < min_subjects:
                continue
            selected = samples[
                samples.cell_type.isin([first, second]) & samples.subject.isin(shared)
            ]
            first_ids = selected.loc[selected.cell_type == first, "sample_id"].tolist()
            second_ids = selected.loc[selected.cell_type == second, "sample_id"].tolist()
            if len(first_ids) != len(shared) or len(second_ids) != len(shared):
                continue
            contrasts.append(
                {
                    "contrast_id": f"cell_type__{first}__{second}",
                    "effect": "cell_type",
                    "stratum": "all",
                    "level_a": str(first),
                    "level_b": str(second),
                    "samples_a": first_ids,
                    "samples_b": second_ids,
                    "paired_subjects": shared,
                }
            )
    return contrasts


def stratified_subject_folds(
    samples: pd.DataFrame, n_folds: int = 2, seed: int = 0
) -> pd.DataFrame:
    """Assign biological subjects to condition-balanced reproducibility folds."""
    subject_conditions = samples.groupby("subject")["condition"].agg(
        lambda values: sorted(set(map(str, values)))
    )
    if subject_conditions.map(len).ne(1).any():
        raise ValueError("each subject must have exactly one condition")
    subjects = pd.DataFrame(
        {
            "subject": subject_conditions.index.astype(str),
            "condition": subject_conditions.map(lambda values: values[0]).to_numpy(),
        }
    )
    if int(n_folds) < 2:
        raise ValueError("n_folds must be at least two")
    rng = np.random.default_rng(seed)
    fold = np.full(len(subjects), -1, dtype=int)
    for _, indices in subjects.groupby("condition", sort=True).groups.items():
        indices = np.asarray(list(indices), dtype=int)
        if len(indices) < int(n_folds):
            raise ValueError("each condition must contain at least n_folds subjects")
        indices = rng.permutation(indices)
        fold[indices] = np.arange(len(indices)) % int(n_folds)
    subjects["fold"] = fold
    return subjects.sort_values(["fold", "condition", "subject"]).reset_index(drop=True)


def normalize_pvalue_table(
    table: pd.DataFrame,
    *,
    method: str,
    contrast: dict,
    feature_column: str,
    pvalue_column: str,
    effect_column: str | None = None,
) -> pd.DataFrame:
    """Map a comparator result onto the common per-contrast schema."""
    result = pd.DataFrame(
        {
            "method": method,
            "contrast_id": contrast["contrast_id"],
            "effect": contrast["effect"],
            "stratum": contrast["stratum"],
            "level_a": contrast["level_a"],
            "level_b": contrast["level_b"],
            "feature_id": table[feature_column].astype(str),
            "p_value": pd.to_numeric(table[pvalue_column], errors="coerce"),
        }
    )
    result["effect_size"] = (
        pd.to_numeric(table[effect_column], errors="coerce")
        if effect_column is not None
        else np.nan
    )
    valid = result.p_value.notna()
    result["q_value"] = np.nan
    if valid.any():
        result.loc[valid, "q_value"] = benjamini_hochberg(result.loc[valid, "p_value"])
    result["significant"] = result.q_value <= 0.05
    result["criterion"] = "BH FDR <= 0.05"
    return result


def simes_omnibus(pairwise: pd.DataFrame) -> pd.DataFrame:
    """Collapse pairwise tests to a feature-level cell-type omnibus score."""
    rows = []
    for (method, feature), local in pairwise.groupby(["method", "feature_id"], sort=False):
        pvalues = np.sort(local.p_value.dropna().to_numpy(dtype=float))
        if not len(pvalues):
            continue
        simes = min(1.0, np.min(pvalues * len(pvalues) / np.arange(1, len(pvalues) + 1)))
        metadata = {}
        for column in ("gene_id", "gene_name"):
            if column in local:
                values = local[column].dropna().astype(str).unique()
                metadata[column] = values[0] if len(values) == 1 else np.nan
        rows.append({
            "method": method,
            "feature_id": feature,
            "p_value": simes,
            "n_pairwise": len(pvalues),
            **metadata,
        })
    result = pd.DataFrame(rows)
    if result.empty:
        result = pd.DataFrame(
            columns=["method", "feature_id", "p_value", "n_pairwise"]
        )
    result["q_value"] = np.nan
    for method, indices in result.groupby("method").groups.items():
        result.loc[indices, "q_value"] = benjamini_hochberg(
            result.loc[indices, "p_value"]
        )
    result["significant"] = result.q_value <= 0.05
    return result
