"""PTPRC (CD45) exon 4-6 choice per cell type in GSE320507: the junction leaving exon 3 goes to exon 4 (RA), 5 (RB),
6 (RC) or 7 (RO). Shares from STARsolo junction counts pooled over each cell type's pseudobulks (short reads) and from
Kinnex long-read molecules with that junction (long reads, CellTypist labels as in the long-read assessment).
Coordinates are the exon ends and starts of GENCODE v49 PTPRC (0-based starts, 1-based ends); STAR introns are 1-based."""
import sys
from pathlib import Path
import numpy as np, pandas as pd
from scipy import sparse
sys.path.insert(0, ".")
from tealeaf.sc.longread_molecules import cell_labels, parse_chain

B = Path.home() / "tealeaf_runs/pbmc/GSE320507/junction_benchmark"
EXON3_END = 198692373
TARGETS = {"exon 4 (RA)": 198696711, "exon 5 (RB)": 198699563, "exon 6 (RC)": 198702386, "exon 7 (RO)": 198703297}  # 0-based starts
TYPES = ["Classical_monocytes", "Non-classical_monocytes", "Naive_B_cells", "Memory_B_cells", "Tcm/Naive_helper_T_cells", "Regulatory_T_cells", "Tem/Trm_cytotoxic_T_cells", "Tem/Temra_cytotoxic_T_cells", "MAIT_cells", "NK_cells", "CD16+_NK_cells"]
junctions = pd.read_csv(B / "pseudobulk_junctions_junctions.tsv.gz", sep="\t", usecols=["chromosome", "start", "end"])
counts = sparse.load_npz(B / "pseudobulk_junctions_counts.npz").tocsc()
samples = pd.read_csv(B / "pseudobulk_junctions_samples.tsv", sep="\t")
short = {}
for name, start0 in TARGETS.items():
    hit = np.flatnonzero((junctions.chromosome.eq("chr1") & junctions.start.eq(EXON3_END + 1) & junctions.end.eq(start0)).to_numpy())
    col = counts[:, hit[0]].toarray().ravel() if len(hit) else np.zeros(len(samples))
    short[name] = pd.Series(col).groupby(samples.cell_type.to_numpy()).sum()
short = pd.DataFrame(short)
molecules = pd.read_parquet(".cache/pbmc/kinnex_parse/molecules/molecules.parquet", columns=["cell", "gene", "chain"])
molecules = molecules.loc[molecules.gene.astype(str).str.startswith("ENSG00000081237")]
labels = cell_labels(".cache/pbmc/kinnex_parse/celltypist_labels.tsv.gz", "Immune_All_Low_majority")
molecules = molecules.loc[molecules.cell.isin(labels.index)].assign(cell_type=lambda d: d.cell.map(labels.cell_type))
rows = []
for chain, ct in zip(molecules.chain, molecules.cell_type):
    exons = parse_chain(chain)
    for (s0, e0), (s1, e1) in zip(exons, exons[1:]):
        if e0 == EXON3_END:
            target = next((n for n, s in TARGETS.items() if s1 == s), None)
            if target:
                rows.append((ct, target))
long = pd.DataFrame(rows, columns=["cell_type", "target"]).groupby(["cell_type", "target"]).size().unstack(fill_value=0).reindex(columns=list(TARGETS), fill_value=0)
pd.set_option("display.width", 250)
for name, table in (("short reads (junction counts)", short), ("long reads (molecules)", long)):
    t = table.reindex(TYPES).dropna(how="all")
    shares = t.div(t.sum(axis=1), axis=0).round(3)
    shares["reads"] = t.sum(axis=1).astype(int)
    print(f"== {name}: share of exon 3 junctions to each exon"); print(shares.to_string())
out = pd.concat({"short": short.reindex(TYPES), "long": long.reindex(TYPES)}, names=["source"])
out.to_csv("analyses/ptprc/exon3_choice.tsv", sep="\t")
