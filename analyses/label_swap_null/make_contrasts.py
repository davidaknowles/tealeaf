"""Within-subject label-swap null on GSE320507: the 21 cell-type pairs among the seven most abundant cell types,
each subject's two labels exchanged with probability 1/2 (permute_paired_contrasts). Writes the unpermuted and
permuted contrasts, the level pairs, and the swapped subjects per pair for run_local_path_tests.py --swapped-pairs."""
import itertools, json, sys
from pathlib import Path
import pandas as pd
sys.path.insert(0, ".")
from tealeaf.sc.junction_benchmark import permute_paired_contrasts
B = Path.home() / "tealeaf_runs/pbmc/GSE320507/junction_benchmark"; O = Path(".cache/label_swap")
samples = pd.read_csv(B / "pseudobulk_junctions_samples.tsv", sep="\t", dtype={"subject": str})
top = samples.groupby("cell_type").n_cells.sum().sort_values(ascending=False).index[:7]
pairs = {tuple(sorted(p)) for p in itertools.combinations(top, 2)}
contrasts = [c for c in json.loads((B / "contrasts.json").read_text()) if tuple(sorted((c["level_a"], c["level_b"]))) in pairs]
assert len(contrasts) == len(pairs), len(contrasts)
(O / "contrasts.json").write_text(json.dumps(contrasts, indent=2) + "\n")
permuted = permute_paired_contrasts(contrasts, samples, seed=1)
(O / "contrasts_swapped.json").write_text(json.dumps(permuted, indent=2) + "\n")
cell = samples.set_index("sample_id").cell_type.to_dict(); subject = samples.set_index("sample_id").subject.to_dict()
rows = [dict(level_a=c["level_a"], level_b=c["level_b"], subject=subject[s]) for c in permuted for s in c["samples_a"] if cell[s] != c["level_a"]]
pd.DataFrame(rows).to_csv(O / "swapped.tsv", sep="\t", index=False)
pd.DataFrame(sorted(pairs), columns=["level_a", "level_b"]).to_csv(O / "level_pairs.tsv", sep="\t", index=False)
print(len(contrasts), "contrasts;", pd.DataFrame(rows).groupby(["level_a", "level_b"]).size().describe().round(1).to_dict(), "swapped subjects per pair;", [len(c["paired_subjects"]) for c in permuted][:5])
