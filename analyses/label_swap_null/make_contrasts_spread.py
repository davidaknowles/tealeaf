"""Within-subject label-swap null on the two halves of GSE320507, for 21 pairs of seven cell types spread over the range
of abundance. One swap per (pair, subject), drawn by permute_paired_contrasts on the full-data contrasts (seed 1) and
applied in both halves, for every method: junction-benchmark fold contrasts (LeafCutter, scQuint), rMATS fold contrasts
(samples and count-pass indices), and the Tealeaf table (level_a, level_b, subject) for run_local_path_tests.py --swapped-pairs."""
import itertools, json, sys
from pathlib import Path
import pandas as pd
sys.path.insert(0, ".")
from tealeaf.sc.junction_benchmark import permute_paired_contrasts
B = Path.home() / "tealeaf_runs/pbmc/GSE320507/junction_benchmark"; O = Path(".cache/label_swap2")
TYPES = ["Tcm/Naive_helper_T_cells", "Classical_monocytes", "Memory_B_cells", "Regulatory_T_cells", "pDC", "CRTAM+_gamma-delta_T_cells", "Plasma_cells"]
pairs = {tuple(sorted(p)) for p in itertools.combinations(TYPES, 2)}
key = lambda c: tuple(sorted((c["level_a"], c["level_b"])))
samples = pd.read_csv(B / "pseudobulk_junctions_samples.tsv", sep="\t", dtype={"subject": str})
full = [c for c in json.loads((B / "contrasts.json").read_text()) if key(c) in pairs]
assert len(full) == 21, len(full)
permuted = permute_paired_contrasts(full, samples, seed=1)
cell = samples.set_index("sample_id").cell_type.to_dict(); subject = samples.set_index("sample_id").subject.to_dict()
swaps = {key(c): {subject[s] for s in c["samples_a"] if cell[s] != c["level_a"]} for c in permuted}
pd.DataFrame([dict(level_a=a, level_b=b, subject=u) for (a, b), us in swaps.items() for u in sorted(us)]).to_csv(O / "swapped.tsv", sep="\t", index=False)
pd.DataFrame(sorted(pairs), columns=["level_a", "level_b"]).to_csv(O / "level_pairs.tsv", sep="\t", index=False)

def swap(c):
    c = dict(c); s = swaps[key(c)]
    for i, sample in enumerate(c["samples_a"]):
        if subject[sample] in s:
            for field in ("samples", "indices"):
                if f"{field}_a" in c:
                    c[f"{field}_a"], c[f"{field}_b"] = list(c[f"{field}_a"]), list(c[f"{field}_b"])
                    c[f"{field}_a"][i], c[f"{field}_b"][i] = c[f"{field}_b"][i], c[f"{field}_a"][i]
    return c

for k in (0, 1):
    junction = [c for c in json.loads((B / f"reproducibility/fold{k}/contrasts.json").read_text()) if key(c) in pairs]
    rmats = [dict(c, original_index=c["index"], index=i) for i, c in enumerate(c for c in json.loads((B / f"rmats/split/fold{k}_contrasts.json").read_text()) if key(c) in pairs)]  # run_rmats_comparison.py selects by "index"
    for name, cs in (("junction", junction), ("rmats", rmats)):
        for c in cs:
            assert [subject[s] for s in c["samples_a"]] == [subject[s] for s in c["samples_b"]]
        (O / f"fold{k}_{name}_contrasts.json").write_text(json.dumps([swap(c) for c in cs], indent=1) + "\n")
        (O / f"fold{k}_{name}_contrasts_real.json").write_text(json.dumps(cs, indent=1) + "\n")
    print(f"fold {k}: {len(junction)} junction contrasts, {len(rmats)} rMATS contrasts, subjects per pair {sorted(len(c['paired_subjects']) for c in junction)}")
print("swapped subjects per pair:", sorted(len(s) for s in swaps.values()))
