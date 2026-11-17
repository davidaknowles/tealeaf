"""Combine the existing Tealeaf models with their local-anchor block models (suffix _L) into new model names, so the
evaluation scripts can read them unchanged: merged fit tables per cohort (BH recomputed within each combined model),
shard_0/summary.json of the existing model, targeted coverage-null tables, long-read mappings and block files.
Brain: internal_native_v4sh + _L -> internal_native_v4shL (tj likewise); blood: pbmc_<GSE>_<structure>_sh + _L ->
pbmc_<GSE>_<structure>_shL."""
import gzip, json, shutil, sys
from pathlib import Path
import pandas as pd
sys.path.insert(0, ".")
from tealeaf.sc.junction_benchmark import benjamini_hochberg

CACHE, NULL = Path(".cache/binary_grid"), Path(".cache/targeted_null")
pairs = [("internal_native_v4sh", "internal_native_v4sh_L", "internal_native_v4shL"), ("tj_native_v4sh", "tj_native_v4sh_L", "tj_native_v4shL")]
for g in ("GSE320507", "GSE327907", "GSE334840"):
    for s in ("internal", "terminal"):
        pairs.append((f"pbmc_{g}_{s}_sh", f"pbmc_{g}_{s}_sh_L", f"pbmc_{g}_{s}_shL"))
for old, new, out in pairs:
    for c in ("0", "1", "full"):
        a, b = CACHE / f"{old}_{c}", CACHE / f"{new}_{c}"
        if not (b / "merged/paired_path.tsv").exists():
            print("missing", b); continue
        t = pd.concat([pd.read_csv(a / "merged/paired_path.tsv", sep="\t", low_memory=False), pd.read_csv(b / "merged/paired_path.tsv", sep="\t", low_memory=False)], ignore_index=True)
        assert not t.test_id.duplicated().any(), (old, c)
        t["fdr"] = benjamini_hochberg(t.p_value.to_numpy())
        (CACHE / f"{out}_{c}/merged").mkdir(parents=True, exist_ok=True)
        (CACHE / f"{out}_{c}/shard_0").mkdir(exist_ok=True)
        t.to_csv(CACHE / f"{out}_{c}/merged/paired_path.tsv", sep="\t", index=False, na_rep="NA")
        shutil.copy(a / "shard_0/summary.json", CACHE / f"{out}_{c}/shard_0/summary.json")
        nulls = [NULL / f"{m}_{c}/null.tsv.gz" for m in (old, new)]
        if all(p.exists() for p in nulls):
            (NULL / f"{out}_{c}").mkdir(parents=True, exist_ok=True)
            pd.concat([pd.read_csv(p, sep="\t") for p in nulls], ignore_index=True).to_csv(NULL / f"{out}_{c}/null.tsv.gz", sep="\t", index=False, na_rep="NA")
        print(out, c, len(t), "tests", int((t.fdr <= .05).sum()), "BH", "with null" if all(p.exists() for p in nulls) else "", flush=True)
# long-read mappings
for old, new, out in ((".cache/internal_priming/compare/v4sh", ".cache/local_anchor/lr_brain_L", ".cache/local_anchor/lr_brain_v4shL"), ("analyses/pbmc/GSE320507_lr_sh", "analyses/pbmc/GSE320507_lr_sh_L", "analyses/pbmc/GSE320507_lr_shL")):
    Path(out).mkdir(parents=True, exist_ok=True)
    m = pd.concat([pd.read_csv(f"{d}/lr_mapping.tsv.gz", sep="\t", low_memory=False) for d in (old, new)], ignore_index=True)
    m.to_csv(f"{out}/lr_mapping.tsv.gz", sep="\t", index=False, na_rep="NA")
    print(out, len(m), "mapped tests")
# block files: existing collection blocks plus local blocks (keys are distinct)
for old, new, out in ((".cache/native_v4/collection_parity_p2_v4/all/blocks.json.gz", ".cache/native_v4/collection_local/all/blocks.json.gz", ".cache/native_v4/blocks_combined_all.json.gz"),
                      (".cache/native_v4/collection_parity_p2_v4/junction/blocks.json.gz", ".cache/native_v4/collection_local/junction/blocks.json.gz", ".cache/native_v4/blocks_combined_junction.json.gz")) + tuple(
                     (f".cache/pbmc/{g}/collection/{k}/blocks.json.gz", f".cache/pbmc/{g}/collection_local/{k}/blocks.json.gz", f".cache/pbmc/{g}/blocks_combined_{k}.json.gz") for g in ("GSE320507", "GSE327907", "GSE334840") for k in ("all", "junction")):
    a, b = json.load(gzip.open(old, "rt")), json.load(gzip.open(new, "rt"))
    assert not set(a) & set(b)
    with gzip.open(out, "wt") as handle:
        json.dump({**a, **b}, handle)
    print(out, len(a), "+", len(b))
