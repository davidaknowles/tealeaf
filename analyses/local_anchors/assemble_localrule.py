"""Assemble the results of the consistent local-anchor rule for one dataset from three sources: the current blocks
("current"), the added local-anchor blocks ("local", models and collection suffix _L / collection_local) and the new
blocks ("new", _N / collection_lrnew), matched by content in <run>/localrule_match.tsv (match_blocks.py). Each source
contributes only its matched blocks, renamed to the local-rule block ids and keys:
- fit tables per cohort and structure -> models <out>_<cohort> (BH recomputed within each model), with targeted
  coverage-null tables where all sources have them;
- long-read mappings -> <lr out>/lr_mapping.tsv.gz;
- collections (all/ and junction/: blocks, counts, opportunities) -> <run>/collection_localrule (and, for GSE320507,
  the donor-group collection disease_localrule and models pbmc_GSE320507_disease_<structure>_lr);
- label-swap null fits of GSE320507 (.cache/label_swap_all/seed*/fold*/tealeaf_<structure>{,_L,_N} and their
  coverage nulls) -> tealeaf_<structure>_LR and tnull_tealeaf_<structure>_LR.
Usage: assemble_localrule.py <brain|GSE320507|GSE327907|GSE334840> [hier]
With hier, the blocks of the hierarchical rule (<run>/hier_match.tsv; gene-wide blocks plus differing local-rule blocks)
are assembled instead, into models suffixed hl (brain *_native_v4hl, blood *_hl), collection_hier, disease_hier,
long-read mappings *_lr_hl / lr_brain_v4hl and swapped fits tealeaf_<structure>_HL."""
import gzip, json, shutil, sys
from pathlib import Path
import pandas as pd
sys.path.insert(0, ".")
from tealeaf.sc.junction_benchmark import benjamini_hochberg

DS = sys.argv[1]
HIER = len(sys.argv) > 2 and sys.argv[2] == "hier"
OUT_TAG, MATCH, COLL_OUT, SWAP_OUT = ("hl", "hier_match.tsv", "collection_hier", "HL") if HIER else ("lr", "localrule_match.tsv", "collection_localrule", "LR")
CACHE, NULL = Path(".cache/binary_grid"), Path(".cache/targeted_null")
if DS == "brain":
    RUN = Path(".cache/native_v4")
    COLL = {"current": RUN / "collection_parity_p2_v4", "local": RUN / "collection_local", "new": RUN / "collection_lrnew"}
    MODELS = {s: {"current": f"{p}_native_v4sh", "local": f"{p}_native_v4sh_L", "new": f"{p}_native_v4sh_N", "out": f"{p}_native_v4{OUT_TAG}"} for s, p in (("internal", "internal"), ("terminal", "tj"))}
    LR = {"current": ".cache/internal_priming/compare/v4sh", "local": ".cache/local_anchor/lr_brain_L", "new": ".cache/local_anchor/lr_brain_N", "out": f".cache/local_anchor/lr_brain_v4{OUT_TAG}"}
else:
    RUN = Path(f".cache/pbmc/{DS}")
    COLL = {"current": RUN / "collection", "local": RUN / "collection_local", "new": RUN / "collection_lrnew"}
    MODELS = {s: {"current": f"pbmc_{DS}_{s}_sh", "local": f"pbmc_{DS}_{s}_sh_L", "new": f"pbmc_{DS}_{s}_sh_N", "out": f"pbmc_{DS}_{s}_{OUT_TAG}"} for s in ("internal", "terminal")}
    LR = {"current": f"analyses/pbmc/{DS}_lr_sh", "local": f"analyses/pbmc/{DS}_lr_sh_L", "new": f"analyses/pbmc/{DS}_lr_sh_N", "out": f"analyses/pbmc/{DS}_lr_{OUT_TAG}"}
match = pd.read_csv(RUN / MATCH, sep="\t", dtype=str, keep_default_na=False)
fill = match.source.eq("new") & match.old_key.eq("")  # local rule: new blocks were fitted under their own ids
match.loc[fill, ["old_key", "old_block_id"]] = match.loc[fill, ["key", "block_id"]].to_numpy()
rename_id = {src: dict(zip(g.old_block_id, g.block_id)) for src, g in match.groupby("source")}
rename_key = {src: dict(zip(g.old_key, g.key)) for src, g in match.groupby("source")}


def renamed(table, src):
    """Rows of table whose block belongs to src's matched blocks, with block_id and test_id renamed."""
    ids = rename_id[src]
    block = table.test_id.str.split("|", n=1).str[0] if "block_id" not in table else table.block_id
    table = table.loc[block.isin(ids)].copy()
    old = table.test_id.str.split("|", n=1)
    table["test_id"] = old.str[0].map(ids) + "|" + old.str[1]
    if "block_id" in table:
        table["block_id"] = table.block_id.map(ids)
    return table


def fits(models, out, collection):
    for c in ("0", "1", "full"):
        parts = [renamed(pd.read_csv(CACHE / f"{models[src]}_{c}/merged/paired_path.tsv", sep="\t", low_memory=False), src) for src in ("current", "local", "new")]
        t = pd.concat(parts, ignore_index=True)
        assert not t.test_id.duplicated().any()
        t["fdr"] = benjamini_hochberg(t.p_value.to_numpy())
        d = CACHE / f"{out}_{c}"; (d / "merged").mkdir(parents=True, exist_ok=True); (d / "shard_0").mkdir(exist_ok=True)
        t.to_csv(d / "merged/paired_path.tsv", sep="\t", index=False, na_rep="NA")
        settings = json.loads((CACHE / f"{models['current']}_{c}/shard_0/summary.json").read_text())
        settings["candidate_settings"]["local_reads"] = str(collection)
        (d / "shard_0/summary.json").write_text(json.dumps(settings, indent=2) + "\n")
        nulls = [NULL / f"{models[src]}_{c}/null.tsv.gz" for src in ("current", "local", "new")]
        note = ""
        if all(p.exists() for p in nulls):
            (NULL / f"{out}_{c}").mkdir(parents=True, exist_ok=True)
            pd.concat([renamed(pd.read_csv(p, sep="\t"), src) for p, src in zip(nulls, ("current", "local", "new"))], ignore_index=True).to_csv(NULL / f"{out}_{c}/null.tsv.gz", sep="\t", index=False, na_rep="NA")
            note = "with coverage null"
        assert all(len(p) for p, src in zip(parts, ("current", "local", "new")) if (match.source == src).any()), f"{out} {c}: a source with matched blocks contributed no tests"
        print(out, c, len(t), "tests", int((t.fdr <= .05).sum()), "BH", " ".join(f"{src} {len(p)}" for src, p in zip(("current", "local", "new"), parts)), note, flush=True)


def collection(sources, out):
    for k in ("all", "junction"):
        (out / k).mkdir(parents=True, exist_ok=True)
        blocks, opportunities = {}, {}
        for src, root in sources.items():
            b = json.load(gzip.open(root / k / "blocks.json.gz", "rt"))
            o = json.load(gzip.open(root / k / "opportunities.json.gz", "rt"))
            for old, new in rename_key[src].items():
                if old in b:
                    blocks[new] = dict(b[old], block_id=rename_id[src][b[old]["block_id"]])
                    if old in o:
                        opportunities[new] = o[old]
        with gzip.open(out / k / "blocks.json.gz", "wt") as h:
            json.dump(blocks, h)
        with gzip.open(out / k / "opportunities.json.gz", "wt") as h:
            json.dump(opportunities, h)
        shutil.copy(sources["current"] / k / "filters.json", out / k / "filters.json")
        first = True
        with gzip.open(out / k / "counts.tsv.gz", "wt") as h:
            for src, root in sources.items():
                keys = rename_key[src]
                for chunk in pd.read_csv(root / k / "counts.tsv.gz", sep="\t", dtype={"subject": str}, chunksize=5_000_000):
                    chunk = chunk.loc[chunk.path_key.isin(keys)].copy()
                    chunk["path_key"] = chunk.path_key.map(keys)
                    chunk.to_csv(h, sep="\t", index=False, header=first); first = False
        print(out / k, len(blocks), "blocks", flush=True)


out_coll = RUN / COLL_OUT
collection(COLL, out_coll)
for s, models in MODELS.items():
    fits(models, models["out"], out_coll / ("all" if s == "internal" else "junction"))
Path(LR["out"]).mkdir(parents=True, exist_ok=True)
pd.concat([renamed(pd.read_csv(f"{LR[src]}/lr_mapping.tsv.gz", sep="\t", low_memory=False), src) for src in ("current", "local", "new")], ignore_index=True).to_csv(f"{LR['out']}/lr_mapping.tsv.gz", sep="\t", index=False, na_rep="NA")
print(LR["out"], "long-read mapping written", flush=True)
if DS == "GSE320507":
    disease_out = RUN / ("disease_hier" if HIER else "disease_localrule")
    collection({"current": RUN / "disease", "local": RUN / "disease_local", "new": RUN / "disease_lrnew"}, disease_out)
    shutil.copy(RUN / "disease/level_pairs.tsv", disease_out / "level_pairs.tsv")
    for s in ("internal", "terminal"):
        m = {"current": f"pbmc_GSE320507_disease_{s}_sh", "local": f"pbmc_GSE320507_disease_{s}_sh_L", "new": f"pbmc_GSE320507_disease_{s}_sh_N"}
        fits(m, f"pbmc_GSE320507_disease_{s}_{OUT_TAG}", disease_out / ("all" if s == "internal" else "junction"))
    for draw in sorted(Path(".cache/label_swap_all").glob("seed*")):
        for k in (0, 1):
            for s in ("internal", "terminal"):
                base = draw / f"fold{k}"
                parts = [renamed(pd.read_csv(base / f"tealeaf_{s}{x}/merged/paired_path.tsv", sep="\t", low_memory=False), src) for x, src in (("", "current"), ("_L", "local"), ("_N", "new"))]
                (base / f"tealeaf_{s}_{SWAP_OUT}/merged").mkdir(parents=True, exist_ok=True)
                pd.concat(parts, ignore_index=True).to_csv(base / f"tealeaf_{s}_{SWAP_OUT}/merged/paired_path.tsv", sep="\t", index=False, na_rep="NA")
                nulls = [renamed(pd.read_csv(base / f"tnull_tealeaf_{s}{x}/null.tsv.gz", sep="\t"), src) for x, src in (("", "current"), ("_L", "local"), ("_N", "new"))]
                (base / f"tnull_tealeaf_{s}_{SWAP_OUT}").mkdir(parents=True, exist_ok=True)
                pd.concat(nulls, ignore_index=True).to_csv(base / f"tnull_tealeaf_{s}_{SWAP_OUT}/null.tsv.gz", sep="\t", index=False, na_rep="NA")
                print(draw.name, k, s, "swap fits", sum(map(len, parts)), flush=True)
