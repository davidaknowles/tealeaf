"""Match the blocks of the consistent local-anchor rule (<run>/blocks_localrule.json.gz) to blocks already collected
and fitted: the current collection (<collection>/all/blocks.json.gz) and the added local-anchor collection
(collection_local). Blocks match when gene, strand, annotated path chains and novel path chains agree as sets (the
precursor window follows from them); a matched block reuses the earlier counts and fits. Writes
<run>/localrule_match.tsv (new key and block id, source, old key and block id) and <run>/blocks_localrule_new.json.gz
(the unmatched blocks, to collect and fit).
Usage: match_blocks.py <run dir> <current collection dir>"""
import gzip, json, sys
from collections import Counter
from pathlib import Path
import pandas as pd

run, current = Path(sys.argv[1]), Path(sys.argv[2])


def content(block):
    n, a = len(block["signatures"]), block["n_annotated"]
    chains = [tuple(tuple(e) for e in path) for path in block["paths"][:n]]
    return (block["gene_id"], block["strand"], frozenset(chains[:a]), frozenset(chains[a:]))


known = {}
for source, path in (("current", current / "all/blocks.json.gz"), ("local", run / "collection_local/all/blocks.json.gz")):
    for key, block in json.load(gzip.open(path, "rt")).items():
        known.setdefault(content(block), (source, key, block["block_id"]))
new = json.load(gzip.open(run / "blocks_localrule.json.gz", "rt"))
rows, unmatched = [], {}
for key, block in new.items():
    hit = known.get(content(block))
    if hit is None:
        unmatched[key] = block
        rows.append((key, block["block_id"], "new", "", ""))
    else:
        rows.append((key, block["block_id"], *hit))
table = pd.DataFrame(rows, columns=["key", "block_id", "source", "old_key", "old_block_id"])
table.to_csv(run / "localrule_match.tsv", sep="\t", index=False)
with gzip.open(run / "blocks_localrule_new.json.gz", "wt") as handle:
    json.dump(unmatched, handle)
print(run, len(new), "blocks:", Counter(table.source).most_common(), "| current blocks not reused:", len({k for s, k, _ in known.values() if s == "current"} - set(table.old_key)))
