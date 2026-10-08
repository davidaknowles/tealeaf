#!/usr/bin/env python3
"""Selected-panel finite-count fit pilot, never a full replication benchmark."""

import argparse
import hashlib
import json
from pathlib import Path
import time

import numpy as np
import pandas as pd

from extra_scripts.audit_event_local_read_support import file_hash
from tealeaf.sc.conditional_read_odds import ConditionalReadOdds, conditional_read_odds_test
from tealeaf.sc.local_read_mixed import LocalReadMixed, local_read_mixed_test


PRIMERS = ("poly(dT)", "random hexamer")
MODELS = ("conditional", "unconditional")
VARIANTS = ("all local markers", "junction markers")


def marker_lookup(signatures, *, junction_only):
    """Count keys once per class; unresolved conflicts remain outside both."""
    masks = (4, 8) if junction_only else (5, 10)
    frame = signatures.copy()
    included, excluded = [(frame.signature.to_numpy(dtype=int) & mask) != 0 for mask in masks]
    frame["included"] = frame["count"] * (included & ~excluded)
    frame["excluded"] = frame["count"] * (excluded & ~included)
    columns = ["feature_id", "subject", "cell_type", "primer"]
    return frame.groupby(columns)[["included", "excluded"]].sum().to_dict("index")


def count_tensor(feature, subjects, levels, lookup):
    """Keep every declared subject/primer/type, including all-zero strata."""
    if len(set(subjects)) != len(subjects) or len(levels) != 2 or levels[0] == levels[1]:
        raise ValueError("unique subjects and two different cell-type levels required")
    values = np.zeros((len(subjects), len(PRIMERS), 2, 2), dtype=np.int64)
    for u, subject in enumerate(subjects):
        for p, primer in enumerate(PRIMERS):
            for c, level in enumerate(levels):
                record = lookup.get((feature, subject, level, primer), {})
                values[u, p, c] = [record.get(key, 0) for key in ("included", "excluded")]
    return values


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", type=Path, required=True)
    parser.add_argument("--shard-index", type=int, required=True)
    parser.add_argument("--shard-count", type=int, default=16)
    parser.add_argument("--output-root", type=Path, required=True)
    args = parser.parse_args()
    if not 0 <= args.shard_index < args.shard_count:
        raise ValueError("declared shard required")
    folder = args.output_root / f"shard_{args.shard_index}"
    if folder.exists():
        raise ValueError("new shard required, never overwrite an earlier pilot")
    paths = [args.source / name for name in ("manifest.json", "diagnostics.tsv.gz", "subject_support.tsv.gz", "run_support_signatures.tsv.gz")]
    manifest = json.loads(paths[0].read_text())
    if not len(manifest["shards"]) == 8 or any(not row["complete"] or "library" not in row["input"] for row in manifest["shards"]):
        raise ValueError("all eight library-union sources must be complete")
    cases = pd.read_csv(paths[1], sep="\t").sort_values(["fold", "test_id"]).reset_index(drop=True)
    if cases.duplicated(["fold", "test_id"]).any() or len(cases) != manifest["requested_tests"] or len(cases) != manifest["diagnostics"]:
        raise ValueError("complete unchanged selected family required")
    support = pd.read_csv(paths[2], sep="\t", dtype={"subject": str})
    if support.duplicated(["fold", "test_id", "subject"]).any():
        raise ValueError("unique subject support rows required")
    signatures = pd.read_csv(paths[3], sep="\t", dtype={"subject": str})
    lookups = {variant: marker_lookup(signatures, junction_only=variant == "junction markers") for variant in VARIANTS}
    rows = []
    requested = cases.iloc[args.shard_index::args.shard_count]
    for case in requested.itertuples(index=False):
        subjects = sorted(support.loc[support.fold.eq(case.fold) & support.test_id.eq(case.test_id), "subject"])
        if not subjects:
            raise ValueError("every case must retain its original subject rows")
        for variant, lookup in lookups.items():
            counts = count_tensor(case.feature_id, subjects, case.test_id.split("|")[-2:], lookup)
            digest = hashlib.sha256(counts.tobytes()).hexdigest()
            for model in MODELS:
                start = time.monotonic()
                record = dict(fold=case.fold, test_id=case.test_id, panel=case.panel, feature_id=case.feature_id, model=model, variant=variant, counts_sha256=digest, n_local_included_keys=int(counts[..., 0].sum()), n_local_excluded_keys=int(counts[..., 1].sum()), requested_subjects=len(subjects), p_value=1., converged=False, error="")
                try:
                    result = conditional_read_odds_test(ConditionalReadOdds(counts), nodes=21) if model == "conditional" else local_read_mixed_test(LocalReadMixed(counts), nodes=11)
                    record.update(result)
                except (ValueError, np.linalg.LinAlgError) as exc:
                    record["error"] = str(exc)
                record["runtime_seconds"] = time.monotonic() - start
                rows.append(record)
        print(f"{len(rows)} requested model/marker fits complete", flush=True)
    folder.mkdir(parents=True)
    pd.DataFrame(rows).to_csv(folder / "tests.tsv.gz", sep="\t", index=False)
    receipt = dict(source_hashes={str(path): file_hash(path) for path in paths}, shard_index=args.shard_index, shard_count=args.shard_count, requested_cases=len(requested), whole_selected_cases=len(cases), models=MODELS, variants=VARIANTS, requested_fits=len(requested) * len(MODELS) * len(VARIANTS), completed_fits=len(rows), scope="frozen strong-tail and matched weak-control numerical/local-effect pilot, not full-family FDR, discovery counts, own-ranked LR or unbiased split replication", failure_policy="all requested cases retained, missing/numerical fits remain p=1", production_changes=False)
    (folder / "manifest.json").write_text(json.dumps(receipt, indent=2) + "\n")
    print(json.dumps(receipt, indent=2), flush=True)


if __name__ == "__main__":
    main()
