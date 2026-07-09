#!/usr/bin/env python3
"""Check whether Salmon weighted EC rows can be joined to alevin-fry EC rows.

The join key is the sorted transcript-id set for an equivalence class. Salmon
weighted dumps can contain multiple rows with the same transcript set when
range factorization is active, so this script reports both missing and
ambiguous joins.
"""

from __future__ import annotations

import argparse
import gzip
import json
from collections import Counter, defaultdict
from pathlib import Path
from typing import Iterable, TextIO


def smart_open(path: Path) -> TextIO:
    if path.suffix == ".gz":
        return gzip.open(path, "rt")
    return path.open()


def read_lines(path: Path) -> list[str]:
    with smart_open(path) as handle:
        return [line.rstrip("\n") for line in handle]


def parse_salmon_eqclasses(path: Path) -> tuple[list[str], dict[tuple[int, ...], list[tuple[float, ...]]]]:
    by_key: dict[tuple[int, ...], list[tuple[float, ...]]] = defaultdict(list)
    with smart_open(path) as handle:
        num_targets = int(next(handle))
        num_eq = int(next(handle))
        names = [next(handle).rstrip("\n") for _ in range(num_targets)]

        for row_num, line in enumerate(handle, start=1):
            fields = line.split()
            if not fields:
                continue
            group_size = int(fields[0])
            tx_ids = tuple(sorted(int(value) for value in fields[1 : 1 + group_size]))
            remainder = fields[1 + group_size :]
            if len(remainder) == group_size + 1:
                weights = tuple(float(value) for value in remainder[:group_size])
            elif len(remainder) == 1:
                weights = ()
            else:
                raise ValueError(
                    f"Could not parse Salmon EC row {row_num}: expected count or weights+count"
                )
            by_key[tx_ids].append(weights)

    observed_eq = sum(len(values) for values in by_key.values())
    if observed_eq != num_eq:
        raise ValueError(f"Salmon header says {num_eq} ECs, parsed {observed_eq}")
    return names, by_key


def used_ec_columns_from_mtx(path: Path) -> set[int]:
    used: set[int] = set()
    with smart_open(path) as handle:
        for line in handle:
            if line.startswith("%"):
                continue
            fields = line.split()
            if len(fields) == 3:
                break
        for line in handle:
            fields = line.split()
            if len(fields) >= 2:
                used.add(int(fields[1]) - 1)
    return used


def iter_alevin_ecclasses(path: Path) -> Iterable[tuple[int, tuple[int, ...]]]:
    with smart_open(path) as handle:
        next(handle)
        next(handle)
        for line in handle:
            fields = line.split()
            if not fields:
                continue
            yield int(fields[-1]), tuple(sorted(int(value) for value in fields[:-1]))


def weights_differ(rows: list[tuple[float, ...]], tol: float) -> bool:
    nonempty = [row for row in rows if row]
    if len(nonempty) < 2:
        return False
    first = nonempty[0]
    return any(
        len(row) != len(first) or any(abs(a - b) > tol for a, b in zip(first, row))
        for row in nonempty[1:]
    )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--alevin-dir", required=True, type=Path)
    parser.add_argument("--salmon-eqclasses", required=True, type=Path)
    parser.add_argument("--count-matrix", type=Path, default=None)
    parser.add_argument(
        "--restrict-to-count-matrix",
        action="store_true",
        help="only check alevin EC columns observed in geqc_counts.mtx; can be slow for large matrices",
    )
    parser.add_argument("--output-json", type=Path, default=None)
    parser.add_argument("--max-alevin-rows", type=int, default=None)
    parser.add_argument("--weight-tol", type=float, default=1e-8)
    args = parser.parse_args()

    alevin_features = read_lines(args.alevin_dir / "quants_mat_cols.txt")
    salmon_names, salmon_by_key = parse_salmon_eqclasses(args.salmon_eqclasses)

    used_ecs = None
    count_matrix = args.count_matrix or args.alevin_dir / "geqc_counts.mtx"
    if args.restrict_to_count_matrix and count_matrix.is_file():
        used_ecs = used_ec_columns_from_mtx(count_matrix)

    salmon_key_multiplicity = Counter({key: len(rows) for key, rows in salmon_by_key.items()})
    duplicate_keys = [key for key, count in salmon_key_multiplicity.items() if count > 1]
    duplicate_weight_conflicts = [
        key for key in duplicate_keys if weights_differ(salmon_by_key[key], args.weight_tol)
    ]

    stats = {
        "name_order_match": salmon_names == alevin_features,
        "alevin_num_features": len(alevin_features),
        "salmon_num_features": len(salmon_names),
        "salmon_rows": sum(len(rows) for rows in salmon_by_key.values()),
        "salmon_unique_transcript_sets": len(salmon_by_key),
        "salmon_duplicate_transcript_sets": len(duplicate_keys),
        "salmon_duplicate_weight_conflicts": len(duplicate_weight_conflicts),
        "restricted_to_observed_alevin_ecs": used_ecs is not None,
        "observed_alevin_ecs": len(used_ecs) if used_ecs is not None else None,
        "alevin_rows_checked": 0,
        "matched_unique": 0,
        "matched_ambiguous": 0,
        "missing": 0,
        "missing_examples": [],
        "ambiguous_examples": [],
    }

    for ec_id, key in iter_alevin_ecclasses(args.alevin_dir / "gene_eqclass.txt.gz"):
        if used_ecs is not None and ec_id not in used_ecs:
            continue
        stats["alevin_rows_checked"] += 1
        matches = salmon_by_key.get(key)
        if matches is None:
            stats["missing"] += 1
            if len(stats["missing_examples"]) < 10:
                stats["missing_examples"].append({"ec_id": ec_id, "tx_ids": list(key)})
        elif len(matches) == 1:
            stats["matched_unique"] += 1
        else:
            stats["matched_ambiguous"] += 1
            if len(stats["ambiguous_examples"]) < 10:
                stats["ambiguous_examples"].append(
                    {"ec_id": ec_id, "tx_ids": list(key), "salmon_rows": len(matches)}
                )

        if args.max_alevin_rows is not None and stats["alevin_rows_checked"] >= args.max_alevin_rows:
            break

    text = json.dumps(stats, indent=2)
    print(text)
    if args.output_json is not None:
        args.output_json.write_text(text + "\n")


if __name__ == "__main__":
    main()
