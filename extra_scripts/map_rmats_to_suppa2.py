#!/usr/bin/env python3
"""Map rMATS event coordinates to SUPPA2 event coordinates."""

from __future__ import annotations

import argparse
from collections import defaultdict
import re
from pathlib import Path

import pandas as pd


EVENT_TYPES = {"SE": "SE", "A3SS": "A3", "A5SS": "A5", "MXE": "MX", "RI": "RI"}
START_FIELDS = {"exonStart_0base", "longExonStart_0base", "shortES", "flankingES", "1stExonStart_0base", "2ndExonStart_0base", "riExonStart_0base", "upstreamES", "downstreamES"}


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--rmats-events", required=True, type=Path)
    parser.add_argument("--suppa-ioe", required=True, type=Path)
    parser.add_argument("--tests", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--mapping-output", required=True, type=Path)
    return parser.parse_args()


def suppa_index(path):
    index = defaultdict(list)
    for line in path.read_text().splitlines()[1:]:
        fields = line.split("\t")
        if len(fields) < 3:
            continue
        gene = fields[1].split(".")[0]
        event_id = fields[2]
        event_type = event_id.split(";", 1)[1].split(":", 1)[0]
        coordinates = frozenset(int(value) for value in re.findall(r"(?<![A-Za-z])\d+", event_id.split(";", 1)[1]))
        index[(gene, event_type)].append((coordinates, event_id))
    return index


def row_coordinates(event_type, row):
    fields_by_type = {"SE": ["upstreamEE", "exonStart_0base", "exonEnd", "downstreamES"], "A3SS": ["longExonStart_0base", "longExonEnd", "shortES", "shortEE", "flankingES", "flankingEE"], "A5SS": ["longExonStart_0base", "longExonEnd", "shortES", "shortEE", "flankingES", "flankingEE"], "MXE": ["upstreamEE", "1stExonStart_0base", "1stExonEnd", "2ndExonStart_0base", "2ndExonEnd", "downstreamES"], "RI": ["upstreamEE", "riExonStart_0base", "riExonEnd", "downstreamES"]}
    return {int(row[field]) + (1 if field in START_FIELDS else 0) for field in fields_by_type[event_type]}


def build_mapping(event_dir, ioe):
    index = suppa_index(ioe)
    mapping = []
    for event_type, native_type in EVENT_TYPES.items():
        table = pd.read_csv(event_dir / f"fromGTF.{event_type}.txt", sep="\t", dtype=str)
        for row in table.to_dict("records"):
            gene = str(row["GeneID"]).strip('"').split(";")[0].split(".")[0]
            candidates = [event_id for coordinates, event_id in index[(gene, native_type)] if coordinates.issubset(row_coordinates(event_type, row))]
            if len(candidates) == 1:
                mapping.append({"event_type": event_type, "event_id": str(row["ID"]), "suppa_event_id": candidates[0]})
    return pd.DataFrame(mapping).drop_duplicates(["event_type", "event_id"])


def main():
    args = parse_args()
    mapping = build_mapping(args.rmats_events, args.suppa_ioe)
    args.mapping_output.parent.mkdir(parents=True, exist_ok=True)
    mapping.to_csv(args.mapping_output, sep="\t", index=False)
    tests = pd.read_csv(args.tests, sep="\t", compression="infer", low_memory=False)
    lookup = {(row.event_type, str(row.event_id)): row.suppa_event_id for row in mapping.itertuples(index=False)}
    tests["suppa_event_id"] = [lookup.get((row.event_type, str(row.event_id))) for row in tests.itertuples(index=False)]
    tests = tests[tests.suppa_event_id.notna()].copy()
    tests["feature_id"] = "SUPPA2:" + tests.suppa_event_id.astype(str)
    tests["event_id"] = tests.suppa_event_id
    tests.to_csv(args.output, sep="\t", index=False, compression="gzip")
    print(f"mapped {len(mapping):,} rMATS events and retained {len(tests):,} tests")


if __name__ == "__main__":
    main()
