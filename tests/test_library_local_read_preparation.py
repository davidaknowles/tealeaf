import argparse
import json

import pandas as pd
import pytest

from extra_scripts.audit_event_local_read_support import file_hash
from extra_scripts.audit_library_local_read_support import prepare


def setup(args):
    args.source_root.mkdir()
    args.origins.mkdir()
    groups = dict(AAAA=["subject", "cell", "random hexamer"], TTTT=["subject", "cell", "poly(dT)"])
    source = dict(events=[], barcode_groups=groups, bams=[dict(path=f"run{index}/reads.bam", size=10, mtime_ns=1) for index in range(2)], diagnostic_inputs={}, scope="old")
    for name in ("selected_cases.tsv.gz", "feature_definitions.tsv"):
        (args.source_root / name).write_bytes(b"frozen")
    source["selected_sha256"] = file_hash(args.source_root / "selected_cases.tsv.gz")
    recipe_path = args.source_root / "recipe.json"
    recipe_path.write_text(json.dumps(source))
    for index, packet in enumerate(source["bams"]):
        shard = args.source_root / f"shard_{index}"
        shard.mkdir()
        (shard / "alignments.bai").write_bytes(b"index")
        (shard / "manifest.json").write_text(json.dumps(dict(complete=True, input=packet, recipe_sha256=file_hash(recipe_path))))
    origins = [dict(barcode=barcode, subject=group[0], cell_type=group[1], primer=group[2], library="lib", status="unique") for barcode, group in groups.items()]
    pd.DataFrame(origins).to_csv(args.origins / "barcode_library.tsv.gz", sep="\t", index=False)
    pd.DataFrame([dict(library="lib", run_a="run0", run_b="run1")]).to_csv(args.origins / "library_runs.tsv", sep="\t", index=False)
    (args.origins / "manifest.json").write_text(json.dumps(dict(complete_unique_library_assignment=True, source_recipe_sha256=file_hash(recipe_path))))


def arguments(tmp_path):
    return argparse.Namespace(source_root=tmp_path / "source", origins=tmp_path / "origins", output_dir=tmp_path / "union")


def test_library_preparation_retains_original_cells_cases_and_scanned_indexes(tmp_path):
    args = arguments(tmp_path)
    setup(args)
    prepare(args)
    recipe = json.loads((args.output_dir / "recipe.json").read_text())
    assert recipe["library_union"]
    assert len(recipe["bams"]) == 1
    assert recipe["bams"][0]["barcode_groups"] == recipe["barcode_groups"]
    assert len(recipe["bams"][0]["runs"]) == 2
    assert (args.output_dir / "selected_cases.tsv.gz").read_bytes() == b"frozen"
    assert len(recipe["diagnostic_inputs"]) == 4


@pytest.mark.parametrize("problem", ["missing_barcode", "changed_annotation", "reused_origin", "duplicate_run", "incomplete_scan"])
def test_library_preparation_rejects_changed_cell_scope_and_incomplete_partition(tmp_path, problem):
    args = arguments(tmp_path)
    setup(args)
    if problem in ("missing_barcode", "changed_annotation", "reused_origin"):
        path = args.origins / "barcode_library.tsv.gz"
        table = pd.read_csv(path, sep="\t")
        if problem == "missing_barcode":
            table = table.iloc[:1]
        elif problem == "changed_annotation":
            table.loc[0, "subject"] = "other"
        else:
            table.loc[0, "status"] = "ambiguous_or_reused"
        table.to_csv(path, sep="\t", index=False)
    elif problem == "duplicate_run":
        path = args.origins / "library_runs.tsv"
        table = pd.read_csv(path, sep="\t")
        table.loc[0, "run_b"] = "run0"
        table.to_csv(path, sep="\t", index=False)
    else:
        path = args.source_root / "shard_1/manifest.json"
        receipt = json.loads(path.read_text())
        receipt["complete"] = False
        path.write_text(json.dumps(receipt))
    with pytest.raises(ValueError):
        prepare(args)
    assert not args.output_dir.exists()
