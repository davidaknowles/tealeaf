import argparse
import json

import pandas as pd
import pytest

from extra_scripts.audit_event_local_read_support import collate, file_hash


def fixture(args):
    root, diagnostic, bound = args.output_dir, args.output_dir / "diagnostic", args.bound_root
    root.mkdir()
    feature = "SUPPA2:g;SE:chr1:20-30:40-50:+"
    key = feature + "|cell_type|A|B"
    selected, hashes = [], {}
    for fold in (0, 1):
        folder = diagnostic / f"fold{fold}"
        folder.mkdir(parents=True)
        record = dict(fold=fold, test_id=key, feature_id=feature, panel="strong", n_informative_subjects=4, dominant_subject="s0", maximum_precision_share=.95, fitted_mean=-1.)
        pd.DataFrame([record]).to_csv(folder / "diagnostics.tsv.gz", sep="\t", index=False)
        subjects = pd.DataFrame(dict(test_id=[key] * 4, panel=["strong"] * 4, subject=[f"s{x}" for x in range(4)], retained=[True] * 4, precision_share=[.95, .02, .02, .01], pseudo_effect=[-1.] * 4))
        subjects.to_csv(folder / "subject_influence.tsv.gz", sep="\t", index=False)
        bounds = bound / f"fold{fold}"
        bounds.mkdir(parents=True)
        pd.DataFrame([dict(test_id=key, largest_uniform_lower_bound=.94)]).to_csv(bounds / "diagnostics.tsv.gz", sep="\t", index=False)
        for path in (folder / "diagnostics.tsv.gz", folder / "subject_influence.tsv.gz", bounds / "diagnostics.tsv.gz"):
            hashes[str(path)] = file_hash(path)
        selected.append(record)
    selected_path = root / "selected_cases.tsv.gz"
    pd.DataFrame(selected).to_csv(selected_path, sep="\t", index=False)
    pd.DataFrame([dict(feature_id=feature, event_type="SE", status="ok")]).to_csv(root / "feature_definitions.tsv", sep="\t", index=False)
    packets = [dict(path=f"run{x}/reads.bam", size=100, mtime_ns=5) for x in range(2)]
    recipe = dict(bams=packets, events=[dict(event_id=feature)], diagnostic_root=str(diagnostic), diagnostic_inputs=hashes, selected_sha256=file_hash(selected_path), selection="frozen", scope="diagnostic")
    recipe_path = root / "recipe.json"
    recipe_path.write_text(json.dumps(recipe))
    for index, packet in enumerate(packets):
        shard = root / f"shard_{index}"
        shard.mkdir()
        rows = [dict(feature_id=feature, subject=f"s{x}", cell_type=level, primer="random hexamer", signature=17 if level == "A" else 8, count=40) for x in range(4) for level in ("A", "B")]
        pd.DataFrame(rows).to_csv(shard / "support.tsv.gz", sep="\t", index=False)
        receipt = dict(complete=True, shard=index, input=packet, recipe_sha256=file_hash(recipe_path), requested_events=1)
        (shard / "manifest.json").write_text(json.dumps(receipt))


def test_whole_family_collation_retains_both_folds_and_zero_primer_groups(tmp_path):
    args = argparse.Namespace(output_dir=tmp_path / "source", public_dir=tmp_path / "public", bound_root=tmp_path / "bounds")
    fixture(args)
    collate(args)
    cases = pd.read_csv(args.public_dir / "diagnostics.tsv.gz", sep="\t")
    assert len(cases) == 2
    assert cases.fold.tolist() == [0, 1]
    assert cases.n_subjects_with_five_discriminating_keys_each_type.tolist() == [4, 4]
    assert cases.mean_local_inclusion_difference.tolist() == [-1., -1.]
    assert cases.raw_score_direction_agrees_with_local_mean.all()
    subjects = pd.read_csv(args.public_dir / "subject_support.tsv.gz", sep="\t")
    assert len(subjects) == 8
    assert subjects.included_only_0_DT.eq(0).all()
    assert subjects.included_only_0_RH.eq(80).all()
    assert subjects.excluded_only_1_RH.eq(80).all()
    manifest = json.loads((args.public_dir / "manifest.json").read_text())
    assert manifest["requested_tests"] == manifest["diagnostics"] == 2
    assert manifest["production_changes"] is False


@pytest.mark.parametrize("problem", ["incomplete", "wrong_recipe", "changed_diagnostic", "negative_counts"])
def test_collation_rejects_incomplete_or_changed_sources(tmp_path, problem):
    args = argparse.Namespace(output_dir=tmp_path / "source", public_dir=tmp_path / "public", bound_root=tmp_path / "bounds")
    fixture(args)
    if problem == "changed_diagnostic":
        path = args.output_dir / "diagnostic/fold0/subject_influence.tsv.gz"
        path.write_bytes(b"changed")
    elif problem == "negative_counts":
        path = args.output_dir / "shard_1/support.tsv.gz"
        table = pd.read_csv(path, sep="\t")
        table.loc[0, "count"] = -1
        table.to_csv(path, sep="\t", index=False)
    else:
        path = args.output_dir / "shard_1/manifest.json"
        receipt = json.loads(path.read_text())
        receipt["complete" if problem == "incomplete" else "recipe_sha256"] = False if problem == "incomplete" else "changed"
        path.write_text(json.dumps(receipt))
    with pytest.raises(ValueError):
        collate(args)
    assert not args.public_dir.exists()
