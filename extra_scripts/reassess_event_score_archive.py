"""Reassess a complete event shard from saved scores, without selecting events.

Older shards lacking component archives may be refitted explicitly. Run only
after the source shard has completed; source results are never overwritten.
"""

import argparse
import json
from pathlib import Path
import subprocess
import sys
from types import SimpleNamespace

import numpy as np
import pandas as pd

from extra_scripts.run_suppa2_tealeaf_hybrid import mixed_event_record, MIXED_SCORE_EMPTY_COLUMNS
from tealeaf.sc.path_score_mixed import binary_score_components_from_records


def read_table(path):
    try:
        return pd.read_csv(path, sep="\t")
    except pd.errors.EmptyDataError:
        return pd.DataFrame()


def reassess(source, output, *, information_metric, refit_missing=False):
    if source.resolve() == output.resolve() or output.exists():
        raise ValueError("use a new output directory, never overwrite source or previous assessment")
    settings = json.loads((source / "settings.json").read_text())
    summary = json.loads((source / "summary.json").read_text())
    arguments = dict(settings["arguments"])
    if arguments["inference"] != "mixed-score":
        raise ValueError("only mixed-score archives can be reassessed")
    arguments.update(output_dir=str(output), information_metric=information_metric, export_score_components=True)
    contexts_path, scores_path = source / "score_contexts.tsv.gz", source / "subject_scores.tsv.gz"
    if not contexts_path.exists() or not scores_path.exists():
        if not refit_missing:
            raise FileNotFoundError("source predates score capture, explicitly request missing-archive refit")
        command = [sys.executable, str(Path(__file__).with_name("run_suppa2_tealeaf_hybrid.py"))]
        for name, value in arguments.items():
            flag = "--" + name.replace("_", "-")
            if name == "export_score_components":
                command.append(flag)
            elif value is True:
                command.append(flag)
            elif value is not False and value is not None:
                command.extend((flag, str(value)))
        subprocess.run(command, check=True)
        return "refitted missing archive"
    contexts, scores = read_table(contexts_path), read_table(scores_path)
    original = read_table(source / "paired_path.tsv")
    failures = json.loads((source / "failures.json").read_text())
    original_ids = set(original.test_id) if len(original) else set()
    failed_ids = {row["test_id"] for row in failures}
    if len(original_ids) != len(original) or len(failed_ids) != len(failures) or original_ids & failed_ids or len(original) + len(failures) != summary["tests_in_shard"]:
        raise ValueError("source shard is incomplete or contains duplicate failures")
    context_ids = set(contexts.test_id) if len(contexts) else set()
    if len(context_ids) != len(contexts) or not original_ids <= context_ids or not context_ids <= original_ids | failed_ids:
        raise ValueError("archive does not cover the declared source family")
    score_ids = set(scores.test_id) if len(scores) else set()
    if score_ids != context_ids or (len(scores) and scores.duplicated(["test_id", "subject"]).any()):
        raise ValueError("archive subjects do not match fitted contexts")
    retained_failures = [row for row in failures if row["test_id"] not in context_ids]
    observed, null, usage = [], [], []
    groups = dict(tuple(scores.groupby("test_id", sort=False))) if len(scores) else {}
    args = SimpleNamespace(**arguments)
    for context in contexts.to_dict("records"):
        if context["model_version"] != settings["model_version"] or context["report_pseudocount"] != args.report_pseudocount:
            raise ValueError("archived model/reporting recipe differs from source settings")
        key = context["test_id"]
        records = groups[key].to_dict("records")
        if len(records) != context["n_expected_subjects"] or context["n_samples"] != 2 * len(records):
            raise ValueError("subject archive is incomplete for its declared context")
        components = binary_score_components_from_records(records, score_coordinate=context["score_coordinate"])
        event = SimpleNamespace(event_id=context["block_id"], feature_id=context["feature_id"], event_type=context["event_type"])
        mass = context["baseline_event_mass"]
        try:
            row, draws, reports = mixed_event_record(SimpleNamespace(n_isoforms=context["n_isoforms"]), np.array([0, 1, -1]), np.zeros(context["n_samples"]), components.subject_ids, np.array([mass / 2, mass / 2, 1 - mass]), event, context["gene_id"], (context["level_a"], context["level_b"]), np.array([context["median_gene_umis"]]), context["n_ecs"], args, components=components)
            if row["test_id"] != key:
                raise ValueError("context identity changed during reassessment")
            observed.append(row)
            null.extend(draws)
            usage.extend(reports)
        except (ValueError, np.linalg.LinAlgError) as error:
            retained_failures.append(dict(test_id=key, error=repr(error)))
    if len(observed) + len(retained_failures) != summary["tests_in_shard"]:
        raise ValueError("reassessment lost declared hypotheses")
    output.mkdir(parents=True)
    pd.DataFrame(observed, columns=None if observed else MIXED_SCORE_EMPTY_COLUMNS).to_csv(output / "paired_path.tsv", sep="\t", index=False)
    pd.DataFrame(null, columns=None if null else ("test_id", "block_id", "replicate", "p_value")).to_csv(output / "paired_path_null.tsv.gz", sep="\t", index=False)
    if args.export_path_usage:
        pd.DataFrame(usage).to_csv(output / "path_usage.tsv.gz", sep="\t", index=False)
    settings["arguments"] = arguments
    (output / "settings.json").write_text(json.dumps(settings, indent=2) + "\n")
    (output / "summary.json").write_text(json.dumps({**summary, "completed": len(observed), "failures": len(retained_failures)}, indent=2) + "\n")
    (output / "failures.json").write_text(json.dumps(retained_failures, indent=2) + "\n")
    (output / "reassessment.json").write_text(json.dumps({"source": str(source), "method": "saved complete-subject scores, no count or reporting refit", "information_metric": information_metric, "production_changes": False}, indent=2) + "\n")
    return "reused score archive"


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-dir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--information-metric", choices=("absolute", "reference"), required=True)
    parser.add_argument("--refit-missing", action="store_true")
    args = parser.parse_args()
    print(reassess(args.input_dir, args.output_dir, information_metric=args.information_metric, refit_missing=args.refit_missing), flush=True)


if __name__ == "__main__":
    main()
