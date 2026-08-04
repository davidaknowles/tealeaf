import importlib.util
from pathlib import Path
from types import SimpleNamespace

import pandas as pd


SCRIPT = Path(__file__).parents[1] / "extra_scripts" / "run_junction_comparator.py"
SPEC = importlib.util.spec_from_file_location("run_junction_comparator", SCRIPT)
COMPARATOR = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(COMPARATOR)


def test_normalize_majiq_prefers_conservative_bootstrap_pvalue(tmp_path):
    contrasts = tmp_path / "contrasts.json"
    contrasts.write_text(
        """[
          {
            "contrast_id": "condition__A__control__case",
            "effect": "condition",
            "stratum": "A",
            "level_a": "control",
            "level_b": "case"
          }
        ]\n"""
    )
    raw = tmp_path / "raw.tsv"
    raw.write_text(
        "# MAJIQ metadata\n"
        "lsv_id\tec_idx\tmannwhitneyu-raw_pvalue\t"
        "mannwhitneyu-approximate_pvalue_quantiles_0.950\n"
        "gene1:s:1-2\t0\t0.001\t0.04\n"
        "gene1:s:1-2\t1\t0.002\t0.20\n"
    )
    output = tmp_path / "normalized.tsv"
    COMPARATOR.command_normalize_majiq(
        SimpleNamespace(
            contrasts=contrasts,
            contrast_index=0,
            input=raw,
            output=output,
        )
    )

    result = pd.read_csv(output, sep="\t")
    assert result["feature_id"].is_unique
    assert result["p_value"].tolist() == [0.04, 0.2]
    assert result["majiq_pvalue_column"].unique().tolist() == [
        "mannwhitneyu-approximate_pvalue_quantiles_0.950"
    ]


def test_normalize_leafcutter_omits_untestable_clusters(tmp_path):
    contrasts = tmp_path / "contrasts.json"
    contrasts.write_text(
        """[
          {
            "contrast_id": "condition__A__control__case",
            "effect": "condition",
            "stratum": "A",
            "level_a": "control",
            "level_b": "case"
          }
        ]\n"""
    )
    raw = tmp_path / "raw.tsv"
    raw.write_text("cluster\tp\nclu_1\tNA\nclu_2\tNA\n")
    output = tmp_path / "normalized.tsv"
    COMPARATOR.command_normalize_leafcutter(
        SimpleNamespace(
            contrasts=contrasts,
            contrast_index=0,
            input=raw,
            output=output,
        )
    )

    assert output.read_text() == "\n"


def test_summarize_accepts_method_directories(tmp_path):
    input_dir = tmp_path / "method"
    input_dir.mkdir()
    pd.DataFrame(
        {
            "method": ["test"],
            "contrast_id": ["cell_type__A__B"],
            "effect": ["cell_type"],
            "feature_id": ["f1"],
            "p_value": [0.01],
            "significant": [True],
        }
    ).to_csv(input_dir / "0.tsv", sep="\t", index=False)
    args = SimpleNamespace(
        input=None,
        input_dir=[input_dir],
        output=tmp_path / "all.tsv",
        omnibus_output=tmp_path / "omnibus.tsv",
        summary_output=tmp_path / "summary.tsv",
    )
    COMPARATOR.command_summarize(args)
    summary = pd.read_csv(args.summary_output, sep="\t")
    assert summary.to_dict("records") == [
        {"method": "test", "effect": "cell_type", "tests": 1, "discoveries": 1}
    ]


def test_summarize_skips_blank_no_test_files(tmp_path):
    input_dir = tmp_path / "method"
    input_dir.mkdir()
    (input_dir / "no_tests.tsv").write_text("\n")
    args = SimpleNamespace(
        input=None,
        input_dir=[input_dir],
        output=tmp_path / "all.tsv",
        omnibus_output=tmp_path / "omnibus.tsv",
        summary_output=tmp_path / "summary.tsv",
    )
    COMPARATOR.command_summarize(args)
    assert pd.read_csv(args.summary_output, sep="\t").empty
