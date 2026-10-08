import json
import sys

import numpy as np
import pandas as pd
import pytest

from extra_scripts.summarize_sequence_ec_opportunities import main, summarize


def inputs(root):
    for source in ("parsimony_binary", "original_binary"):
        rows = pd.DataFrame([dict(source=source, gene_id=gene, primer=primer, read_length=150, end_window="full", background_fraction=.01, status="ok", full_positive_support=(gene == "gene1"), retained_molecule_fraction=1. if gene == "gene1" else .5, exact_class_fraction=.75, input_molecules=100., sequence_KL=.11, sequence_KL_lower_bound=.1, baseline_same_rows_KL=.91, baseline_same_rows_KL_lower_bound=.9) for gene in ("gene1", "gene2") for primer in (0, 1)])
        for shard in range(8):
            directory = root / source / f"shard_{shard}"
            directory.mkdir(parents=True)
            rows.iloc[shard::8].to_csv(directory / "genes.tsv.gz", sep="\t", index=False)
            (directory / "settings.json").write_text(json.dumps(dict(source=source, requested_genes=["gene1", "gene2"], read_lengths=[150], terminal_start_windows=["full"], shard_count=8, background_fraction=.01)))


def invoke(monkeypatch, root, output):
    monkeypatch.setattr(sys, "argv", ["collate", "--input-root", str(root), "--output-dir", str(output)])
    main()


def test_partial_counts_never_become_full_model_validation(tmp_path, monkeypatch):
    inputs(tmp_path / "input")
    invoke(monkeypatch, tmp_path / "input", tmp_path / "output")
    result = pd.read_csv(tmp_path / "output/summary.tsv", sep="\t")
    assert result.full_positive_support_rows.eq(1).all()
    assert result.requested_gene_primer_rows.eq(2).all()
    assert result.pooled_retained_molecule_fraction.eq(.75).all()
    np.testing.assert_allclose(result.full_positive_support_median_sequence_KL_certificate_gap, .01)


def test_string_false_is_not_full_support_and_invalid_boolean_is_rejected():
    table = pd.DataFrame([dict(gene_id="gene", status="ok", full_positive_support="False", retained_molecule_fraction=.5, input_molecules=100., sequence_KL=.11, baseline_same_rows_KL=.91, sequence_KL_lower_bound=.1, baseline_same_rows_KL_lower_bound=.9)])
    assert summarize(table)["full_positive_support_rows"] == 0
    table["full_positive_support"] = "unknown"
    with pytest.raises(ValueError, match="boolean"):
        summarize(table)


@pytest.mark.parametrize("defect", ["missing", "mixed_background"])
def test_incomplete_or_inconsistent_family_is_not_published(tmp_path, monkeypatch, defect):
    inputs(tmp_path / "input")
    target = tmp_path / "input/original_binary/shard_0/genes.tsv.gz"
    row = pd.read_csv(target, sep="\t")
    if defect == "missing":
        row = row.iloc[:0]
    else:
        row["background_fraction"] = .05
    row.to_csv(target, sep="\t", index=False)
    with pytest.raises(ValueError):
        invoke(monkeypatch, tmp_path / "input", tmp_path / "output")
    assert not (tmp_path / "output/summary.tsv").exists()
