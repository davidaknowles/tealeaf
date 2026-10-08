import json
import sys

import pandas as pd
import pytest

from extra_scripts.summarize_catalog_ec_misfit import main


def shards(root):
    rows = pd.DataFrame([dict(source="parsimony_binary", gene_id=gene, primer=primer, molecules=100., positive_counts=True, integral_counts=True, fit_converged=True, KL=.5, KL_lower_bound=.49, KL_dual_gap=.01, log_probability_bound=-30.) for gene in ("gene1", "gene2") for primer in (0, 1)])
    for index in range(8):
        directory = root / "parsimony_binary" / f"shard_{index}"
        directory.mkdir(parents=True)
        rows.iloc[index::8].to_csv(directory / "genes.tsv.gz", sep="\t", index=False)
        (directory / "settings.json").write_text(json.dumps(dict(shard_count=8, requested_genes=["gene1", "gene2"])))
    return rows


def invoke(monkeypatch, root, output):
    monkeypatch.setattr(sys, "argv", ["collate", "--input-root", str(root), "--output-dir", str(output), "--sources", "parsimony_binary"])
    main()


def test_complete_family_retains_all_gene_primer_rows(tmp_path, monkeypatch):
    shards(tmp_path / "input")
    invoke(monkeypatch, tmp_path / "input", tmp_path / "output")
    result = pd.read_csv(tmp_path / "output/summary.tsv", sep="\t", dtype={"primer": str})
    assert result.loc[result.primer.eq("all"), "gene_primer_rows"].item() == 4
    assert result.loc[result.primer.eq("all"), "incompatible_0.01_rows"].item() == 4


@pytest.mark.parametrize("defect", ["missing", "duplicate", "fractional_bound"])
def test_partial_or_invalid_family_is_not_published(tmp_path, monkeypatch, defect):
    root = tmp_path / "input"
    rows = shards(root)
    target = root / "parsimony_binary/shard_0/genes.tsv.gz"
    if defect == "missing":
        rows.iloc[:0].to_csv(target, sep="\t", index=False)
    elif defect == "duplicate":
        rows.iloc[[0, 0]].to_csv(target, sep="\t", index=False)
    else:
        rows.loc[0, "integral_counts"] = False
        rows.iloc[:1].to_csv(target, sep="\t", index=False)
    with pytest.raises(ValueError):
        invoke(monkeypatch, root, tmp_path / "output")
    assert not (tmp_path / "output/summary.tsv").exists()
