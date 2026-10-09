import json
import sys

import pandas as pd
import pytest

from extra_scripts import prepare_full_library_read_catalog as driver


def run_fixture(tmp_path, monkeypatch, *, bad_qc=False):
    library, qc, output = [tmp_path / name for name in ("library", "qc", "full")]
    library.mkdir()
    qc.mkdir()
    groups = {f"barcode{index}": ["s", "cell", "poly(dT)"] for index in range(8)}
    packets = [dict(library=f"lib{index}", runs=[], barcode_groups={f"barcode{index}": groups[f"barcode{index}"]}) for index in range(8)]
    (library / "recipe.json").write_text(json.dumps(dict(library_union=True, bams=packets, barcode_groups=groups, diagnostic_inputs={"original recipe": "source"})))
    (qc / "manifest.json").write_text(json.dumps(dict(exact_cached_group_and_primer_total_match=not bad_qc, source_read_recipe_sha256="source", retained_production_barcodes=8)))
    pd.DataFrame([dict(barcode=barcode, subject=value[0], cell_type=value[1], primer=value[2]) for barcode, value in groups.items()]).to_csv(qc / "retained_barcodes.tsv.gz", sep="\t", index=False)
    catalog = tmp_path / "events.tsv"
    pd.DataFrame([dict(feature_id=key, gene_id="g", event_type="SE", event_id=f"g;SE:chr1:20-30:40-50:+", included="missing" if key == "unavailable" else "inc", excluded="exc") for key in ("old_supported", "new_without_global_support", "unavailable")]).to_csv(catalog, sep="\t", index=False)
    reference = tmp_path / "reference.tsv"
    pd.DataFrame([dict(source="original_binary", gene_id="g", feature_id="old_supported", n_contrasts=6)]).to_csv(reference, sep="\t", index=False)
    reference_manifest = tmp_path / "reference.json"
    reference_manifest.write_text(json.dumps(dict(candidate_settings=dict(subject_fold=None, min_gene_umis=25, min_celltype_mice=4))))
    gtf = tmp_path / "annotation.gtf"
    gtf.write_text("annotation fixture")
    annotation = dict(g=dict(chromosome="chr1", strand="+", transcripts=dict(inc=((10, 20), (30, 40), (50, 60)), exc=((10, 20), (50, 60)))))
    monkeypatch.setattr(driver, "read_gtf_exons", lambda _: annotation)
    monkeypatch.setattr(sys, "argv", ["prepare", "--library-root", str(library), "--cell-qc", str(qc), "--reference-family", str(reference), "--reference-family-manifest", str(reference_manifest), "--event-catalog", str(catalog), "--gtf", str(gtf), "--output-dir", str(output)])
    driver.main()
    return output


def test_full_local_catalog_keeps_events_without_global_support_and_unavailable_geometry(tmp_path, monkeypatch):
    output = run_fixture(tmp_path, monkeypatch)
    recipe = json.loads((output / "recipe.json").read_text())
    assert set(recipe["declared_events"]) == {"old_supported", "new_without_global_support", "unavailable"}
    assert set(row["event_id"] for row in recipe["events"]) == {"old_supported", "new_without_global_support"}
    family = pd.read_csv(output / "declared_family.tsv.gz", sep="\t")
    assert family.prior_global_isoform_supported.sum() == 1
    definitions = pd.read_csv(output / "feature_definitions.tsv", sep="\t")
    assert len(definitions) == 3
    assert definitions.status.eq("unavailable").sum() == 1
    assert recipe["production_cell_qc"]["exact_cached_group_and_primer_total_match"]


def test_full_catalog_requires_verified_production_qc(tmp_path, monkeypatch):
    with pytest.raises(ValueError, match="production QC"):
        run_fixture(tmp_path, monkeypatch, bad_qc=True)
    assert not (tmp_path / "full").exists()
