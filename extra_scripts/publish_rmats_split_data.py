#!/usr/bin/env python3
"""Publish completed paired split metrics without substituting full-data tests."""

import argparse
import hashlib
from pathlib import Path
import re

import pandas as pd


def direction_cell(row):
    if row.n_events == 0:
        return r"--- (0 events)"
    if row.n_direction_evaluable == 0:
        return r"--- (unevaluable)"
    return rf"\shortstack{{{100 * row.agreement:.1f}\%\\{int(row.agree):,}/{int(row.n_direction_evaluable):,}}}"


def render_table1_row(metrics, directions):
    tea = metrics.loc["Tealeaf"]
    comparator = metrics.loc["rMATS (paired)"]
    for field in ("shared_gene_pairs", "shared_genes"):
        if tea[field] != comparator[field]:
            raise ValueError("Tealeaf and rMATS shared universes differ")
    cells = ["rMATS (paired JCEC)", f"{int(tea.shared_gene_pairs):,}", f"{int(tea.shared_genes):,}"]
    for method, row in (("Tealeaf", tea), ("rMATS (paired)", comparator)):
        cells.extend([str(int(row.replicated_bh)), f"{row.heldout_nominal_replication:.3f}" if pd.notna(row.heldout_nominal_replication) else "---", f"{row.spearman_logp:.3f}" if pd.notna(row.spearman_logp) else "---", direction_cell(directions.loc[method])])
    return " & ".join(cells) + r"\\"


def replace_once(text, old, new):
    if text.count(old) != 1:
        raise ValueError(f"expected one manuscript target, {old[:100]}")
    return text.replace(old, new)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo-root", type=Path, required=True)
    parser.add_argument("--expected-document-sha", required=True)
    args = parser.parse_args()
    audit = args.repo_root / "analyses/split_coverage_direction"
    output = args.repo_root / "analyses/comparator_suppa_rmats/rmats"
    metrics = pd.read_csv(audit / "table1_gene_metrics.tsv", sep="\t")
    metrics = metrics.loc[metrics.comparison.eq("rMATS (paired)")].set_index("method", verify_integrity=True)
    directions = pd.read_csv(audit / "direction_summary.tsv", sep="\t")
    directions = directions.loc[directions.comparison.eq("rMATS (paired)") & directions.selection.eq("event_BH_union") & directions.coverage_bin.eq("all")].set_index("method", verify_integrity=True)
    row = render_table1_row(metrics, directions)
    (output / "table1_row.tex").write_text(row + "\n")
    metrics.to_csv(output / "split_data_tealeaf_comparison_metrics.tsv", sep="\t", na_rep="NA")
    directions.to_csv(output / "split_data_direction_agreement.tsv", sep="\t", na_rep="NA")
    document = args.repo_root / "docs/differential.tex"
    raw = document.read_bytes()
    if hashlib.sha256(raw).hexdigest() != args.expected_document_sha:
        raise ValueError("manuscript changed while fitting, generated row saved but document not overwritten")
    text = raw.decode()
    targets = re.findall(r"^rMATS \(native, historical\).*?$", text, re.MULTILINE)
    if len(targets) != 1:
        raise ValueError("expected one historical rMATS Table 1 row")
    text = replace_once(text, targets[0], row)
    text = replace_once(text, "Asterisks mark unavailable historical rMATS split event tables, not zero agreement; full-data effects are not substituted.", "The rMATS row uses paired PAIRADISE statistical refits on corrected single-end JCEC counts in both subject halves, not full-data effects.")
    text = replace_once(text, "The rMATS row is explicitly the historical unpaired split audit. The corrected paired-JCEC rMATS rerun is integrated below as a full-data event comparator, rather than being presented as a paired split result that was not rerun.", "The paired-JCEC rMATS row uses independent statistical refits on the same two subject-half manifests and the same shared-gene conjunction procedure. The separate full-data paired-JCEC results remain the source for the long-read event audit.")
    text = replace_once(text, "The historical rMATS split event tables are unavailable, so its aggregate row is retained without fabricated directional values or substituted full-data effects.", "The paired rMATS split rerun replaces the historical unpaired row. Its native mean-PSI difference is converted from level-a-minus-level-b to level-b-minus-level-a before the common directional audit; rounded zero effects are unevaluable, not disagreements.")
    text = replace_once(text, "The split reproducibility row for rMATS above is retained as a historical unpaired audit; a new split paired rerun would be required to replace it with a directly matched conjunction analysis.", "The split reproducibility row for rMATS uses independent paired-JCEC refits and a directly matched conjunction analysis, rather than reusing full-data p-values.")
    text = replace_once(text, "The new paired-JCEC rMATS run restores 1,497 full-data event FDR calls and is included in the event-level table and long-read comparison above, but those calls are not substituted for the historical split row because the paired split manifests were not rerun.", "The full-data paired-JCEC rMATS run restores 1,497 event FDR calls and is included in the long-read comparison above. The separately refitted paired subject-half comparison replaces the historical row in the main reproducibility table.")
    coverage = pd.read_csv(audit / "coverage_correlations.tsv", sep="\t")
    coverage = coverage.loc[coverage.comparison.eq("rMATS (paired)") & coverage.unit.eq("gene")].set_index(["method", "fold"], verify_integrity=True)
    cells = []
    for method in ("Tealeaf", "rMATS (paired)"):
        values = [coverage.loc[(method, fold), "rho_p_coverage"] for fold in (0, 1)]
        cells.append(r"\(" + " / ".join(f"{value:.3f}" if pd.notna(value) else r"\text{NA}" for value in values) + r"\)")
    coverage_row = "rMATS paired JCEC & " + " & ".join(cells) + r"\\"
    anchor = r"scQuint & \(-.425 / -.338\) & \(-.069 / -.080\)\\"
    text = replace_once(text, anchor, anchor + "\n" + coverage_row)
    audit_readme = audit / "README.md"
    if audit_readme.exists():
        description = audit_readme.read_text()
        pattern = r"The historical rMATS row is retained separately,.*?paired full-data rMATS results are not a substitute for split effects\."
        description, count = re.subn(pattern, "The paired rMATS split outputs now replace the historical unpaired row. Both its Tealeaf and comparator direction cells use the same eligible matched universe and union-of-event-BH selection as the other rows. The complete paired catalogue outputs and comparison summaries are in the sibling comparator audit's `rmats` directory; full-data effects are not substituted for split estimates.", description)
        if count != 1:
            raise ValueError("expected one historical rMATS audit description")
    notebook = args.repo_root / "LABNOTEBOOK.md"
    if notebook.exists():
        notes = notebook.read_text()
        tea, comparator = metrics.loc["Tealeaf"], metrics.loc["rMATS (paired)"]
        findings = f"Completed all 380 paired split contrasts. On {int(tea.shared_gene_pairs):,} shared gene--pair hypotheses and {int(tea.shared_genes):,} genes, Tealeaf and rMATS have {int(tea.replicated_bh)} and {int(comparator.replicated_bh)} replicated BH genes, held-half replication {tea.heldout_nominal_replication:.3f} and {comparator.heldout_nominal_replication:.3f}, and split log-p correlation {tea.spearman_logp:.3f} and {comparator.spearman_logp:.3f}. "
        findings += "Direction agreement is " + ", ".join(f"{method} {int(directions.loc[method].agree)}/{int(directions.loc[method].n_direction_evaluable)} among {int(directions.loc[method].n_events)} selected common events" for method in ("Tealeaf", "rMATS (paired)")) + ". The archived full-data long-read analysis remains unchanged."
        notes = replace_once(notes, "Final split metrics are pending, the archived full-data long-read analysis remains unchanged.", findings)
    document.write_text(text)
    if audit_readme.exists():
        audit_readme.write_text(description)
    if notebook.exists():
        notebook.write_text(notes)
    (output / "README.md").write_text("# Paired rMATS subject-split comparison\n\nBoth independent 24-subject halves use the production cell-type contrast manifests and paired PAIRADISE statistical refits. The retained corrected single-end count pass supplies JCEC inclusion/skipping counts and one common event catalogue, without rescanning alignments. Each contrast's sample indices are validated and reordered into identical biological-subject order on both sides. No full-data effects or p-values are substituted for split outputs.\n\n`split_data_fold0_tests.tsv.gz` and `split_data_fold1_tests.tsv.gz` retain all catalogue rows, including events with unavailable native p-values. Finite eligible event tests enter Simes aggregation within gene--pair and then gene, followed by two-fold conjunction and BH on the gene--pair universe shared with current production Tealeaf. Native per-contrast rMATS FDR is retained but is not a screening step for Table 1.\n\n`split_data_tealeaf_comparison_metrics.tsv` gives both methods' matched replication metrics. `split_data_direction_agreement.tsv` uses event BH within each complete matched fold, followed by common feature identities and the union of significance in either fold. Native rMATS delta is A-minus-B and is converted to B-minus-A for alignment. Missing or rounded-zero effects remain selected but are excluded from evaluable denominators. Tealeaf uses verified recovered path-effect vectors.\n\n`split_data_contrast_summary.tsv` and `split_data_manifest.json` record completeness, pairing and provenance. `table1_row.tex` is generated from the matched metrics and directional audit. Historical unpaired aggregate metrics remain separate archived results; the new split fits do not change the independent full-data long-read analysis.\n")
    print(row, flush=True)


if __name__ == "__main__":
    main()
