# Omnibus cell-type event annotation

`event_catalog.tsv` contains one row per BH-significant block from the production subject-blocked omnibus cell-type analysis after testing-level EB selected fixed-total concentration \(A=64\). `event_type_summary.tsv` and `event_type_block_summary.tsv` classify the 791 discoveries from their supported local paths. Because the omnibus family has one hypothesis per block, the two summaries are identical apart from their count-column names.

`literature_matches.tsv` intersects the omnibus discoveries with the nine events in the literature audit at `analyses/celltype_event_annotation/literature_examples.tsv`. Eight pass omnibus BH, `Ntrk2`, `Grin1`, `Dlg4`, `App`, `Gria2`, `Scn2a`, `Rbfox3`, and `Mef2d`; `Enah` does not. Pairwise results are not used to select the event classes or literature examples, although they remain available as a secondary localization analysis.
