# Omnibus cell-type event annotation

`event_catalog.tsv` contains one row per BH-significant block from the production subject-blocked omnibus cell-type analysis. `event_type_summary.tsv` and `event_type_block_summary.tsv` classify the 584 discoveries from their supported local paths. Because the omnibus family has one hypothesis per block, the two summaries are identical apart from their count-column names.

`literature_matches.tsv` intersects the omnibus discoveries with the five events in the pre-specified literature audit at `analyses/celltype_event_annotation/literature_examples.tsv`. Four pass omnibus BH, `Ntrk2`, `Grin1`, `Dlg4`, and `App`; `Enah` does not. Pairwise results are not used to select the event classes or literature examples, although they remain available as a secondary localization analysis.
