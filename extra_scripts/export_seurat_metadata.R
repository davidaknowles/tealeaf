#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) {
  stop("usage: export_seurat_metadata.R INPUT.rds[.gz] OUTPUT.tsv.gz")
}

object <- readRDS(args[[1]])
metadata <- object[[]]
metadata <- cbind(cell_id = rownames(metadata), metadata)
connection <- gzfile(args[[2]], open = "wt")
write.table(
  metadata,
  file = connection,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = TRUE,
  na = ""
)
close(connection)
