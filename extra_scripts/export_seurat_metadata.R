#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) {
  stop("usage: export_seurat_metadata.R INPUT.rds[.gz] OUTPUT.tsv.gz")
}

input <- if (grepl("\\.gz$", args[[1]], ignore.case = TRUE)) {
  outer <- gzfile(args[[1]], open = "rb")
  magic <- readBin(outer, what = "raw", n = 2)
  close(outer)
  outer <- gzfile(args[[1]], open = "rb")
  if (identical(magic, as.raw(c(0x1f, 0x8b)))) gzcon(outer) else outer
} else {
  args[[1]]
}
object <- readRDS(input)
if (inherits(input, "connection")) {
  close(input)
}
metadata <- object[[]]
metadata <- cbind(cell_id = rownames(metadata), metadata)
connection <- gzfile(args[[2]], open = "wt")
write.table(
  metadata,
  file = connection,
  sep = "\t",
  quote = TRUE,
  row.names = FALSE,
  col.names = TRUE,
  na = ""
)
close(connection)
