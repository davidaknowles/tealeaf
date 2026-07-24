#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (!(length(args) %in% c(2, 4))) {
  stop(
    paste(
      "usage: export_seurat_metadata.R INPUT.rds[.gz] OUTPUT.tsv.gz",
      "[EMBEDDING.tsv.gz REDUCTION]"
    )
  )
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

if (length(args) == 4) {
  reduction <- args[[4]]
  if (!(reduction %in% names(object@reductions))) {
    stop(paste("missing reduction", reduction))
  }
  embedding <- object@reductions[[reduction]]@cell.embeddings
  if (is.null(rownames(embedding)) || anyDuplicated(rownames(embedding))) {
    stop("embedding row names are absent or duplicated")
  }
  embedding <- data.frame(
    cell_id = rownames(embedding),
    embedding,
    check.names = FALSE
  )
  connection <- gzfile(args[[3]], open = "wt")
  write.table(
    embedding,
    file = connection,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    col.names = TRUE,
    na = ""
  )
  close(connection)
}
