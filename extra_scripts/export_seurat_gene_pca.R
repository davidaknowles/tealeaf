#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (!(length(args) %in% c(2, 3, 4))) {
  stop("usage: export_seurat_gene_pca.R INPUT.rds[.gz] OUTPUT.tsv.gz [NFEATURES] [NPCS]")
}

nfeatures <- if (length(args) >= 3) as.integer(args[[3]]) else 2000L
npcs <- if (length(args) >= 4) as.integer(args[[4]]) else 30L
if (is.na(nfeatures) || nfeatures < 2L || is.na(npcs) || npcs < 2L) {
  stop("NFEATURES and NPCS must be integers of at least two")
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

suppressPackageStartupMessages(library(Matrix))
assay <- object@assays[["RNA"]]
features <- head(assay@var.features, nfeatures)
if (length(features) < 2L) {
  means <- Matrix::rowMeans(assay@data)
  second_moments <- Matrix::rowMeans(assay@data ^ 2)
  variances <- (ncol(assay@data) / (ncol(assay@data) - 1)) *
    pmax(second_moments - means ^ 2, 0)
  features <- rownames(assay@data)[
    head(order(variances, decreasing = TRUE), nfeatures)
  ]
}
expression <- Matrix::t(assay@data[features, , drop = FALSE])
centers <- Matrix::colMeans(expression)
second_moments <- Matrix::colMeans(expression ^ 2)
scales <- sqrt(
  pmax(
    (nrow(expression) / (nrow(expression) - 1)) *
      (second_moments - centers ^ 2),
    0
  )
)
scales[!is.finite(scales) | scales == 0] <- 1

rank <- min(npcs, ncol(expression) - 1L)
working_rank <- min(rank + 10L, ncol(expression))
set.seed(0)
omega <- matrix(rnorm(ncol(expression) * working_rank), ncol(expression))
center_over_scale <- centers / scales
right_multiply <- function(vectors) {
  scaled <- sweep(vectors, 1, scales, `/`)
  as.matrix(expression %*% scaled) -
    outer(rep(1, nrow(expression)), as.vector(centers %*% scaled))
}
left_multiply <- function(vectors) {
  cross <- as.matrix(Matrix::crossprod(expression, vectors))
  sweep(cross, 1, scales, `/`) -
    outer(center_over_scale, colSums(vectors))
}

q <- qr.Q(qr(right_multiply(omega)))
for (iteration in seq_len(2L)) {
  q <- qr.Q(qr(right_multiply(left_multiply(q))))
}
small <- t(left_multiply(q))
decomposition <- svd(small, nu = rank, nv = 0L)
embedding <- sweep(
  q %*% decomposition$u[, seq_len(rank), drop = FALSE],
  2,
  decomposition$d[seq_len(rank)],
  `*`
)
rownames(embedding) <- rownames(expression)
embedding <- data.frame(
  cell_id = rownames(embedding),
  embedding,
  check.names = FALSE
)
connection <- gzfile(args[[2]], open = "wt")
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
