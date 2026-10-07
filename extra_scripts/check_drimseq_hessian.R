# Reference implementation diagnostic, optional first argument is an R library.
args <- commandArgs(trailingOnly = TRUE)
if (length(args)) .libPaths(c(args[1], .libPaths()))
library(DRIMSeq)
x <- cbind(intercept = 1, group = rep(c(0, 1), each = 6))
y <- rbind(path1 = c(40, 45, 38, 44, 42, 39, 60, 65, 58, 64, 62, 59),
           path2 = 100 - c(40, 45, 38, 44, 42, 39, 60, 65, 58, 64, 62, 59))
colnames(y) <- paste0("sample", seq_len(ncol(y)))
fit <- DRIMSeq:::dm_fitRegression(y, x, prec = 20)
parameters <- c(t(fit$b[-nrow(fit$b), , drop = FALSE]))
analytic <- DRIMSeq:::dm_Hessian_regG_prop(y = t(y), prec = 20, prop = t(fit$fit), x = x)
numeric <- optimHess(parameters, DRIMSeq:::dm_lik_regG, DRIMSeq:::dm_score_regG, x = x, prec = 20, y = y)
print(packageVersion("DRIMSeq"))
print(list(analytic = analytic, numeric = numeric,
           analytic_information_eigenvalues = eigen(-analytic, symmetric = TRUE)$values,
           numeric_information_eigenvalues = eigen(-numeric, symmetric = TRUE)$values,
           stock_adjustment = DRIMSeq:::dm_CRadjustmentRegression(y, x, 20, fit$fit),
           numeric_adjustment = as.numeric(determinant(ncol(y) * (-numeric), logarithm = TRUE)$modulus) / 2))
