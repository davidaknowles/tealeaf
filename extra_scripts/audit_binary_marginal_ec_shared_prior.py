#!/usr/bin/env python3
"""Run the identical EC pilot with shared prior quadrature/probability grids."""

from extra_scripts.audit_binary_marginal_ec import main
from tealeaf.sc.path_marginal_quadrature import shared_prior_quadrature


if __name__ == "__main__":
    main(likelihood_transform=shared_prior_quadrature, quadrature_backend="shared_prior_gauss_jacobi")
