#!/usr/bin/env python3
"""Run the identical EC pilot with prior-only Gauss-Jacobi quadrature."""

from extra_scripts.audit_binary_marginal_ec import main
from tealeaf.sc.path_marginal_quadrature import prior_quadrature


if __name__ == "__main__":
    main(likelihood_transform=prior_quadrature, quadrature_backend="prior_gauss_jacobi")
