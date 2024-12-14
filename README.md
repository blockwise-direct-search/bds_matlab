# Blockwise Direct Search (BDS)

BDS is a package for solving nonlinear optimization problems without using derivatives. The current version can handle unconstrained problems. 

## What is BDS?

BDS is a derivative-free package using blockwise direct-search methods. The current version is implemented in MATLAB, and it is being implemented in other programming languages.

See [Haitian LI's presentation](https://lht97.github.io/documents/DFOS2024.pdf) on BDS for more information.

## How to install BDS?

1. Clone this repository. You should then get a folder named `bds` containing this README file and the
[`setup.m`](https://github.com/blockwise-direct-search/bds/blob/main/setup.m) file.

2. In the command window of MATLAB, change your directory to the above-mentioned folder, and execute

```matlab
setup
```

If the above succeeds, then the package `bds` is installed and ready to use. Try `help bds` for more information.

We do not support MATLAB R2017a or earlier. If there exists any problems, please open an issue by
https://github.com/blockwise-direct-search/bds/issues.

## The coverage of unit test (offered by [Codecov](https://about.codecov.io/))

[![Codecov](https://img.shields.io/codecov/c/github/blockwise-direct-search/bds?style=for-the-badge&logo=codecov)](https://app.codecov.io/github/blockwise-direct-search/bds)

## Test of BDS.
The tests are **automated** by
[GitHub Actions](https://docs.github.com/en/actions).
- [![Check Spelling](https://github.com/blockwise-direct-search/bds/actions/workflows/spelling.yml/badge.svg)](https://github.com/blockwise-direct-search/bds/actions/workflows/spelling.yml)
- [![Unit test of BDS](https://github.com/blockwise-direct-search/bds/actions/workflows/unit_test.yml/badge.svg)](https://github.com/blockwise-direct-search/bds/actions/workflows/unit_test.yml)
- [![Coverage test of BDS](https://github.com/blockwise-direct-search/bds/actions/workflows/unit_test_coverage.yml/badge.svg)](https://github.com/blockwise-direct-search/bds/actions/workflows/unit_test_coverage.yml)
- [![Stress test of BDS](https://github.com/blockwise-direct-search/bds/actions/workflows/stress_test.yml/badge.svg)](https://github.com/blockwise-direct-search/bds/actions/workflows/stress_test.yml)
- [![Parallel test of BDS](https://github.com/blockwise-direct-search/bds/actions/workflows/parallel_test.yml/badge.svg)](https://github.com/blockwise-direct-search/bds/actions/workflows/parallel_test.yml)
- [![Recursive test of BDS](https://github.com/blockwise-direct-search/bds/actions/workflows/recursive_test.yml/badge.svg)](https://github.com/blockwise-direct-search/bds/actions/workflows/recursive_test.yml)
- [![Verify norma](https://github.com/blockwise-direct-search/bds/actions/workflows/verify_norma.yml/badge.svg)](https://github.com/blockwise-direct-search/bds/actions/workflows/verify_norma.yml)
- [![Performance profiles of cbds and default cbds, big](https://github.com/blockwise-direct-search/bds/actions/workflows/profile_cbds_default_cbds_big.yml/badge.svg)](https://github.com/blockwise-direct-search/bds/actions/workflows/profile_cbds_default_cbds_big.yml)
- [![Performance profiles of cbds and default cbds, small](https://github.com/blockwise-direct-search/bds/actions/workflows/profile_cbds_default_cbds_small.yml/badge.svg)](https://github.com/blockwise-direct-search/bds/actions/workflows/profile_cbds_default_cbds_small.yml)
- [![Performance profiles of cbds and ds, big](https://github.com/blockwise-direct-search/bds/actions/workflows/profile_cbds_ds_big.yml/badge.svg)](https://github.com/blockwise-direct-search/bds/actions/workflows/profile_cbds_ds_big.yml)
- [![Performance profiles of cbds and ds, small](https://github.com/blockwise-direct-search/bds/actions/workflows/profile_cbds_ds_small.yml/badge.svg)](https://github.com/blockwise-direct-search/bds/actions/workflows/profile_cbds_ds_small.yml)
- [![Performance profiles of default cbds and ds, big](https://github.com/blockwise-direct-search/bds/actions/workflows/profile_default_cbds_ds_big.yml/badge.svg)](https://github.com/blockwise-direct-search/bds/actions/workflows/profile_default_cbds_ds_big.yml)
- [![Performance profiles of default cbds and ds, small](https://github.com/blockwise-direct-search/bds/actions/workflows/profile_default_cbds_ds_small.yml/badge.svg)](https://github.com/blockwise-direct-search/bds/actions/workflows/profile_default_cbds_ds_small.yml)
- [![Performance profiles of cbds and lam, big](https://github.com/blockwise-direct-search/bds/actions/workflows/profile_cbds_lam_big.yml/badge.svg)](https://github.com/blockwise-direct-search/bds/actions/workflows/profile_cbds_lam_big.yml)
- [![Performance profiles of cbds and lam, small](https://github.com/blockwise-direct-search/bds/actions/workflows/profile_cbds_lam_small.yml/badge.svg)](https://github.com/blockwise-direct-search/bds/actions/workflows/profile_cbds_lam_small.yml)
- [![Performance profiles of default cbds and lam, big](https://github.com/blockwise-direct-search/bds/actions/workflows/profile_default_cbds_lam_big.yml/badge.svg)](https://github.com/blockwise-direct-search/bds/actions/workflows/profile_default_cbds_lam_big.yml)
- [![Performance profiles of default cbds and lam, small](https://github.com/blockwise-direct-search/bds/actions/workflows/profile_default_cbds_lam_small.yml/badge.svg)](https://github.com/blockwise-direct-search/bds/actions/workflows/profile_default_cbds_lam_small.yml)
- [![Performance profiles of cbds and BFGS, big](https://github.com/blockwise-direct-search/bds/actions/workflows/profile_cbds_bfgs_big.yml/badge.svg)](https://github.com/blockwise-direct-search/bds/actions/workflows/profile_cbds_bfgs_big.yml)
- [![Performance profiles of cbds and BFGS, small](https://github.com/blockwise-direct-search/bds/actions/workflows/profile_cbds_bfgs_small.yml/badge.svg)](https://github.com/blockwise-direct-search/bds/actions/workflows/profile_cbds_bfgs_small.yml)
- [![Performance profiles of default cbds and BFGS, big](https://github.com/blockwise-direct-search/bds/actions/workflows/profile_default_cbds_bfgs_big.yml/badge.svg)](https://github.com/blockwise-direct-search/bds/actions/workflows/profile_default_cbds_bfgs_big.yml)
- [![Performance profiles of default cbds and BFGS, small](https://github.com/blockwise-direct-search/bds/actions/workflows/profile_default_cbds_bfgs_small.yml/badge.svg)](https://github.com/blockwise-direct-search/bds/actions/workflows/profile_default_cbds_bfgs_small.yml)
- [![Performance profiles of cbds and PRIMA, big](https://github.com/blockwise-direct-search/bds/actions/workflows/profile_cbds_newuoa_big.yml/badge.svg)](https://github.com/blockwise-direct-search/bds/actions/workflows/profile_cbds_newuoa_big.yml)
- [![Performance profiles of cbds and PRIMA, small](https://github.com/blockwise-direct-search/bds/actions/workflows/profile_cbds_newuoa_small.yml/badge.svg)](https://github.com/blockwise-direct-search/bds/actions/workflows/profile_cbds_newuoa_small.yml)
- [![Performance profiles of default cbds and PRIMA, big](https://github.com/blockwise-direct-search/bds/actions/workflows/profile_default_cbds_newuoa_big.yml/badge.svg)](https://github.com/blockwise-direct-search/bds/actions/workflows/profile_default_cbds_newuoa_big.yml)
- [![Performance profiles of default cbds and PRIMA, small](https://github.com/blockwise-direct-search/bds/actions/workflows/profile_default_cbds_newuoa_small.yml/badge.svg)](https://github.com/blockwise-direct-search/bds/actions/workflows/profile_default_cbds_newuoa_small.yml)
- [![Performance profiles of cbds and simplex, big](https://github.com/blockwise-direct-search/bds/actions/workflows/profile_cbds_simplex_big.yml/badge.svg)](https://github.com/blockwise-direct-search/bds/actions/workflows/profile_cbds_simplex_big.yml)
- [![Performance profiles of cbds and simplex, small](https://github.com/blockwise-direct-search/bds/actions/workflows/profile_cbds_simplex_small.yml/badge.svg)](https://github.com/blockwise-direct-search/bds/actions/workflows/profile_cbds_simplex_small.yml)
- [![Performance profiles of default cbds and simplex, big](https://github.com/blockwise-direct-search/bds/actions/workflows/profile_default_cbds_simplex_big.yml/badge.svg)](https://github.com/blockwise-direct-search/bds/actions/workflows/profile_default_cbds_simplex_big.yml)
- [![Performance profiles of default cbds and simplex, small](https://github.com/blockwise-direct-search/bds/actions/workflows/profile_default_cbds_simplex_small.yml/badge.svg)](https://github.com/blockwise-direct-search/bds/actions/workflows/profile_default_cbds_simplex_small.yml)
- [![Performance profiles of NLOPT, big](https://github.com/blockwise-direct-search/bds/actions/workflows/profile_nlopt_big.yml/badge.svg)](https://github.com/blockwise-direct-search/bds/actions/workflows/profile_nlopt_big.yml)
- [![Performance profiles of NLOPT, small](https://github.com/blockwise-direct-search/bds/actions/workflows/profile_nlopt_small.yml/badge.svg)](https://github.com/blockwise-direct-search/bds/actions/workflows/profile_nlopt_small.yml)
- [![Performance profiles of cbds and bfo, big](https://github.com/blockwise-direct-search/bds/actions/workflows/profile_cbds_bfo_big.yml/badge.svg)](https://github.com/blockwise-direct-search/bds/actions/workflows/profile_cbds_bfo_big.yml)
- [![Performance profiles of cbds and bfo, small](https://github.com/blockwise-direct-search/bds/actions/workflows/profile_cbds_bfo_small.yml/badge.svg)](https://github.com/blockwise-direct-search/bds/actions/workflows/profile_cbds_bfo_small.yml)
