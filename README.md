# Blockwise Direct Search (BDS)

BDS is a package for solving nonlinear optimization problems without using derivatives. The current version can handle unconstrained problems. 

## What is BDS?

BDS is a derivative-free package using blockwise direct-search methods. The current version is implemented in MATLAB, and it is being implemented in other programming languages.

See [Haitian LI's presentation](https://lht97.github.io/documents/DFOS2024.pdf) on BDS for more information.

## How to install BDS?

1. Clone this repository. You should then get a folder named `bds_matlab` containing this README file and the
[`setup.m`](https://github.com/blockwise-direct-search/bds/blob/main/setup.m) file.

2. In the command window of MATLAB, change your directory to the above-mentioned folder, and execute

```matlab
setup
```

If the above succeeds, then the package `bds` is installed and ready to use. Try `help bds` for more information.

We do not support MATLAB R2017a or earlier. If there exists any problems, please open an issue by
https://github.com/blockwise-direct-search/bds_matlab/issues.

## The coverage of unit test (offered by [Codecov](https://about.codecov.io/))

[![Codecov](https://img.shields.io/codecov/c/github/blockwise-direct-search/bds_matlab?style=for-the-badge&logo=codecov)](https://app.codecov.io/github/blockwise-direct-search/bds_matlab)

## Test of BDS.
The tests are **automated** by
[GitHub Actions](https://docs.github.com/en/actions).
- [![Check Spelling](https://github.com/blockwise-direct-search/bds_matlab/actions/workflows/spelling.yml/badge.svg)](https://github.com/blockwise-direct-search/bds_matlab/actions/workflows/spelling.yml)
- [![Unit test of BDS](https://github.com/blockwise-direct-search/bds_matlab/actions/workflows/unit_test.yml/badge.svg)](https://github.com/blockwise-direct-search/bds_matlab/actions/workflows/unit_test.yml)
- [![Coverage test of BDS](https://github.com/blockwise-direct-search/bds_matlab/actions/workflows/unit_test_coverage.yml/badge.svg)](https://github.com/blockwise-direct-search/bds_matlab/actions/workflows/unit_test_coverage.yml)
- [![Stress test of BDS](https://github.com/blockwise-direct-search/bds_matlab/actions/workflows/stress_test.yml/badge.svg)](https://github.com/blockwise-direct-search/bds_matlab/actions/workflows/stress_test.yml)
- [![Parallel test of BDS](https://github.com/blockwise-direct-search/bds_matlab/actions/workflows/parallel_test.yml/badge.svg)](https://github.com/blockwise-direct-search/bds_matlab/actions/workflows/parallel_test.yml)
- [![Recursive test of BDS](https://github.com/blockwise-direct-search/bds_matlab/actions/workflows/recursive_test.yml/badge.svg)](https://github.com/blockwise-direct-search/bds_matlab/actions/workflows/recursive_test.yml)
- [![Test default, small](https://github.com/blockwise-direct-search/bds_matlab/actions/workflows/test_default_small.yml/badge.svg)](https://github.com/blockwise-direct-search/bds_matlab/actions/workflows/test_default_small.yml)
- [![Test default, big](https://github.com/blockwise-direct-search/bds_matlab/actions/workflows/test_default_big.yml/badge.svg)](https://github.com/blockwise-direct-search/bds_matlab/actions/workflows/test_default_big.yml)
- [![Verify norma](https://github.com/blockwise-direct-search/bds_matlab/actions/workflows/verify_norma.yml/badge.svg)](https://github.com/blockwise-direct-search/bds_matlab/actions/workflows/verify_norma.yml)
- [![Test badly_scaled](https://github.com/blockwise-direct-search/bds_matlab/actions/workflows/test_badly_scaled.yml/badge.svg)](https://github.com/blockwise-direct-search/bds_matlab/actions/workflows/test_badly_scaled.yml)
- [![Test StepTolerance](https://github.com/blockwise-direct-search/bds_matlab/actions/workflows/test_StepTolerance.yml/badge.svg)](https://github.com/blockwise-direct-search/bds_matlab/actions/workflows/test_StepTolerance.yml)
- [![Test forcing function](https://github.com/blockwise-direct-search/bds_matlab/actions/workflows/test_forcing_function.yml/badge.svg)](https://github.com/blockwise-direct-search/bds_matlab/actions/workflows/test_forcing_function.yml)
- [![Test replacement_delay, small](https://github.com/blockwise-direct-search/bds_matlab/actions/workflows/test_replacement_delay_small.yml/badge.svg)](https://github.com/blockwise-direct-search/bds_matlab/actions/workflows/test_replacement_delay_small.yml)
- [![Test replacement_delay, big](https://github.com/blockwise-direct-search/bds_matlab/actions/workflows/test_replacement_delay_big.yml/badge.svg)](https://github.com/blockwise-direct-search/bds_matlab/actions/workflows/test_replacement_delay_big.yml)
- [![Test the technique of initial step size using initial point](https://github.com/blockwise-direct-search/bds_matlab/actions/workflows/test_initial_point_technique.yml/badge.svg)](https://github.com/blockwise-direct-search/bds_matlab/actions/workflows/test_initial_point_technique.yml)
- [![Test reduction factor of BDS, small](https://github.com/blockwise-direct-search/bds_matlab/actions/workflows/profile_small_reduction_factor.yml/badge.svg)](https://github.com/blockwise-direct-search/bds_matlab/actions/workflows/profile_small_reduction_factor.yml)
- [![Test reduction factor of BDS, big](https://github.com/blockwise-direct-search/bds_matlab/actions/workflows/profile_big_reduction_factor.yml/badge.svg)](https://github.com/blockwise-direct-search/bds_matlab/actions/workflows/profile_big_reduction_factor.yml)
- [![Performance profiles of PADS, small](https://github.com/blockwise-direct-search/bds_matlab/actions/workflows/profile_pads_small.yml/badge.svg)](https://github.com/blockwise-direct-search/bds_matlab/actions/workflows/profile_pads_small.yml)
- [![Performance profiles of PADS, big](https://github.com/blockwise-direct-search/bds_matlab/actions/workflows/profile_pads_big.yml/badge.svg)](https://github.com/blockwise-direct-search/bds_matlab/actions/workflows/profile_pads_big.yml)
- [![Performance profiles of scbds, small](https://github.com/blockwise-direct-search/bds_matlab/actions/workflows/profile_scbds_small.yml/badge.svg)](https://github.com/blockwise-direct-search/bds_matlab/actions/workflows/profile_scbds_small.yml)
- [![Performance profiles of scbds, big](https://github.com/blockwise-direct-search/bds_matlab/actions/workflows/profile_scbds_big.yml/badge.svg)](https://github.com/blockwise-direct-search/bds_matlab/actions/workflows/profile_scbds_big.yml)
- [![Performance profiles of PBDS, small](https://github.com/blockwise-direct-search/bds_matlab/actions/workflows/profile_pbds_small.yml/badge.svg)](https://github.com/blockwise-direct-search/bds_matlab/actions/workflows/profile_pbds_small.yml)
- [![Performance profiles of PBDS, big](https://github.com/blockwise-direct-search/bds_matlab/actions/workflows/profile_pbds_big.yml/badge.svg)](https://github.com/blockwise-direct-search/bds_matlab/actions/workflows/profile_pbds_big.yml)
- [![Performance profiles of RBDS, small](https://github.com/blockwise-direct-search/bds_matlab/actions/workflows/profile_rbds_small.yml/badge.svg)](https://github.com/blockwise-direct-search/bds_matlab/actions/workflows/profile_rbds_small.yml)
- [![Performance profiles of RBDS, big](https://github.com/blockwise-direct-search/bds_matlab/actions/workflows/profile_rbds_big.yml/badge.svg)](https://github.com/blockwise-direct-search/bds_matlab/actions/workflows/profile_rbds_big.yml)
- [![Performance profiles of lam, big](https://github.com/blockwise-direct-search/bds_matlab/actions/workflows/profile_lam_big.yml/badge.svg)](https://github.com/blockwise-direct-search/bds_matlab/actions/workflows/profile_lam_big.yml)
- [![Performance profiles of lam, small](https://github.com/blockwise-direct-search/bds_matlab/actions/workflows/profile_lam_small.yml/badge.svg)](https://github.com/blockwise-direct-search/bds_matlab/actions/workflows/profile_lam_small.yml)
- [![Performance profiles of BFGS, big](https://github.com/blockwise-direct-search/bds_matlab/actions/workflows/profile_bfgs_big.yml/badge.svg)](https://github.com/blockwise-direct-search/bds_matlab/actions/workflows/profile_bfgs_big.yml)
- [![Performance profiles of BFGS, small](https://github.com/blockwise-direct-search/bds_matlab/actions/workflows/profile_bfgs_small.yml/badge.svg)](https://github.com/blockwise-direct-search/bds_matlab/actions/workflows/profile_bfgs_small.yml)
- [![Performance profiles of PRIMA, big](https://github.com/blockwise-direct-search/bds_matlab/actions/workflows/profile_prima_big.yml/badge.svg)](https://github.com/blockwise-direct-search/bds_matlab/actions/workflows/profile_prima_big.yml)
- [![Performance profiles of PRIMA, small](https://github.com/blockwise-direct-search/bds_matlab/actions/workflows/profile_prima_small.yml/badge.svg)](https://github.com/blockwise-direct-search/bds_matlab/actions/workflows/profile_prima_small.yml)
- [![Performance profiles of simplex, big](https://github.com/blockwise-direct-search/bds_matlab/actions/workflows/profile_simplex_big.yml/badge.svg)](https://github.com/blockwise-direct-search/bds_matlab/actions/workflows/profile_simplex_big.yml)
- [![Performance profiles of simplex, small](https://github.com/blockwise-direct-search/bds_matlab/actions/workflows/profile_simplex_small.yml/badge.svg)](https://github.com/blockwise-direct-search/bds_matlab/actions/workflows/profile_simplex_small.yml)
- [![Performance profiles of NLOPT, big](https://github.com/blockwise-direct-search/bds_matlab/actions/workflows/profile_nlopt_big.yml/badge.svg)](https://github.com/blockwise-direct-search/bds_matlab/actions/workflows/profile_nlopt_big.yml)
- [![Performance profiles of NLOPT, small](https://github.com/blockwise-direct-search/bds_matlab/actions/workflows/profile_nlopt_small.yml/badge.svg)](https://github.com/blockwise-direct-search/bds_matlab/actions/workflows/profile_nlopt_small.yml)
- [![Performance profiles of BFO, big](https://github.com/blockwise-direct-search/bds_matlab/actions/workflows/profile_bfo_big.yml/badge.svg)](https://github.com/blockwise-direct-search/bds_matlab/actions/workflows/profile_bfo_big.yml)
- [![Performance profiles of BFO, small](https://github.com/blockwise-direct-search/bds_matlab/actions/workflows/profile_bfo_small.yml/badge.svg)](https://github.com/blockwise-direct-search/bds_matlab/actions/workflows/profile_bfo_small.yml)
