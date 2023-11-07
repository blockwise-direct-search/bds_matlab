# Blockwise Direct Search (BDS)

BDS is a package for solving nonlinear optimization problems without using derivatives. The current version can handle unconstrained problems. 

## What is BDS?

BDS is a derivative-free package using blockwise direct-search methods. The current version is implemented in MATLAB and it will be implemented in other programming languages in the future.

See [Haitian LI's presentation](https://lht97.github.io/documents/ICNONLA2023.pdf) on BDS for more information.

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

## Test of BDS.
The tests are **automated** by
[GitHub Actions](https://docs.github.com/en/actions). 
- [![Check Spelling](https://github.com/blockwise-direct-search/bds/actions/workflows/spelling.yml/badge.svg)](https://github.com/blockwise-direct-search/bds/actions/workflows/spelling.yml)
- [![Unit test of BDS](https://github.com/blockwise-direct-search/bds/actions/workflows/unit_test.yml/badge.svg)](https://github.com/blockwise-direct-search/bds/actions/workflows/unit_test.yml)
- [![Stress test of BDS](https://github.com/blockwise-direct-search/bds/actions/workflows/stress_test.yml/badge.svg)](https://github.com/blockwise-direct-search/bds/actions/workflows/stress_test.yml)
- [![Parallel test of BDS](https://github.com/blockwise-direct-search/bds/actions/workflows/parallel_test_matlab.yml/badge.svg)](https://github.com/blockwise-direct-search/bds/actions/workflows/parallel_test_matlab.yml)
- [![Recursive test of BDS](https://github.com/blockwise-direct-search/bds/actions/workflows/recursive_test_matlab.yml/badge.svg)](https://github.com/blockwise-direct-search/bds/actions/workflows/recursive_test_matlab.yml)
- [![Test StepTolerance](https://github.com/blockwise-direct-search/bds/actions/workflows/test_StepTolerance.yml/badge.svg)](https://github.com/blockwise-direct-search/bds/actions/workflows/test_StepTolerance.yml)
- [![Test forcing function](https://github.com/blockwise-direct-search/bds/actions/workflows/test_forcing_function.yml/badge.svg)](https://github.com/blockwise-direct-search/bds/actions/workflows/test_forcing_function.yml)
- [![Test the technique of initial step size using initial point](https://github.com/blockwise-direct-search/bds/actions/workflows/test_initial_point_technique.yml/badge.svg)](https://github.com/blockwise-direct-search/bds/actions/workflows/test_initial_point_technique.yml)
- [![Test reduction factor of BDS, small](https://github.com/blockwise-direct-search/bds/actions/workflows/profile_small_reduction_factor.yml/badge.svg)](https://github.com/blockwise-direct-search/bds/actions/workflows/profile_small_reduction_factor.yml)
- [![Test reduction factor of BDS, big](https://github.com/blockwise-direct-search/bds/actions/workflows/profile_big_reduction_factor.yml/badge.svg)](https://github.com/blockwise-direct-search/bds/actions/workflows/profile_big_reduction_factor.yml)
- [![Performance profiles of PBDS, big](https://github.com/blockwise-direct-search/bds/actions/workflows/profile_pbds_big.yml/badge.svg)](https://github.com/blockwise-direct-search/bds/actions/workflows/profile_pbds_big.yml)
- [![Performance profiles of PBDS, small](https://github.com/blockwise-direct-search/bds/actions/workflows/profile_pbds_small.yml/badge.svg)](https://github.com/blockwise-direct-search/bds/actions/workflows/profile_pbds_small.yml)
- [![Performance profiles of RBDS, small](https://github.com/blockwise-direct-search/bds/actions/workflows/profile_rbds_small.yml/badge.svg)](https://github.com/blockwise-direct-search/bds/actions/workflows/profile_rbds_small.yml)
- [![Performance profiles of lam, big](https://github.com/blockwise-direct-search/bds/actions/workflows/profile_lam_big.yml/badge.svg)](https://github.com/blockwise-direct-search/bds/actions/workflows/profile_lam_big.yml)
- [![Performance profiles of lam, small](https://github.com/blockwise-direct-search/bds/actions/workflows/profile_lam_small.yml/badge.svg)](https://github.com/blockwise-direct-search/bds/actions/workflows/profile_lam_small.yml)
- [![Performance profiles of DSPD, big](https://github.com/blockwise-direct-search/bds/actions/workflows/profile_dspd_big.yml/badge.svg)](https://github.com/blockwise-direct-search/bds/actions/workflows/profile_dspd_big.yml)
- [![Performance profiles of DSPD, small](https://github.com/blockwise-direct-search/bds/actions/workflows/profile_dspd_small.yml/badge.svg)](https://github.com/blockwise-direct-search/bds/actions/workflows/profile_dspd_small.yml)
- [![Performance profiles of BFGS, big](https://github.com/blockwise-direct-search/bds/actions/workflows/profile_bfgs_big.yml/badge.svg)](https://github.com/blockwise-direct-search/bds/actions/workflows/profile_bfgs_big.yml)
- [![Performance profiles of BFGS, small](https://github.com/blockwise-direct-search/bds/actions/workflows/profile_bfgs_small.yml/badge.svg)](https://github.com/blockwise-direct-search/bds/actions/workflows/profile_bfgs_small.yml)
- [![Performance profiles of PRIMA, big](https://github.com/blockwise-direct-search/bds/actions/workflows/profile_prima_big.yml/badge.svg)](https://github.com/blockwise-direct-search/bds/actions/workflows/profile_prima_big.yml)
- [![Performance profiles of PRIMA, small](https://github.com/blockwise-direct-search/bds/actions/workflows/profile_prima_small.yml/badge.svg)](https://github.com/blockwise-direct-search/bds/actions/workflows/profile_prima_small.yml)
- [![Performance profiles of simplex, big](https://github.com/blockwise-direct-search/bds/actions/workflows/profile_simplex_big.yml/badge.svg)](https://github.com/blockwise-direct-search/bds/actions/workflows/profile_simplex_big.yml)
- [![Performance profiles of simplex, small](https://github.com/blockwise-direct-search/bds/actions/workflows/profile_simplex_small.yml/badge.svg)](https://github.com/blockwise-direct-search/bds/actions/workflows/profile_simplex_small.yml)
- [![Performance profiles of NLOPT, big](https://github.com/blockwise-direct-search/bds/actions/workflows/profile_nlopt_big.yml/badge.svg)](https://github.com/blockwise-direct-search/bds/actions/workflows/profile_nlopt_big.yml)
- [![Performance profiles of NLOPT, small](https://github.com/blockwise-direct-search/bds/actions/workflows/profile_nlopt_small.yml/badge.svg)](https://github.com/blockwise-direct-search/bds/actions/workflows/profile_nlopt_small.yml)
- [![Performance profiles of BFO, big](https://github.com/blockwise-direct-search/bds/actions/workflows/profile_bfo_big.yml/badge.svg)](https://github.com/blockwise-direct-search/bds/actions/workflows/profile_bfo_big.yml)
- [![Performance profiles of BFO, small](https://github.com/blockwise-direct-search/bds/actions/workflows/profile_bfo_small.yml/badge.svg)](https://github.com/blockwise-direct-search/bds/actions/workflows/profile_bfo_small.yml)