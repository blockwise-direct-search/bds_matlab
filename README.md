# Blockwise Direct Search (BDS)

BDS is a package for solving nonlinear optimization problems without using derivatives. The current version can handle unconstrained problems. 

## What is BDS?

BDS is a derivative-free package using blockwise direct-search methods. The current version is implemented in MATLAB and it will be implemented in other programming languages in the future.

See [Haitian LI's presentation](https://lht97.github.io/documents/ICNONLA2023.pdf) on BDS for more information.

## How to test BDS?
The tests are **automated** by
[GitHub Actions](https://docs.github.com/en/actions). 

- [Tests](https://github.com/blockwise-direct-search/bds/actions) at [blockwise-direct-search/bds](https://github.com/blockwise-direct-search/bds/)

    - [![Check Spelling](https://github.com/blockwise-direct-search/bds/actions/workflows/spelling.yml/badge.svg)](https://github.com/blockwise-direct-search/bds/actions/workflows/spelling.yml)
    - [![Unit test of BDS](https://github.com/blockwise-direct-search/bds/actions/workflows/unit_test.yml/badge.svg)](https://github.com/blockwise-direct-search/bds/actions/workflows/unit_test.yml)
    - [![Stress test of BDS](https://github.com/blockwise-direct-search/bds/actions/workflows/stress_test.yml/badge.svg)](https://github.com/blockwise-direct-search/bds/actions/workflows/stress_test.yml)
    - [![Test sufficient decrease factor of BDS, small](https://github.com/blockwise-direct-search/bds/actions/workflows/profile_small_sd.yml/badge.svg)](https://github.com/blockwise-direct-search/bds/actions/workflows/profile_small_sd.yml)
    - [![Test sufficient decrease factor of BDS, big](https://github.com/blockwise-direct-search/bds/actions/workflows/profile_big_sd.yml/badge.svg)](https://github.com/blockwise-direct-search/bds/actions/workflows/profile_big_sd.yml)
    - [![Performance profiles of PBDS, big](https://github.com/blockwise-direct-search/bds/actions/workflows/profile_pbds_big.yml/badge.svg)](https://github.com/blockwise-direct-search/bds/actions/workflows/profile_pbds_big.yml)
    - [![Performance profiles of PBDS, small](https://github.com/blockwise-direct-search/bds/actions/workflows/profile_pbds_small.yml/badge.svg)](https://github.com/blockwise-direct-search/bds/actions/workflows/profile_pbds_small.yml)
    - [![Performance profiles of RBDS, big](https://github.com/blockwise-direct-search/bds/actions/workflows/profile_rbds_big.yml/badge.svg)](https://github.com/blockwise-direct-search/bds/actions/workflows/profile_rbds_big.yml)
    - [![Performance profiles of RBDS, small](https://github.com/blockwise-direct-search/bds/actions/workflows/profile_rbds_small.yml/badge.svg)](https://github.com/blockwise-direct-search/bds/actions/workflows/profile_rbds_small.yml)
    - [![Performance profiles of LAM, big](https://github.com/blockwise-direct-search/bds/actions/workflows/profile_lam_big.yml/badge.svg)](https://github.com/blockwise-direct-search/bds/actions/workflows/profile_lam_big.yml)
    - [![Performance profiles of LAM, small](https://github.com/blockwise-direct-search/bds/actions/workflows/profile_lam_small.yml/badge.svg)](https://github.com/blockwise-direct-search/bds/actions/workflows/profile_lam_small.yml)
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