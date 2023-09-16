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
    - [![Performance profiles of BDS with noise, big](https://github.com/blockwise-direct-search/bds/actions/workflows/profile_big_noise.yml/badge.svg)](https://github.com/blockwise-direct-search/bds/actions/workflows/profile_big_noise.yml)
    - [![Performance profiles of BDS with noise, small](https://github.com/blockwise-direct-search/bds/actions/workflows/profile_small_noise.yml/badge.svg)](https://github.com/blockwise-direct-search/bds/actions/workflows/profile_small_noise.yml)
    - [![Performance profiles of BDS, big](https://github.com/blockwise-direct-search/bds/actions/workflows/profile_big.yml/badge.svg)](https://github.com/blockwise-direct-search/bds/actions/workflows/profile_big.yml)
    - [![Performance profiles of BDS, small](https://github.com/blockwise-direct-search/bds/actions/workflows/profile_small.yml/badge.svg)](https://github.com/blockwise-direct-search/bds/actions/workflows/profile_small.yml)
    - [![Stress test of BDS](https://github.com/blockwise-direct-search/bds/actions/workflows/stress_test.yml/badge.svg)](https://github.com/blockwise-direct-search/bds/actions/workflows/stress_test.yml)
    - [![Performance profiles of BDS and BFGS, with noise](https://github.com/blockwise-direct-search/bds/actions/workflows/profile_bfgs_noise.yml/badge.svg)](https://github.com/blockwise-direct-search/bds/actions/workflows/profile_bfgs_noise.yml)
    - [![Performance profiles of BDS and BFGS](https://github.com/blockwise-direct-search/bds/actions/workflows/profile_bfgs.yml/badge.svg)](https://github.com/blockwise-direct-search/bds/actions/workflows/profile_bfgs.yml)
    - [![Performance profiles of BDS and PRIMA, with noise](https://github.com/blockwise-direct-search/bds/actions/workflows/profile_prima_noise.yml/badge.svg)](https://github.com/blockwise-direct-search/bds/actions/workflows/profile_prima_noise.yml)
    - [![Performance profiles of BDS and PRIMA](https://github.com/blockwise-direct-search/bds/actions/workflows/profile_prima.yml/badge.svg)](https://github.com/blockwise-direct-search/bds/actions/workflows/profile_prima.yml)
    - [![Performance profiles of BDS and simplex, with noise](https://github.com/blockwise-direct-search/bds/actions/workflows/profile_simplex_noise.yml/badge.svg)](https://github.com/blockwise-direct-search/bds/actions/workflows/profile_simplex_noise.yml)
    - [![Performance profiles of BDS and simplex](https://github.com/blockwise-direct-search/bds/actions/workflows/profile_simplex.yml/badge.svg)](https://github.com/blockwise-direct-search/bds/actions/workflows/profile_simplex.yml)
    - [![Performance profiles of BDS and NLOPT, with noise](https://github.com/blockwise-direct-search/bds/actions/workflows/profile_nlopt_noise.yml/badge.svg)](https://github.com/blockwise-direct-search/bds/actions/workflows/profile_nlopt_noise.yml)
    - [![Performance profiles of BDS and NLOPT](https://github.com/blockwise-direct-search/bds/actions/workflows/profile_nlopt.yml/badge.svg)](https://github.com/blockwise-direct-search/bds/actions/workflows/profile_nlopt.yml)