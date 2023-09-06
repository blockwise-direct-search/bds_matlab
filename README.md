# Blockwise Direct Search (BDS)

BDS is a package for solving nonlinear optimization problems without using derivatives. The current version can handle unconstrained problems. 

### What is BDS?

BDS is a derivative-free package using blockwise direct-search methods. The current version is implemented in MATLAB and it will be implemented in other programming languages in the future.

See [Haitian LI's talk](https://lht97.github.io/documents/ICNONLA2023.pdf) on BDS at [The 14th International Conference on Numerical Optimization and Numerical Linear Algebra](http://lsec.cc.ac.cn/~icnonla23/) for more information.

### How to test BDS?
The tests are **automated** by
[GitHub Actions](https://docs.github.com/en/actions). 

- [Tests](https://github.com/blockwise-direct-search/bds/actions) at [blockwise-direct-search/bds](https://github.com/blockwise-direct-search/bds/)

    - [![Check Spelling](https://github.com/blockwise-direct-search/bds/actions/workflows/spelling.yml/badge.svg)](https://github.com/blockwise-direct-search/bds/actions/workflows/spelling.yml)
    - [![Performance profiles of BDS with noise, big](https://github.com/blockwise-direct-search/bds/actions/workflows/test_bds_big_with_noise.yml/badge.svg)](https://github.com/blockwise-direct-search/bds/actions/workflows/test_bds_big_with_noise.yml)
    - [![Performance profiles of BDS with noise, small](https://github.com/blockwise-direct-search/bds/actions/workflows/test_bds_small_with_noise.yml/badge.svg)](https://github.com/blockwise-direct-search/bds/actions/workflows/test_bds_small_with_noise.yml)
    - [![Performance profiles of BDS, big](https://github.com/blockwise-direct-search/bds/actions/workflows/test_bds_big.yml/badge.svg)](https://github.com/blockwise-direct-search/bds/actions/workflows/test_bds_big.yml)
    - [![Performance profiles of BDS, small](https://github.com/blockwise-direct-search/bds/actions/workflows/test_bds_small.yml/badge.svg)](https://github.com/blockwise-direct-search/bds/actions/workflows/test_bds_small.yml)
    - [![Stress test of BDS](https://github.com/blockwise-direct-search/bds/actions/workflows/stress_test.yml/badge.svg)](https://github.com/blockwise-direct-search/bds/actions/workflows/stress_test.yml)