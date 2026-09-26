# Performance architecture

## Retained fast paths

The local implementation retains three useful computational patterns:

1. binary weighted joint masses are computed by dense matrix multiplication rather than one Python task per pair;
2. multiallelic cgMLST NMI uses dense integer codes and Numba-parallel triangular pair enumeration;
3. expensive statistical layers operate on filtered candidate sets rather than rerunning every downstream procedure on every pair.

## Known bottlenecks

- the full dense binary joint matrix requires quadratic memory in the number of loci;
- multiallelic contingency construction still scales approximately with samples times pairs;
- SHC-preserving permutation is expensive and currently pair-oriented;
- historical scripts perform repeated dataframe serialization;
- the original WGS mapping frontend depends on external wrappers not included here.

## Exact-equivalent optimization boundary

Permitted engineering changes include dense integer coding, precomputation, blockwise triangular enumeration, compiled contingency kernels, shared read-only matrices, cached legal permutations, batched multiprocessing, and chunked output. Candidate screening, fewer permutations, adaptive stopping, restricted high-order backgrounds, and hierarchical hypothesis testing change the statistical program and require a separately validated method version.

