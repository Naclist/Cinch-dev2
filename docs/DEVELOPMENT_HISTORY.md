# Development history and repository design

Cinch dev2 combines two bodies of work without claiming they were previously equivalent.

## Retained from the local workflow

- the complete historical and extension stage inventory;
- fast matrix and Numba information kernels;
- weighted MI v2 and pairwise EpiDis;
- physical-distance diagnostics;
- stable hierarchical cluster permutation, BH correction, and ARACNE pruning;
- SHC-conditioned Diff-GWES;
- allele-diversity versus same-allele phylogenetic-distance analysis;
- aggregate validation results and formula/code audit documents.

## Release conventions adopted from the public repository style

- a `pyproject.toml` package;
- one installable CLI;
- machine-readable freeze and workflow contracts;
- CI tests and integrity verification;
- a static documentation site;
- deterministic figure generation;
- clear separation of source, documentation, result summaries, and excluded data.

## Explicitly not adopted

- a pure-Python five-seed mapper as a production allele caller;
- a `--threads` option that only records provenance;
- per-pair Python `np.unique` and contingency allocation as the main all-pairs engine;
- documentation of a non-callable presence state that the mapper cannot emit;
- silent replacement of sample-weighted inference with unweighted counts;
- recurrence/Neff filtering as a substitute for SHC permutation p-values.

## Data exclusion

The repository does not contain virtual environments, caches, smoke outputs, demonstration genomes, synthetic toy profiles, complete research pair tables, or raw chat exports. Their exclusion is intentional and documented, not evidence that the corresponding methods were removed.

