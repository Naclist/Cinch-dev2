# Cinch algorithm results documentary

Frozen Kp analysis date: 2026-07-13.

## Pair universe

| Item | Frozen value |
|---|---:|
| Genomes | 1,000 |
| Full profile loci | 41,312 |
| Pair-universe loci | 3,279 |
| Valid pair rows | 2,609,249 |

## Weighted MI v2 validation

All 2,609,249 pair keys are unique and all required fields are populated. MI ranges from approximately `3.70e-14` to `0.693138` nats; NMI ranges from approximately `6.29e-14` to `1.0`. State-frequency sums differ from one by at most `2.22e-16`. The stored validation flag is true.

## Pairwise EpiDis distribution

| Statistic | Value |
|---|---:|
| Valid EpiDis pairs | 2,609,249 |
| Minimum | 2.3095417906209117e-07 |
| Median | 0.1178323812623859 |
| Mean | 0.15238904206705234 |
| Maximum | 0.9999934522320495 |
| Maximum directional absolute difference | 1.8206777849789626e-10 |
| Q3 + 3 IQR threshold | 0.684741012905513 |
| IQR outliers | 15,888 |
| Outlier fraction | 0.006089108398623512 |

## Metric comparison

| Comparison | Pearson | Spearman | Main implication |
|---|---:|---:|---|
| EpiDis vs MI v2 | 0.912491 | ~1.0 | Same rank, nonlinear square-root scale |
| EpiDis squared vs MI v2 bits | 1.0 | ~1.0 | Algebraically equivalent in frozen implementation |
| EpiDis vs old weighted MI | 0.585771 | 0.360846 | Old weights and pseudocount materially change rank |
| EpiDis vs old weighted NMI | 0.606028 | 0.340807 | Entropy normalization and old construction change rank |

## Distance result

The highest EpiDis pairs are heavily enriched for physically adjacent loci. Examples include `ref_gene_692-ref_gene_693` with order distance 1 and BP distance 14, and `ref_gene_692-ref_gene_694` with order distance 2 and BP distance 477. Raw upper-tail membership is therefore not evidence of remote epistasis.

## SHC-conditioned Diff-GWES validation

| Check | Result |
|---|---:|
| True-data eligible hypotheses | 351 |
| True-data significant hypotheses | 0 |
| True-data minimum q | 0.5383435582822086 |
| Synthetic truth triples | 20 |
| Synthetic truth triples passing | 20 |
| Predeclared negative controls passing | 0 of 20 |
| Non-truth synthetic hypotheses passing | 0 |

The toy benchmark supports implementation sensitivity and specificity under its own synthetic design. It does not validate biological truth in the observed Kp data. The true-data result does not support a significant background-dependent interaction under the frozen eligibility, split, permutation, and correction rules.

## Reproducibility gaps

- The exact command that created `0427_1.ldlike_pairs.tsv.gz` is not stored beside the pair file.
- The profile contract collapses missing/non-callable values with absence.
- Only 952 of 3,279 loci have positions in the single-reference flattened coordinate map used by frozen Diff-GWES; local-block preservation is incomplete for unplaced accessory loci.
- The public audited EpiDive repository did not contain the thesis pairwise EpiDis implementation, so the thesis formula is the score authority.

