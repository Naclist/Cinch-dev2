# Outputs and interpretation

## Score outputs

- `MI`: weighted mutual information in nats unless explicitly labelled otherwise.
- `NMI`: `2*MI/(H(X)+H(Y))`.
- `EpiDis`: square root of weighted MI in bits under the frozen construction.
- physical distances: mean order and base-pair separation where available.
- SHC fields: informative clusters, conditional MI, permutation null mean, z score, p value, and BH q value.
- ARACNE field: whether an edge survived post-test triangle pruning.
- Diff-GWES fields: discovery and held-out validation effects for each background.

## Claims that outputs do not support by themselves

No single score establishes causality, direct biochemical interaction, positive selection, functional cooperation, ecological competition, or clinical relevance. Distance filtering reduces one confounder; sample weighting and SHC permutation address different aspects of population structure; ARACNE addresses network redundancy.

## Frozen aggregate result tables

Small aggregate tables under `docs/source-results` are retained for audit and figure regeneration. Raw profiles, complete pair matrices, genomes, toy data, and smoke-test data are intentionally absent.

