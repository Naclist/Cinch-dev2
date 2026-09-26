# weighted MI v2 field dictionary

- `mi`: sample-reweighted mutual information of binary locus presence states, in natural-log units (nats), estimated without a pseudocount.
- `nmi`: symmetric normalized MI, `2 * mi / (H(X) + H(Y))`, using the same weighted empirical probabilities.
- `gene1_frequency_0`, `gene1_frequency_1`: thesis equation (2-11) `beta` and `alpha` for `gene1`; weighted absence and presence proportions after sample weights are normalized to sum one.
- `gene2_frequency_0`, `gene2_frequency_1`: corresponding weighted state proportions for `gene2`.
- `frequency_weight`: compatibility scalar fixed to `1.0`. Thesis equation (2-11) supplies a two-component marginal mixing vector, not a second scalar multiplier for MI.
- `frequency_weighted_MI`: equal to `mi`. The alpha/beta locus frequencies are already inside the joint and marginal probabilities defining MI; multiplying MI again would double-weight frequency and would not reproduce the audited formula.
- `sample_weight_scheme`: provenance label for the sample weights.
- `sample_weight_tau`: Hamming-distance neighbourhood cutoff. Similarity >=0.90 is distance <=0.10 in the frozen upstream calculation.
- `sample_weight_sum_raw`: sum of unnormalized inverse-neighbour weights; reported as an effective-size diagnostic.
- `sample_weight_sum_normalized`: sum after thesis equation (2-7) normalization; expected to be 1.

The retained source fields (`gene1`, `gene2`, `r2_adj`, distances, and related columns) keep their definitions from `kp_universe_spec.md`. This is weighted MI v2, not EpiDis.
