# EpiDive / EpiDis implementation-grade method audit

## Audit identity and evidence boundary

- [KNOWN, HIGH] Thesis source: Xu Wenxuan, *Epistasis Divergence (EpiDis)-Based Genome-wide Pairwise and High-order Epistasis Analysis in Bacteria*, PDF dated 2026, audited at printed pages 25-52, especially equations (2-3)-(2-16) and (3-1)-(3-8). Printed page numbers, not PDF object indices, are used below.
- [KNOWN, HIGH] Repository source: `Riffleriver/EpiDive_consis`, default branch `main`, commit `2101406bb26136e9eaa7b84d535b3a0bd5dc508b` (2025-04-21). Audited files: `README.md`, `R_EpiDis_matrix10.py`, `R_consis_snp.py`, and `R_classify_snps_near_faiss.py`.
- [KNOWN, HIGH] Current-project source: `Leca/code/profile_wg_cg_similarity.py`, `Leca/code/spydrpick_locus_wg_nmi.py`, and `Leca/code/hiercc_shc_phylo_correct_mi.py`.
- [KNOWN, HIGH] The GitHub repository does **not** implement the EpiDis formulas described in the thesis. Its file named `R_EpiDis_matrix10.py` only converts each SNP row to a major-state indicator. This gap is an audit result, not an implementation detail to fill by assumption.

## 1. Problems solved by EpiDis and EpiDive

- [KNOWN, HIGH] **EpiDis** is the pairwise statistic. It measures how much the distribution of locus `j` changes between the two states of conditioning locus `i`, after sample reweighting and conditioning-state proportion weighting (thesis pp. 32-37, eqs. 2-3 and 2-8 to 2-14).
- [KNOWN, HIGH] **EpiDive** is the analysis framework around EpiDis: binary encoding, pairwise genome-wide scanning, empirical outlier selection, distance plotting, sparse-network construction, CopSCC/Tarjan decomposition, and a differential-background extension for higher-order interactions (thesis pp. 40-52, eqs. 3-1 to 3-8).
- [KNOWN, HIGH] The thesis defines high-order interaction as a change in a pair's EpiDis across background subsets, not as exhaustive enumeration of triples or larger tuples (thesis pp. 47-52, sections 3.3-3.4).
- [UNRESOLVED, HIGH] The public repository cannot reproduce either the pairwise EpiDis scan or the high-order EpiDive workflow because it contains no divergence, weighting, IQR, network, CopSCC, or differential-background implementation.

## 2. Input data and state encoding

### Thesis definition

- [KNOWN, HIGH] The pairwise input is a loci-by-samples genotype matrix. SNP states are reduced to binary: the most frequent state at a locus is `1` (major state), and all other observed states are merged into `0` (minor state) (thesis p. 34, eqs. 2-4 and 2-5).
- [KNOWN, HIGH] For pangenome presence/absence data, presence is `1` and absence is `0` (thesis pp. 34-35, section 2.3.1).
- [KNOWN, HIGH] Therefore `1` has two source-specific meanings: major allele for SNP input, presence for gene-content input. These meanings are not interchangeable.
- [KNOWN, HIGH] The thesis states that gaps/missing states in SNP alignments are filled using a population-frequency-weighted imputation procedure (thesis p. 34).
- [UNRESOLVED, HIGH] The thesis does not provide an explicit formula, sampling rule, deterministic tie rule, or random seed for that imputation.

### Repository behavior

- [KNOWN, HIGH] `R_EpiDis_matrix10.py::_transform_row` computes `row.value_counts().idxmax()` and encodes equality to that value as `1`; every other value becomes `0` (lines 78-81). It does not explicitly exclude missing values or gaps before choosing the major state.
- [KNOWN, HIGH] Ties are resolved by `pandas.Series.value_counts().idxmax()` according to returned ordering; the repository does not document this as a biological rule.
- [KNOWN, HIGH] `R_consis_snp.py` groups exact binary row patterns by MD5, while `R_classify_snps_near_faiss.py` groups or queries rows by Hamming distance. Neither computes an association statistic.

### Kp implementation decision

- [KNOWN, HIGH] For Kp gene-content presence/absence analysis, use `profile[locus] > 0` as `1` and `0`, missing, or non-positive calls as `0` only if zero is already defined as absence/non-call in the profile contract.
- [INFERRED, HIGH] To avoid conflating absence with missingness, a production implementation should accept an explicit callable mask and report callable counts per locus pair. This is a required Kp adaptation because the current profile encoding collapses both to zero.
- [KNOWN, HIGH] Do not recode Kp presence/absence loci by majority state. Doing so would turn common presence into `1` but rare presence into `0`, changing the biological meaning across loci.

## 3. Sample-level weighting

### Thesis formula

[KNOWN, HIGH] On thesis p. 35, equation (2-6), for sample `k`:

```text
w_k^raw = 1 / count{ l : (1/n) * sum_{r=1}^n Hamming(G_{k,r}, G_{l,r}) < tau }
```

[KNOWN, HIGH] Here `n` is the number of loci, `G_{k,r}` is sample `k` at locus `r`, and the expression inside the set is the genome-wide mean Hamming difference proportion between samples `k` and `l`.

[KNOWN, HIGH] Equation (2-7) normalizes the weights:

```text
w_k = w_k^raw / sum_{s=1}^m w_s^raw
sum_{k=1}^m w_k = 1
```

- [KNOWN, HIGH] The thesis states `tau in [0.1, 0.25]` and calls it a similarity threshold, although the inequality is written on a **difference** proportion (`distance < tau`) (thesis p. 35).
- [INFERRED, HIGH] The naming is internally confusing: the mathematical object is a distance cutoff even though the prose calls it similarity.
- [KNOWN, HIGH] The neighborhood count includes the focal sample because its self-distance is zero and satisfies any positive `tau`.

### Repository correspondence

- [UNRESOLVED, HIGH] No repository file computes equation (2-6) or normalizes equation (2-7).

### Current Kp weighted-MI correspondence

- [KNOWN, HIGH] The current project uses `w_k = 1 / N_k`, where `N_k` is the number of samples with a selected similarity matrix value `>= 0.90` (`profile_wg_cg_similarity.py`, `weights_from_similarity`, lines 82-84 and 151-160).
- [KNOWN, HIGH] The current MI scan uses `wg_content_weight` derived from wg gene-content Jaccard similarity, not the thesis's binary-genotype mean Hamming distance (`spydrpick_locus_wg_nmi.py`, lines 221-225).
- [INFERRED, HIGH] These schemes are conceptually related inverse-neighborhood weights but are not formula-equivalent. A faithful EpiDis implementation must expose `distance metric`, `tau`, `self included`, and `normalization` as provenance fields.

## 4. Conditional distributions and locus-level weighting

[KNOWN, HIGH] For a conditioning locus `i`, thesis p. 35 equation (2-8) defines weighted masses:

```text
d_p = sum_{k: x_{i,k}=1} w_k
d_q = sum_{k: x_{i,k}=0} w_k
```

[KNOWN, HIGH] Thesis p. 36 equation (2-9) defines the conditional probability that target locus `j` is in state `1`:

```text
p_1 = [sum_{k: x_{i,k}=1} x_{j,k} w_k] / d_p
q_1 = [sum_{k: x_{i,k}=0} x_{j,k} w_k] / d_q
p_0 = 1 - p_1
q_0 = 1 - q_1
```

[KNOWN, HIGH] Equation (2-10) forms binary distributions and adds a numerical constant component-wise:

```text
p = (p_1, p_0) + epsilon
q = (q_1, q_0) + epsilon
epsilon = 1e-27
```

[KNOWN, HIGH] Equation (2-11) defines locus-level/background-proportion weights:

```text
alpha = size_ratio_p = d_p / (d_p + d_q)
beta  = size_ratio_q = 1 - alpha
```

- [KNOWN, HIGH] This is not a separate per-locus multiplicative penalty such as `4 f(1-f)`. It is the empirical weighted proportion of the two conditioning backgrounds.
- [KNOWN, HIGH] Because normalized sample weights sum to one and all samples are assigned to one background, `d_p + d_q = 1` when there are no missing calls.
- [INFERRED, HIGH] Equations (2-9) to (2-13) are written directionally (`i` conditions `j`), but with exact normalized empirical probabilities their weighted JSD equals `I(X_i;X_j)` and is therefore symmetric. The repository provides no executable test of this identity.

## 5. Mixture distribution and EpiDis formula

[KNOWN, HIGH] Thesis p. 36 equation (2-12) defines the mixture distribution:

```text
u = alpha * p + beta * q
```

[KNOWN, HIGH] Equation (2-13) is a base-2, proportion-weighted Jensen-Shannon divergence:

```text
e_JS = alpha * sum_{s in {1,0}} p_s log2(p_s / u_s)
     + beta  * sum_{s in {1,0}} q_s log2(q_s / u_s)
```

[KNOWN, HIGH] Equation (2-14) defines EpiDis:

```text
d_EpiDis(x_i, x_j) = sqrt(e_JS)
```

- [KNOWN, HIGH] With valid normalized distributions and base-2 logarithms, weighted JSD is in `[0,1]`, hence EpiDis is in `[0,1]` (thesis pp. 33 and 36).
- [KNOWN, HIGH] EpiDis is zero when the two conditional distributions are identical; statistical independence implies identical conditionals when both conditioning states have positive mass.
- [KNOWN, HIGH] The square root turns JSD into Jensen-Shannon distance, a metric on probability distributions.
- [UNRESOLVED, HIGH] The thesis adds `epsilon` to both entries of `p` and `q` but does not state that the vectors are renormalized afterward. Without renormalization, they sum to `1 + 2 epsilon`, so the strict probability-distribution proof is only numerically approximate.
- [UNRESOLVED, HIGH] If `d_p=0` or `d_q=0`, equation (2-9) divides by zero. The thesis does not state whether monomorphic conditioning loci are removed, skipped, or regularized.
- [UNRESOLVED, HIGH] No repository code implements equations (2-8) to (2-14).

## 6. Symmetry, range, and boundary cases

- [KNOWN, HIGH] The thesis explicitly claims `D_EpiDis(X_i,X_j)=D_EpiDis(X_j,X_i)` (p. 33).
- [INFERRED, HIGH] With exact normalized probabilities, equation (2-13) is mutual information expressed as weighted JSD, so the symmetry claim is mathematically consistent even though the formula is written conditionally.
- [UNRESOLVED, HIGH] Literal component-wise addition of `epsilon` without renormalization can introduce a tiny directional discrepancy. A faithful audit implementation must compute both directions on contingency-table test cases and report the discrepancy before collapsing to one symmetric value.
- [KNOWN, HIGH] The claimed theoretical range is `[0,1]` (pp. 33, 36). The extreme value `1` requires maximally separated conditionals with compatible mixture weights under base-2 JSD.
- [KNOWN, HIGH] Constant target locus `j` gives identical conditionals and EpiDis zero when both conditioning groups exist.
- [UNRESOLVED, HIGH] Constant conditioning locus `i` is undefined under the displayed formula because one denominator is zero.
- [UNRESOLVED, HIGH] Missing-value handling after subgrouping and reweighting is not fully specified.

## 7. Missing values, zero frequencies, and pseudocounts

- [KNOWN, HIGH] Thesis SNP input: missing/gap calls are said to use population-frequency-weighted imputation (p. 34), but no algorithm is supplied.
- [KNOWN, HIGH] Thesis numerical zero handling: add `epsilon=10^-27` to each binary conditional-distribution component (p. 36, eq. 2-10).
- [KNOWN, HIGH] The thesis does not use a Jeffreys count pseudocount in the displayed EpiDis formulas.
- [KNOWN, HIGH] Current weighted MI instead adds `0.5` to every cell of the weighted 2x2 table before normalization (`spydrpick_locus_wg_nmi.py`, lines 39-48 and default at line 198).
- [INFERRED, HIGH] These two stabilizers are not equivalent: `epsilon` prevents logarithms of zero after conditional estimation, whereas `0.5` changes all joint-cell estimates and therefore MI/NMI values materially at low effective counts.

## 8. Original outlier and IQR thresholds

[KNOWN, HIGH] Pairwise EpiDis threshold, thesis p. 37 equation (2-16):

```text
IQR = Q3 - Q1
T_high = Q3 + k * 1.5 * IQR
k in {1,2,3}; default k = 2
significant if d_EpiDis(i,j) > T_high
```

- [KNOWN, HIGH] Thus the default multiplier on IQR is `3.0`, not `2.0`.
- [KNOWN, HIGH] The threshold is calculated on the full set of pairwise EpiDis values, not per-locus maxima, according to equation (2-15) and the accompanying text.
- [KNOWN, HIGH] High-order differential threshold, thesis p. 51 equations (3-6) to (3-8), applies two-sided fences to the distribution of `diff_C`: `Q3 + k*1.5*IQR` or `Q1 - k*1.5*IQR`, with `k in {1,2,3,4,5}` and default `k=3`.
- [UNRESOLVED, HIGH] The repository implements none of these thresholds.
- [KNOWN, HIGH] The current weighted-MI workflow instead applies `Q3+1.5*IQR` and `Q3+3*IQR` to **per-locus maximum non-LD NMI**, after an order-distance elbow filter (`spydrpick_locus_wg_nmi.py`, lines 239-261). This is a different null population and must remain separately named.

## 9. EpiDive network construction and high-order extension

- [KNOWN, HIGH] Pairwise network: retain locus pairs above the EpiDis threshold as weighted edges, store a sparse adjacency representation, and decompose with CopSCC using Tarjan SCC (thesis pp. 43-45, eq. 3-1).
- [INFERRED, HIGH] A pairwise EpiDis statistic is naturally undirected if symmetric, but SCC is meaningful only for directed graphs. The thesis says the graph is directed without specifying edge orientation. This is an internal methodological gap.
- [UNRESOLVED, HIGH] The repository contains no sparse EpiDis network, Tarjan algorithm, CopSCC, edge orientation rule, module scoring, or network export.
- [KNOWN, HIGH] High-order EpiDive/Diff-GWES partitions samples by a background variable `C`, recomputes EpiDis within each subgroup, and compares subgroup values (thesis pp. 48-51).
- [KNOWN, HIGH] Equations (3-2) and (3-3): `d1(i,j)=EpiDis_G1(i,j)`, `d2(i,j)=EpiDis_G2(i,j)`, and `diff_C(i,j)=d1(i,j)-d2(i,j)`.
- [KNOWN, HIGH] Equations (3-4) and (3-5) also define subgroup-to-global deviations `d1-dG` and `d2-dG`, each in `[-1,1]` if EpiDis is in `[0,1]`.
- [KNOWN, HIGH] Multi-category backgrounds are to be one-hot expanded into mutually exclusive subsets (thesis p. 50).
- [UNRESOLVED, HIGH] The thesis does not specify minimum subgroup size, handling of monomorphic loci within subsets, multiple-testing control, or whether sample-neighborhood weights are recomputed within each subgroup versus inherited and renormalized. The prose says the weighting strategy continues in each subgroup but is not algorithmically precise.

## 10. Thesis-to-GitHub correspondence

- [KNOWN, HIGH] Binary major-state encoding maps to `R_EpiDis_matrix10.py::_transform_row`.
- [KNOWN, HIGH] Exact-pattern grouping maps to `R_consis_snp.py::classify_snps_by_pattern`; this step is not an EpiDis formula in the thesis.
- [KNOWN, HIGH] Approximate-pattern grouping maps to `R_classify_snps_near_faiss.py::build_index`, `query_single`, `IndexBinaryFlat.range_search`, and SciPy `connected_components`; this is Hamming-neighbor clustering, not EpiDis or CopSCC.
- [UNRESOLVED, HIGH] Equations (2-6) through (2-16) and (3-1) through (3-8) have no GitHub implementation at audited commit `2101406...`.

## 11. Thesis steps absent from the repository

- [KNOWN, HIGH] Sample-level inverse-neighborhood weighting and normalization.
- [KNOWN, HIGH] Weighted conditional distributions and locus/background proportion weights.
- [KNOWN, HIGH] Mixture distribution, weighted JSD, square-root EpiDis.
- [KNOWN, HIGH] `epsilon=1e-27` zero handling.
- [KNOWN, HIGH] Pairwise IQR outlier selection and GWES distance plotting.
- [KNOWN, HIGH] Sparse EpiDis network, edge weighting, CopSCC/Tarjan decomposition.
- [KNOWN, HIGH] Differential-background EpiDis, two-sided IQR threshold, and high-order modules.
- [KNOWN, HIGH] GPU/block-matrix EpiDis implementation claimed in the thesis.

## 12. Repository steps not clearly specified in the thesis

- [KNOWN, HIGH] Exact binary-pattern MD5 grouping and group ranking (`R_consis_snp.py`).
- [KNOWN, HIGH] Faiss binary Hamming range search, zero-padding to a multiple of eight bits, undirected adjacency construction, and connected-component clustering (`R_classify_snps_near_faiss.py`).
- [KNOWN, HIGH] User-selected `allowed_diff` mismatch threshold for pattern neighborhoods.
- [KNOWN, HIGH] Chunked joblib major-state transformation with all available CPU cores (`R_EpiDis_matrix10.py`).
- [INFERRED, HIGH] These repository functions look like preprocessing/search utilities for consistency patterns, but there is no code evidence that they produce the thesis's EpiDis statistic.

## 13. Mathematical relationship to MI, weighted MI, NMI, JSD, and SpydrPick

- [KNOWN, HIGH] For two discrete variables, mutual information can be written as a generalized Jensen-Shannon divergence between conditional distributions weighted by the conditioning-variable marginal: `I(X;Y)=sum_x P(x) KL(P(Y|x)||P(Y))`, where `P(Y)=sum_x P(x)P(Y|x)`.
- [INFERRED, HIGH] The thesis equation (2-13), before the square root, is therefore algebraically the weighted mutual information between two binary variables if `p`, `q`, and mixture `u` are valid normalized weighted empirical distributions.
- [INFERRED, HIGH] Under that reading, `EpiDis = sqrt(MI_bits)` for the directional construction, not a fundamentally different dependence functional. Its reported `[0,1]` range follows because binary MI in bits is at most one.
- [KNOWN, HIGH] Weighted MI changes empirical probabilities by sample weights; EpiDis uses the same kind of weighted conditional probabilities and additionally takes a square root after expressing MI as weighted JSD.
- [KNOWN, HIGH] NMI further divides MI by a function of marginal entropies. The current project uses `2 MI / (H(X)+H(Y))` (`spydrpick_locus_wg_nmi.py`, lines 53-57). EpiDis has no entropy normalization in the thesis formula.
- [KNOWN, HIGH] Jensen-Shannon divergence is exactly the form in equation (2-13); square-root JSD is Jensen-Shannon distance.
- [KNOWN, HIGH] SpydrPick uses sample reweighting, MI, an empirical high-MI cutoff, and ARACNE-style pruning. The current project intentionally follows that pattern but uses locus-level presence/absence and its own distance/elbow and SHC filters.
- [INFERRED, HIGH] The strongest mathematical distinction between audited EpiDis and weighted MI is presentation/scale (`sqrt`, base 2, conditional/JSD form) plus its stated high-order subgroup-difference framework. The displayed pairwise core does not establish a distinct information quantity.
- [UNRESOLVED, HIGH] The exact symmetry of the **implemented** statistic still requires testing because the thesis's literal epsilon step does not state renormalization and the repository cannot settle the implementation convention.

## 14. Required adaptations for Kp presence/absence data

1. [KNOWN, HIGH] Preserve semantic encoding: presence `1`, absence `0`; do not majority-recode loci.
2. [INFERRED, HIGH] Separate absence from missing/non-call with a callable mask; otherwise conditional masses and EpiDis are biased for incompletely called loci.
3. [KNOWN, HIGH] Filter loci with no variation and pair directions with `d_p=0` or `d_q=0`; record them as undefined rather than zero.
4. [KNOWN, HIGH] Implement thesis sample weights exactly as a named option using binary-genotype mean Hamming distance and `tau`; keep existing wg-content Jaccard `>=0.90` weights as a separate comparator.
5. [KNOWN, HIGH] Normalize sample weights to sum one before equations (2-8) to (2-14).
6. [KNOWN, HIGH] Implement `epsilon=1e-27` exactly for thesis replication, while also reporting whether conditional vectors are renormalized. Run both interpretations if necessary because the thesis is silent.
7. [KNOWN, HIGH] Compute both `i -> j` and `j -> i` EpiDis on a diagnostic subset. Do not force symmetry until equality is empirically and algebraically established.
8. [KNOWN, HIGH] Apply the thesis pairwise threshold to all valid EpiDis pairs: `Q3 + 3*IQR` for default `k=2`. Keep the current per-locus-maximum MI threshold in a separately labeled comparison.
9. [KNOWN, HIGH] Keep physical/order-distance filtering and SHC phylogenetic permutation outside the faithful EpiDis score. Report raw EpiDis, distance-filtered EpiDis, and SHC-filtered EpiDis as separate stages.
10. [KNOWN, HIGH] For high-order analysis, define the background variable explicitly, recompute/renormalize subgroup weights, require minimum subgroup size and polymorphism, and document these additions as Kp-specific because the thesis does not specify them.
11. [INFERRED, HIGH] Do not use SCC on a symmetric undirected EpiDis graph without an orientation rule. For Kp, use connected components/community detection for symmetric edges, or define and justify directional edges from `i -> j` scores.
12. [KNOWN, HIGH] Export provenance: profile checksum, locus filter, callable rule, distance metric, `tau`, weight normalization, epsilon rule, log base, directionality rule, IQR population, IQR multiplier, distance cutoff, SHC definition, permutation count, and random seed.

## Implementation stop conditions

- [KNOWN, HIGH] A claim of “faithful repository reproduction” is forbidden for the audited commit because the core algorithm is absent.
- [KNOWN, HIGH] A claim of symmetric EpiDis is not implementation-ready until the directional formula is reconciled with the symmetry claim.
- [KNOWN, HIGH] A claim that EpiDis corrects phylogeny is too strong for Kp: inverse-neighborhood weighting reduces clonal redundancy but does not condition on a phylogenetic covariance model or perform a phylogenetic null test.
- [KNOWN, HIGH] Statistical association, EpiDis outlier status, or network membership is not evidence of biological epistasis without additional validation.
