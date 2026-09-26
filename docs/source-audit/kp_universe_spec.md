# Kp universe and existing-pipeline freeze

Frozen on: 2026-07-13

Scope: frozen files under `<workspace>/Leca`; no new association method is specified here.

## 1. Profile universe and encoding

- [COMPUTED, HIGH] Primary profile: `Leca/source/profile/0427_1.profile`.
- [COMPUTED, HIGH] Shape: 1,000 genome rows and 41,312 locus columns, plus the first `genome` identifier column.
- [KNOWN, HIGH] The table is an allele profile: positive integers are allele identifiers at a locus. Existing scripts consistently define callable/present as `value > 0` and absent/non-callable as `value <= 0` after filling missing values with zero.
- [UNRESOLVED, MED] The upstream semantic distinction among literal `0`, negative values, missing calls, paralogy, and assembly failure is not documented by the inspected files. Existing downstream code collapses all of them to not-present/not-callable.
- [COMPUTED, HIGH] At prevalence >=0.01, 8,651 loci enter the current wg set; at prevalence >=0.95, 3,632 loci enter the current cg set (`profile_wg_cg_similarity.py` outputs).
- [COMPUTED, HIGH] The existing wgMLST association pair universe contains 3,279 loci, which is narrower than both the full 41,312-locus profile and the 8,651-locus similarity set.

## 2. Representations currently used

- [KNOWN, HIGH] Both representations are used, but for different tasks.
- [KNOWN, HIGH] Presence/absence is `profile > 0`; it is used by `Cinch_v8.py` for wgMLST `r2_adj`, by `spydrpick_locus_wg_nmi.py` for weighted binary MI/NMI, and by the strict12/top20 presence heatmaps.
- [KNOWN, HIGH] Allele identity is used by `profile_wg_cg_similarity.py` to calculate wg/cg allele concordance and by `hiercc_shc_phylo_correct_mi.py` to construct the cgMLST HieCC-like partition.
- [KNOWN, HIGH] The current weighted MI/NMI association score is not an allele-state association score; it is a binary locus-presence association score.

## 3. Existing pair universe

- [KNOWN, HIGH] Source pair table: `Leca/source/profile/0427_1.ldlike_pairs.tsv.gz`, produced by `Leca/code/Cinch_v8.py`.
- [COMPUTED, HIGH] It contains 2,609,249 pair rows plus one header and 3,279 represented loci.
- [KNOWN, HIGH] `Cinch_v8.py` first prevalence-filters binary locus presence and writes a pair only when `n_joint >= min_joint_count`; therefore this is a filtered observed-pair universe, not all combinations among the 41,312 profile loci.
- [UNRESOLVED, MED] The exact command line used for the frozen `0427_1.ldlike_pairs.tsv.gz`, including `min_prev`, `max_prev`, `min_joint_count`, `n_pcs`, tree, positions file, and linear/circular distance mode, is not stored beside the inspected pair file.
- [KNOWN, HIGH] Derived weighted pair table: `Leca/results/gene_gene_wgmlst/spydrpick_locus_nmi/wgmlst_all_weighted_mi_nmi_pairs.tsv.gz`; it preserves the same 2,609,249 rows and appends `weighted_mi` and `weighted_nmi`.

## 4. Frozen pair-field definitions

| Field | Frozen meaning |
|---|---|
| `r2_adj` | [KNOWN, HIGH] Squared Pearson correlation between two binary presence vectors after projecting each centered locus vector away from the first `n_pcs` PCA score subspace in `Cinch_v8.py`. [UNRESOLVED, MED] The frozen run's exact `n_pcs` value is not recoverable from the result table alone. |
| `weighted_mi` | [KNOWN, HIGH] Mutual information in natural-log units for the 2x2 binary presence table, using per-genome `wg_content_weight` and pseudocount 0.5 in `spydrpick_locus_wg_nmi.py`. |
| `weighted_nmi` | [KNOWN, HIGH] `2 * weighted_mi / (H(X) + H(Y))`, using the same weighted, pseudocount-adjusted binary probabilities. |
| `mean_order_dist` | [KNOWN, HIGH] Mean minimum within-contig locus-order separation across genomes where both loci are placed on the same contig/replicon; the source script optionally supports linear or circular order. [UNRESOLVED, MED] The exact distance mode used for the frozen source pair file is not recorded beside it. |
| `mean_bp_dist` | [KNOWN, HIGH] Mean minimum interval-to-interval physical distance in bp across genomes where both loci are placed on the same contig/replicon; the source script optionally supports linear or circular distance. [UNRESOLVED, MED] The exact frozen distance mode is not recorded beside the pair file. |
| `mean_phylo_dist` | [KNOWN, HIGH] In the source `ldlike_pairs` table, mean patristic tree distance among unordered genome pairs jointly carrying both loci. It is computed from the binary presence vectors and the input tree distance matrix. [KNOWN, HIGH] This field was dropped from `wgmlst_all_weighted_mi_nmi_pairs.tsv.gz`; analyses requiring it use the earlier source/filtered tables. |

Additional frozen fields: [KNOWN, HIGH] `n_joint` is the number of genomes carrying both loci; `same_replicon_rate` is the same-contig placement count divided by `n_joint` in `Cinch_v8.py`; `p_adj` is the correlation p-value calculated from adjusted correlation and the configured covariate count.

## 5. Current sample weights

- [KNOWN, HIGH] Generator: `Leca/code/profile_wg_cg_similarity.py`.
- [KNOWN, HIGH] Input: the complete 1,000 x 41,312 allele profile.
- [KNOWN, HIGH] For `wg_content_weight`, loci with prevalence >=0.01 are converted to presence/absence, pairwise genome similarity is Jaccard intersection/union, neighbours are genomes with similarity >=0.90 including self, and each sample receives weight `1 / neighbour_count`.
- [COMPUTED, HIGH] The frozen wg-content effective sample size, defined here as sum of weights, is 101.232972; neighbour counts range from 1 to 400 with median 143.5.
- [KNOWN, HIGH] `spydrpick_locus_wg_nmi.py` and `hiercc_shc_phylo_correct_mi.py` use the `wg_content_weight` column, not the alternative wg-allele or cg-allele weights.
- [UNRESOLVED, HIGH] No independent validation has frozen 0.90 as an optimal biological threshold; it is the current engineering choice and must be sensitivity-tested before treating it as method-invariant.

## 6. Distance-cutoff versions

Two incompatible cutoff families coexist and must remain named explicitly.

### A. Early r2_adj + phylogeny elbow filters

- [COMPUTED, HIGH] Source: `Leca/results/distance_filter/secondary_filter_0427_1/elbow_thresholds.tsv`.
- [COMPUTED, HIGH] Keep `mean_order_dist >= 50.09`.
- [COMPUTED, HIGH] Keep `mean_bp_dist >= 105,585.9`.
- [COMPUTED, HIGH] Keep `mean_phylo_dist >= 0.020894`.
- [KNOWN, HIGH] These were estimated from significant `r2_adj` pairs (`p_adj <= 0.05`) using the 0.90 trend quantile; they define the old r2 PCA edge universes.

### B. Later weighted-MI 99.5% upper-tail segmented elbows

- [COMPUTED, HIGH] Source: `Leca/results/gene_gene_wgmlst/spydrpick_locus_nmi/mi_nmi_distance_comparison/*_mi_aracne_elbow_summary.tsv`.
- [COMPUTED, HIGH] Order-distance cutoff: `mean_order_dist >= 102.85`; weighted-MI Tukey threshold 0.1479019902.
- [COMPUTED, HIGH] BP-distance cutoff: `mean_bp_dist >= 246,474.5`; weighted-MI Tukey threshold 0.1639267674.
- [KNOWN, HIGH] These use the 99.5th-percentile binned MI trend and segmented linear fits. They define the candidate sets passed to the current SHC correction.

## 7. Current SHC/HieCC-like correction

- [KNOWN, HIGH] Script: `Leca/code/hiercc_shc_phylo_correct_mi.py`.
- [KNOWN, HIGH] Similarity input: `Leca/results/profile_similarity/wg_cg_profile_0427_1/wg_cg_similarity_matrices.npz`, key `cg_allele_concordance_callable`.
- [KNOWN, HIGH] The cg scheme size is fixed to 3,632 loci; distance is `3632 * (1 - cg allele concordance)` and clustering uses single linkage.
- [KNOWN, HIGH] Stable blocks require all within-block partition NMI values >=0.90 and at least three levels; the first stable block passing singleton fraction <=0.20 is considered, then the level with best silhouette is selected.
- [COMPUTED, HIGH] Frozen selection: HC426; stable block HC129-HC430; 142 SHC clusters; selected silhouette approximately 0.698.
- [KNOWN, HIGH] Candidate inputs are the MI upper-tail order and bp outlier tables described above.
- [KNOWN, HIGH] Per candidate: SHC size >=5; both loci must vary within an SHC; at least two informative SHCs and at least two positive-direction SHCs are required; 999 within-SHC label permutations are run; BH `q < 0.05` and conditional-MI z-score >0 define phylogenetic passage; ARACNE then labels direct/indirect edges.
- [COMPUTED, HIGH] Order candidates: 3,941; eligible: 108; q-passing: 16; strict direct: 16.
- [COMPUTED, HIGH] BP candidates: 3,793; eligible: 90; q-passing: 28; strict direct: 25; strict indirect: 3.

## 8. Frozen output locations

- [KNOWN, HIGH] strict12: `Leca/results/gene_gene_wgmlst/spydrpick_locus_nmi/hiercc_shc_correction/strict12_visualization/` (12 loci).
- [KNOWN, HIGH] top20 order/BP Phandango-like outputs: `Leca/results/gene_gene_wgmlst/spydrpick_locus_nmi/top20_phandango_presence/` (20 pairs each; 15 order loci and 20 BP loci).
- [KNOWN, HIGH] Current full-MI-matrix PCA attempt: `Leca/results/gene_gene_wgmlst/spydrpick_locus_nmi/mi_pca/` (3,279 loci). [KNOWN, HIGH] This is not the PCA definition accepted for the old r2 comparison.
- [KNOWN, HIGH] Accepted MI-profile PCA on the old filtered edge universes: `Leca/results/metadata_ordination/order_bp_phylo_compare_0427_1_mi_profile_pca/` (3,274 loci for each panel).
- [KNOWN, HIGH] Original r2-profile PCA: `Leca/results/metadata_ordination/order_bp_phylo_compare_0427_1/` (3,274 loci; order 756,541 edges; BP 770,117 edges).

## 9. Valid EpiDis comparison baselines

- [KNOWN, HIGH] Frozen Kp genome universe: the same 1,000 profile rows in the same order, after explicit sample-ID alignment.
- [KNOWN, HIGH] Frozen candidate/pair universe: the 2,609,249 rows in `0427_1.ldlike_pairs.tsv.gz`, if EpiDis is evaluated only where current MI/r2 results exist.
- [KNOWN, HIGH] Binary weighted MI and NMI in `wgmlst_all_weighted_mi_nmi_pairs.tsv.gz` are valid score-level baselines if EpiDis uses the identical binary coding, genome universe, pair universe, missingness rule, and sample weights.
- [KNOWN, HIGH] Early r2 PCA edge sets are valid structural baselines for comparing profile geometry, but their cutoff-generating metric is `r2_adj`, not MI or EpiDis.
- [KNOWN, HIGH] MI 99.5%-elbow candidate sets and SHC-filtered direct edges are valid pipeline-level baselines for candidate yield, edge overlap, distance distribution, SHC eligibility, and network topology.
- [KNOWN, HIGH] strict12 and top20 outputs are valid named-set/visualization baselines, not unbiased performance benchmarks because they were selected from the existing MI workflow.

## 10. Results that must be recomputed for a fair EpiDis comparison

- [KNOWN, HIGH] EpiDis scores for every pair must be computed from the frozen underlying profile; they cannot be inferred from MI, NMI, or r2 columns.
- [KNOWN, HIGH] Any EpiDis-specific score threshold, null distribution, IQR cutoff, or distance elbow must be re-estimated on EpiDis values; MI/r2 numeric cutoffs cannot be transferred.
- [KNOWN, HIGH] Direct/indirect ARACNE classification must be rerun using EpiDis as edge weight if ARACNE is part of the comparison.
- [KNOWN, HIGH] SHC conditional score, within-SHC permutation null, p/q values, and strict-edge sets must be rerun for EpiDis; current conditional MI results are not EpiDis results.
- [KNOWN, HIGH] EpiDis PCA, clustering, network, top20 and Phandango-like panels must be rebuilt from EpiDis-selected edges.
- [KNOWN, HIGH] If the frozen EpiDis specification uses allele states rather than presence/absence, all current binary MI/r2 comparisons become representation-mismatched; an allele-state MI control must then also be recomputed.
- [KNOWN, HIGH] If EpiDis sample weighting differs from `wg_content_weight`, both weighted MI and EpiDis must be recomputed under a shared weight scheme plus a clearly labelled native-weight comparison.
- [UNRESOLVED, HIGH] Whether EpiDis should use all 41,312 loci, the 8,651 wg set, the 3,279 association loci, or another prevalence-defined set must be fixed by the separate method/source audit before implementation.

## Freeze decision

- [KNOWN, HIGH] No downstream agent may call the current weighted MI an allele association: it is binary presence association.
- [KNOWN, HIGH] No downstream agent may reuse an r2 or MI cutoff as an EpiDis threshold without recomputing the corresponding EpiDis empirical distribution.
- [KNOWN, HIGH] No downstream agent may mix the early r2 elbow cutoffs (`50.09`, `105585.9`, `0.020894`) with the later MI 99.5% cutoffs (`102.85`, `246474.5`) without naming which candidate universe is intended.
