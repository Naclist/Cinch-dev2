# Cinch frozen methodology

## Interpretation boundary

Cinch detects phylogeny-aware statistical association candidates in bacterial wg/cgMLST data. A retained edge or community is not proof of molecular epistasis, fitness interaction, or causality.

## Frozen representation

- Kp profile: 1,000 genomes by 41,312 loci, plus genome identifier.
- Frozen pair universe: 3,279 loci and 2,609,249 observed pairs.
- Presence: allele-profile value greater than zero.
- Absence/non-callable: value less than or equal to zero after missing values are filled with zero.
- Known limitation: the input contract does not distinguish biological absence from missingness, assembly failure, or other non-positive calls.

## Layered workflow

### 1. Mapping and pair-universe construction

`Cinch_v8.py` deduplicates reference CDS sequences, maps genome sequences, constructs allele profiles, and computes pairwise fields. For wgMLST it uses binary presence, raw correlation, PCA-residualized correlation, tree distance, and physical distance. For cgMLST it computes allele-state NMI and a K-means-stratified permutation z-score.

This stage is an upstream candidate universe. PCA residualization is not the final SHC-based phylogenetic test.

### 2. Existing weighted MI comparator

The original comparator computes binary weighted MI/NMI from a 2 by 2 table with inverse-neighbour sample weights and a 0.5 pseudocount, then applies distance summaries and ARACNE-style pruning. Keep this result as a comparator because the pseudocount and weight definition alter pair ranking.

### 3. Weighted MI v2

Weighted MI v2 uses thesis-style Hamming-neighbour sample weights at tau 0.10, normalized to sum one, and no pseudocount. It reports MI in nats and symmetric NMI. The thesis alpha/beta state frequencies already enter the joint and marginal probabilities; no extra scalar locus-frequency multiplier is applied.

### 4. Pairwise EpiDis

The frozen implementation directly evaluates thesis equations 2-8 through 2-14 using base-2 weighted Jensen-Shannon divergence and then takes the square root. Both conditioning directions are retained and numerically audited.

Under the frozen normalized weighted empirical probabilities:

```text
EpiDis^2 = MI_bits = MI_nats / ln(2)
```

Therefore pairwise EpiDis and weighted MI v2 have the same ranking. EpiDis changes the scale and distribution-based threshold, not the underlying pair ordering. This equivalence does not erase the thesis high-order framework; it limits the novelty of the pairwise core.

### 5. Physical-distance layer

Inspect order distance and base-pair distance before interpreting the high-score tail. The frozen Kp EpiDis tail is dominated by short-distance pairs. Do not treat the global IQR threshold as a linkage correction.

### 6. SHC phylogenetic layer

Construct stable hierarchical clusters from the cgMLST-like profile and test association within informative SHCs. Permute within SHC, compute empirical p-values, and apply Benjamini-Hochberg correction. Require a positive effect direction and adequate informative clusters. This is stronger than inverse-neighbour weighting alone because it conditions the null on population structure.

### 7. Network layer

Apply ARACNE-style triangle pruning only after statistical and distance filtering. It removes a uniquely weakest triangle edge as a possible indirect association. It is a topology heuristic, not a significance test.

### 8. SHC-conditioned Diff-GWES

For seed pair A-B and background locus C, compare SHC-conditioned information between C=1 and C=0:

```text
I_c(A;B | SHC) = sum_h P(h | c, informative) I(A;B | h,c)
E_c = sqrt(I_c / ln(2))
Delta(A,B;C) = E_1 - E_0
```

Use the frozen 50:50 within-SHC discovery/validation split, background prevalence 0.05-0.95, effect screen 0.15, 999 within-SHC permutations, and BH correction. Interpret passing triples as phylogeny-conditioned statistical interactions.

## Source authority

Use the thesis for EpiDis formulas and theoretical claims, audited source code for implemented encoding utilities, and Cinch frozen scripts/provenance for Kp-specific decisions. Record disagreements; do not merge them silently.

