# Scientific narrative

## Question

Cinch asks whether two bacterial loci carry more statistical dependence than expected after accounting for physical linkage, sample redundancy, and population structure. The result is a candidate dependency, not a direct measurement of molecular fitness interaction.

## Representation

The frozen local implementation starts from a genome-by-locus allele profile. The presence layer uses `value > 0` as present and collapses every other value into absence/non-callable. Multiallelic cgMLST analysis retains categorical allele identifiers. These are separate statistical representations and must not be silently interchanged.

## Pairwise information

For normalized sample weights `w_s`, binary states `X` and `Y`, and weighted empirical probabilities,

```text
MI_w(X;Y) = sum_xy p_w(x,y) log[p_w(x,y)/(p_w(x)p_w(y))]
NMI_w     = 2 MI_w / (H_w(X) + H_w(Y)).
```

Weighted MI v2 applies no pseudocount. Pairwise EpiDis evaluates the weighted, base-2 Jensen–Shannon construction in both conditioning directions. Under the frozen valid empirical distribution,

```text
EpiDis(X,Y)^2 = MI_bits(X;Y) = MI_nats(X;Y)/ln(2).
```

EpiDis therefore changes scale and thresholding but does not introduce a new pair ranking.

## Confounding layers

Physical order and base-pair distances are inspected before upper-tail interpretation. Inverse-neighbour sample weights reduce repeated-lineage dominance but do not establish independence from phylogeny. Cinch therefore constructs stable hierarchical clusters, restricts testing to informative clusters, permutes within cluster, and applies Benjamini–Hochberg correction.

ARACNE-style triangle pruning occurs only after statistical and distance filtering. It removes a uniquely weakest triangle edge as a possible indirect association; it is not a p-value procedure.

## High-order extension

For a seed pair `A-B` and background locus `C`, the frozen Diff-GWES extension estimates

```text
Delta_C = I(A;B | SHC, C=1) - I(A;B | SHC, C=0).
```

Samples are split within SHC into discovery and validation sets. Background states are permuted within SHC in validation. A passing triple is a phylogeny-conditioned statistical interaction, not proof of biochemical epistasis.

## Endpoint interpretation

Cinch retains raw scores, linkage evidence, population-structure evidence, multiplicity correction, and graph pruning as separate fields. Combining these layers into one opaque score would erase the meaning of each diagnostic and prevent exact-equivalence testing.

