# Frozen SHC-conditioned Diff-GWES specification, version 1.0

## Target estimand

For a seed pair `(A,B)` and binary background locus `C`, estimate whether the strength of the `A-B` association differs between `C=1` and `C=0` after conditioning on SHC.

Within background state `c`, the conditional information is:

```text
I_c(A;B | SHC) = sum_h P(h | c, informative) I(A;B | h,c)
E_c = sqrt(I_c / ln(2))
Delta(A,B;C) = E_1 - E_0
```

`Delta>0` means enhancement in the `C=1` background; `Delta<0` means attenuation. This is a Leca extension of thesis Diff-GWES, not a verbatim thesis statistic.

## Independent discovery and validation

Samples are deterministically split 50:50 inside every SHC. Discovery scans all eligible backgrounds and retains `abs(Delta)>=0.15`, capped at the 25 strongest backgrounds per seed pair. Validation recomputes effects on held-out samples only.

## Eligibility

- Background whole-data prevalence: 0.05 to 0.95.
- At least 25 validation samples in each background state.
- Within each state, at least two informative SHCs.
- An informative SHC-state cell contains at least three samples and both `A` and `B` vary.
- Backgrounds sharing a 100 kb reference block with `A` or `B` are excluded when all loci are placed.

## Population-preserving block null

For each validation permutation, background `C` is permuted within SHC. Loci assigned to the same reference block share deterministic permutation streams. `A` and `B` remain fixed. This preserves SHC-specific background frequency and block-correlated permutation structure while breaking background modification of the seed association.

The two-sided empirical p-value is `(1 + count(abs(Delta_null)>=abs(Delta_obs))) / (R+1)`, with `R=999`. BH correction is applied across independently selected, validation-eligible hypotheses.

Unplaced loci receive singleton block IDs. In the frozen Kp input, 952/3,279 loci map to the single-reference flattened coordinate map; block preservation is therefore complete only for those loci. Unplaced accessory loci retain SHC-specific frequency under permutation but have no defensible local-LD block preservation.

## Toy benchmark

Twenty synthetic `(A,B,C)` triples are appended to a copy of the profile. Within every sufficiently large SHC and independently within discovery and validation samples, `B=A` for `C=1`; the four `A-B` cells are balanced for `C=0`. Twenty cross-triple backgrounds are predeclared as negative controls before validation. Synthetic A-B pairs are assigned a declared simulated order distance of 300 and BP distance of 300,000; these are not observed genomic coordinates.

## Interpretation boundary

Passing triples are phylogeny-conditioned statistical interactions. They are not proof of molecular epistasis or causality. Presence, absence, missingness, and non-callable states remain collapsed by the source profile contract (`value>0` is 1; otherwise 0).
