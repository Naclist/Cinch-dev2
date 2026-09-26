# Cinch and EpiDis comparison

## What is identical

With the same binary encoding, sample weights, valid weighted empirical probabilities, and pair universe, frozen pairwise EpiDis obeys:

```text
EpiDis = sqrt(MI_nats / ln(2))
```

It therefore produces the same pair ranking as weighted MI v2. The full Kp run confirms this with Spearman correlation approximately one and maximum absolute error of about `2.87e-15` after squaring EpiDis and comparing with MI in bits.

## What differs

- EpiDis is expressed as square-root base-2 weighted JSD; MI v2 is reported in nats and also supplies NMI.
- A global IQR cutoff is not invariant under nonlinear transformation, so EpiDis and MI can select different tail counts despite identical rank.
- The older Cinch comparator uses a different neighbour definition and a 0.5 pseudocount, so it is not rank-equivalent to frozen EpiDis.
- The thesis high-order or differential-background framework is not reduced to the pairwise equivalence. Cinch's SHC-conditioned Diff-GWES is an explicit Kp extension, not a verbatim thesis statistic.

## Population structure

The thesis-style inverse-neighbour weighting reduces redundant sampling. It does not eliminate phylogenetic covariance. Cinch adds SHC-conditioned tests and within-SHC permutations as a distinct inference layer.

## Thresholds and filters

EpiDis's `Q3 + 3 IQR` rule identifies a global long tail; it is not a physical-distance or phylogenetic correction. Cinch requires distance diagnostics and SHC-aware inference before biological interpretation.

## Naming rule

Call the direct thesis-formula output `pairwise EpiDis`. Call the no-pseudocount Hamming-weighted information output `weighted MI v2`. Call the full layered process `Cinch`. Do not use the names interchangeably.

