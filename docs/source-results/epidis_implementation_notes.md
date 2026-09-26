# Kp pairwise EpiDis implementation notes

The score directly implements thesis equations 2-8 through 2-14. It does not call MI or JSD library functions. The public repository has no pairwise EpiDis implementation, so there is no code path to reproduce.

Kp decisions: presence is profile value >0; all other values are collapsed to 0; thesis-style Hamming-neighbour weights at tau=0.10 are normalized to sum one; epsilon=1e-27 is added literally without post-addition renormalization; both directions are retained and their mean is reported only after numerical equality testing; the comparison universe is the existing 3,279 loci and 2,609,249 pairs.

Differences from source code: the audited GitHub repository only implements major-state binary encoding and does not implement weights, EpiDis, thresholds, or pairwise scanning. Differences from the thesis are Kp-specific handling of non-positive calls as absence/non-call, the selected tau=0.10 within the thesis range, and explicit skipping of a zero-mass conditioning background.
