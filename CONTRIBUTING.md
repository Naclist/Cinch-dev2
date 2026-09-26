# Contributing

Cinch dev2 is a scientific-computing research preview. Correctness takes priority over speed.

Before proposing a change:

1. identify whether it is exact-equivalent engineering or a statistical-method change;
2. preserve the frozen pair universe, state semantics, weights, thresholds, random seeds, and permutation universe unless the change is explicitly versioned as a new method;
3. add a unit test or golden comparison;
4. run `pytest` and `cinch-dev2 verify`;
5. regenerate figures only with `python scripts/generate_figures.py`;
6. disclose any numerical difference, even if rankings appear unchanged.

Do not present MI, EpiDis, an IQR outlier, an ARACNE-retained edge, or a graph community as causal molecular epistasis without independent experimental evidence.

