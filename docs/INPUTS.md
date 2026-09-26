# Inputs and state semantics

## Required profile inputs

The profile-based workflow expects:

- a tab-separated genome-by-locus allele profile;
- a stable sample identifier in the first column;
- a frozen locus-pair table for exact historical reproduction;
- locus positions for order/base-pair diagnostics;
- a phylogenetic tree for tree-distance and mixing analyses;
- optional metadata and annotation tables.

Every production run should provide explicit paths. The legacy defaults document the accepted Kp analysis environment but are not portable data discovery.

## Frozen presence encoding

The accepted Kp implementation uses:

```text
profile value > 0  -> presence 1
profile value <= 0 -> absence/non-callable 0
missing value      -> filled with 0
```

This encoding cannot distinguish biological absence, low coverage, assembly failure, allele-calling failure, or other non-callable states. A future state model may separate these values, but doing so changes the statistical input and requires a separately versioned validation.

## Multiallelic encoding

cgMLST allele identifiers are nominal categories. Their integer values do not represent numeric distance, SNP distance, or branch length. Dense integer recoding is permitted as an exact engineering optimization; ordinal arithmetic on allele IDs is not.

## External WGS frontend

The retained `Cinch_v8.py` frontend uses external `configure` and `uberBlast` modules. They are not redistributed. Cinch dev2 does not replace them with a limited internal exact-seed search because such a replacement would alter sensitivity and absence calls.

