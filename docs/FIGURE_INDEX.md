# Scientific figure index

This index records the result-bearing figures bundled with Cinch dev2. The repository does not treat a workflow diagram as scientific evidence.

| ID | Scientific question | Frozen source | Scope and interpretation boundary |
|---|---|---|---|
| SCI01 | What is the EpiDis distribution, its numerical relation to weighted MI v2, and its dependence on gene-order/physical distance? | `Leca/results/gene_gene_wgmlst/epidis_pairwise_kp/figures/kp_epidis_overview.png` | 1,000-genome Kp pair universe. Upper-tail membership alone is not evidence of epistasis. |
| SCI02 | How does EpiDis compare with the previous weighted MI and NMI workflow? | `Leca/results/gene_gene_wgmlst/epidis_pairwise_kp/figures/kp_epidis_vs_previous_mi.png` | Method comparison over the same frozen pair universe; correlation does not imply equivalence of thresholds. |
| SCI03 | How do weighted MI and weighted NMI change with gene-order and base-pair distance? | `Leca/results/gene_gene_wgmlst/spydrpick_locus_nmi/mi_nmi_distance_comparison/weighted_mi_nmi_vs_order_and_bp_distance.png` | Descriptive distance comparison. Highlighted tails use metric-specific fences and elbows. |
| SCI04 | How is the smallest stable cluster selected, and how do edges change after within-SHC conditioning? | `Leca/results/gene_gene_wgmlst/spydrpick_locus_nmi/hiercc_shc_correction/hiercc_shc_correction_diagnostics.png` | Calibration and attrition diagnostic. Conditional MI is not a causal effect estimate. |
| SCI05 | How does per-gene allele diversity relate to phylogenetic mixing among genomes sharing the same allele? | `Leca/results/gene_gene_cgmlst/phylo_mixing/allele_diversity_vs_same_allele_tree_distance.png` | x = mean phylogenetic distance among same-allele genomes; y = per-gene Simpson allele diversity. Regimes are descriptive classes. |
| SCI06 | Where are hypotheses removed in SHC-conditioned Diff-GWES, and what survives held-out validation? | `Leca/workflows/shc_diff_gwes/output/figures/shc_diff_gwes_pipeline_attrition.png` | Natural Kp and planted positive-control counts are displayed together but are not the same evidence class. |
| SCI07 | What network topology remains among true-data SHC-strict conditional-MI edges? | `Leca/workflows/shc_diff_gwes/output/figures/true_shc_community_network.png` | Community labels summarize topology; they are not additional hypothesis tests or functional annotations. |
| SCI08 | How does mean MIraw vary across the three public Cinch distance diagnostics? | Public `Naclist/Cinch` commit `98c41ef081be4face6b32cb61b557964537cb353`, `FIG06_DISTANCE_MEAN_CURVES.png` | Three panels are phylogenetic/tree distance, gene-order distance, and physical bp distance. The third panel is not a SNP-difference count. |
| SCI09 | What are the pair-density envelopes, coordinate coverage, and order-bp agreement for those distances? | Public `Naclist/Cinch` commit `98c41ef081be4face6b32cb61b557964537cb353`, `FIG06_DISTANCE_COMPARISON.png` | Public SPN534 development snapshot; descriptive correlations are not independent validation or causal evidence. |

## Distance-name reconciliation

The public Cinch source names the three quantities `phylo_dist`, `order_dist`, and `bp_dist`. In this repository:

- “tree-dist” maps to the public phylogenetic-distance panel;
- “order-dist” maps directly to the public gene-order panel;
- the requested “SNP-dist” label is not used for the third panel because the underlying public variable is physical base-pair coordinate distance, not the number of SNP differences.

Both public composites are retained so none of the three panels is lost.

## Generated validation summaries

- `FIG02_EPIDIS_VALIDATION` is regenerated from `docs/source-results/epidis_correlations_with_mi.tsv`.
- `FIG03_DIFF_GWES_VALIDATION` is regenerated from `docs/source-results/shc_diff_gwes_validation_summary.json` and remains available as an implementation-validation figure, but it is not the primary biological result display.

The copied scientific figures are frozen artifacts. Regeneration belongs to their source analysis scripts in the local Leca workflow snapshots; `scripts/generate_figures.py` only regenerates compact validation summaries from retained aggregate tables.

Raster/SVG bytes are not asserted equal across operating systems, Python versions, Matplotlib versions, or font-rendering libraries. The CI contract checks the numerical tests, frozen source hashes, successful figure generation, and presence of every expected output; scientific values come from the retained source tables rather than image bytes.
