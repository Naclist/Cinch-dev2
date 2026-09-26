# Historical result coverage through the gene phylogenetic-mixing endpoint

The comprehensive runner is `scripts/cinch_pipeline.py`. Its `historical` track covers every major executable result family discussed from the first Leca/Cinch analysis through the requested plot whose Y axis is per-gene allele diversity and whose X axis is mean tree distance among genomes sharing the same allele.

| Stage | Result family | Primary visual/table |
|---|---|---|
| 01 | Significant r2 profile ordination | PCA/UMAP by strength and degree |
| 02 | Physical and phylogenetic distance decay | Scatter/hexbin trends, elbows, filtered candidates |
| 03 | Calibration diagnostics | Raw, adjusted, and filtered QQ plots |
| 04 | Positive/negative sign-aware phylogeny diagnostic | Sign-specific filter comparison |
| 05 | Order versus BP filtering | Two PCAs, K-means views/evaluation, Venn diagram |
| 06 | Gene-metadata association | r2-profile PCoA and gene-metadata PCA/PCoA |
| 07 | Metadata projection | PCA/PCoA arrow overlays |
| 08 | Constrained ordination | NMDS, constrained RDA, K-means, vector tables |
| 09 | Population structure | Genome-presence PCA, scree and cumulative variance |
| 10 | Cluster interpretation | Annotation tables and Phandango-like tree heatmaps |
| 10b | Selected-cluster geography | India/Taiwan enrichment tables, gene-geography NMI, comparison visual |
| 11 | Allele-state association prototype | cgMLST allelic NMI and approximate inference |
| 12 | Gene-level cgMLST association | Gene-gene NMI, distance decay, flattened-genome top links |
| 13 | Phylogenetic mixing endpoint | Simpson diversity versus same-allele mean tree distance |
| 14 | Required lineage-gene exclusion | Repeat stage 12 excluding every pair containing a lineage-strong gene |

The `extension` track then adds the later sample reweighting, original weighted MI/NMI comparator, MI-distance visualization, weighted MI v2, pairwise EpiDis, EpiDis visualization, and SHC correction.

## Endpoint definition

For each locus, positive calls define observed alleles. Y is Simpson diversity `1 - sum_a p_a^2`. X is a pair-count-weighted mean patristic distance over genome pairs sharing an allele whose count meets `--min-allele-count`. The frozen defaults evaluate loci with prevalence at least 0.95 and alleles observed in at least 20 genomes.

The lower 10% of X is labeled `lineage_strong`, the upper 10% is labeled `mixed_high`, and the remaining loci are intermediate. The labels are empirical quantile classes, not universal biological categories.

The Kp frozen result evaluated 3,632 core loci, labeled 364 lineage-strong loci and 364 high-mixing loci, and used an overall mean genome-pair tree distance of approximately 0.006933. The lineage-strong cutoff was approximately 0.0003715 and the high-mixing cutoff approximately 0.0047316.
