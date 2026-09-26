# Allele diversity vs phylogenetic mixing

This workflow builds a per-locus plot with:

- Y axis: Simpson allele diversity, calculated from positive allele calls.
- X axis: pair-weighted mean patristic distance among genomes carrying the same allele.

## Layout

```text
source_code/  analysis script and launcher
input/        allelic profile and Newick tree
output/       result table, summary, PNG and SVG
```

## Input

The profile must be tab-delimited:

```text
genome  locus_1  locus_2  ...
sample1 1        4        ...
sample2 2        4        ...
```

Allele values greater than zero are called alleles. Zero, negative values and
missing values are treated as missing. Tree tip names must match profile genome
names after removal of `.fna`, `.fa`, or `.fasta`; `GCA_` and `GCF_` accessions
are matched by accession.

## Run

```bash
bash source_code/run.sh
```

Optional parameters are passed through:

```bash
bash source_code/run.sh --core-min 0.95 --wg-min 0.01 --min-allele-count 20
```

Direct invocation:

```bash
python3 source_code/allele_diversity_phylo_mixing.py \
  --profile input/allelic_profile.tsv \
  --tree input/phylogeny.nwk \
  --output-dir output
```

## Outputs

- `gene_allele_diversity_phylo_mixing.tsv`: all per-locus metrics.
- `run_summary.tsv`: input and filtering summary.
- `allele_diversity_vs_same_allele_tree_distance.png`
- `allele_diversity_vs_same_allele_tree_distance.svg`
- `allele_diversity_vs_same_allele_tree_distance_cg_wg_panels.png`
- `allele_diversity_vs_same_allele_tree_distance_cg_wg_panels.svg`

The optional two-panel plot is enabled with `--include-wg-panel` and separates
non-overlapping locus sets:

- cgMLST-like: prevalence greater than or equal to `--core-min`.
- wgMLST accessory: prevalence from `--wg-min` up to, but excluding, `--core-min`.
