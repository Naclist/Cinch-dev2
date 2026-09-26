# Data inclusion policy

## Included

- source code required to reconstruct the local Cinch workflow;
- small aggregate validation TSV and JSON files;
- algorithm audits, field dictionaries, frozen specifications and provenance hashes;
- deterministic documentation figures derived from those aggregate files.

## Excluded

- the local virtual environment and Python caches;
- smoke-test and demonstration outputs;
- raw or assembled genomes;
- full wgMLST/cgMLST profiles;
- complete pairwise result matrices;
- synthetic toy profiles and toy result directories;
- private conversation exports and local session manifests;
- unrelated Leca projects and historical invalid/partial outputs.

The exclusions keep the repository small and avoid publishing private or non-redistributable material. They do not change the registered scientific stages. Exact historical reproduction requires the original frozen inputs whose hashes are recorded in the provenance documents.

