---
name: cinch
description: Freeze, audit, run, compare, and hand off the Cinch phylogeny-aware bacterial epistasis workflow built from wg/cgMLST profiles. Use for Cinch, NMI/MI-based epistasis, weighted MI v2, thesis-defined pairwise EpiDis comparison, physical-distance filtering, SHC-conditioned permutation, Diff-GWES, result interpretation, or importing the structured Cinch conversation archive into another AI agent.
---

# Cinch

Treat Cinch as a staged statistical-candidate workflow, not a molecular-causality test.

## Start here

1. Read `references/methodology.md` before changing a statistic or filter.
2. Read `references/usage.md` before running a stage.
3. Read `references/historical-coverage.md` to map every pre-EpiDis result family.
4. Read `references/results-documentary.md` before citing the frozen Kp results.
5. Read `references/epidis-comparison.md` when comparing Cinch with EpiDis.
6. Load `conversations/cinch-conversations.json` when reconstructing prior decisions.
7. Run `python3 scripts/cinch.py verify` before trusting the snapshots.

## Workflow contract

Execute stages in this order unless performing an explicitly scoped audit:

1. Build the wg/cgMLST profile and pair universe.
2. Compute the frozen comparator and thesis-style weighted MI v2.
3. Compute pairwise EpiDis on the identical pair universe and weights.
4. Quantify physical-distance decay before interpreting outliers.
5. Apply SHC-aware conditional testing and multiple-testing correction.
6. Apply ARACNE only as network pruning; never treat it as a significance test.
7. Run SHC-conditioned Diff-GWES only for background-dependent effects.
8. Report statistical candidates and unresolved confounding explicitly.

## Integrity rules

- Keep `profile > 0` as presence and `profile <= 0` as absence/non-callable for the frozen Kp workflow.
- Do not silently reinterpret zero, missingness, paralogy, or assembly failure.
- Use an identical pair universe for metric comparisons.
- Do not call weighted MI v2 EpiDis.
- Record that pairwise EpiDis is `sqrt(MI_bits)` under the frozen valid-probability construction.
- Keep raw scores, distance filtering, SHC inference, and network pruning as separate layers.
- Do not label a pair or community biological epistasis without experimental validation.
- Preserve source hashes and provenance when editing or rerunning a stage.

## Commands

Use `python3 scripts/cinch.py list` to inspect stages, `show STAGE` for details, and
`verify` to check frozen snapshot hashes. Use `export-chats` to rebuild the structured
conversation archive from local Codex session JSONL files.

Use `python3 scripts/cinch_pipeline.py --track historical` to reproduce all analysis
and visualization stages through the allele-diversity/phylogenetic-mixing endpoint.
Use `--track extension` for the later weighted MI, EpiDis, and SHC stages.
