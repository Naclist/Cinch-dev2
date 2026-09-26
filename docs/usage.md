# Cinch dev2 usage

## Install and inspect

```bash
python -m venv .venv
source .venv/bin/activate
python -m pip install -e ".[test]"
pytest
cinch-dev2 doctor
cinch-dev2 list
cinch-dev2 verify
```

`doctor` reports the external `configure` and `uberBlast` modules as optional because they are required only by the legacy genome-mapping frontend.

## Inspect the staged runner

```bash
cinch-dev2 run -- --list-stages
cinch-dev2 run -- --track all --dry-run \
  --profile /data/profile.tsv \
  --pairs /data/pairs.tsv.gz \
  --positions /data/positions.tsv \
  --tree /data/tree.nwk \
  --metadata /data/metadata.tsv \
  --gff /data/reference.gff
```

The numbered main runner currently covers stages 01 through 21, including the 10b substage. Stage 22, SHC-conditioned Diff-GWES, remains a validated standalone snapshot because its frozen discovery/validation inputs were deliberately excluded from this data-free repository.

## Run selected historical stages

```bash
cinch-dev2 run -- \
  --track historical \
  --only 12-cgmlst-gene-nmi 13-gene-phylo-mixing \
  --profile /data/profile.tsv \
  --positions /data/positions.tsv \
  --tree /data/tree.nwk \
  --run-dir runs/phylo-mixing
```

Use `--dry-run` first. The runner writes a command and status manifest and stops on missing inputs, non-zero status, or missing declared outputs.

## Run weighted MI and EpiDis directly

The audited snapshots retain explicit CLIs:

```bash
python -m cinch_dev2.scripts.snapshots.weighted_mi_v2 \
  --profile /data/profile.tsv \
  --pairs /data/pairs.tsv.gz \
  --weights /data/sample_weights.tsv \
  --old-mi /data/comparator.tsv.gz \
  --outdir runs/weighted-mi-v2

python -m cinch_dev2.scripts.snapshots.epidis_pairwise_kp \
  --profile /data/profile.tsv \
  --pairs /data/pairs.tsv.gz \
  --weights /data/sample_weights.tsv \
  --mi-v2 runs/weighted-mi-v2/kp_weighted_mi_v2_all_pairs.tsv.gz \
  --old-mi /data/comparator.tsv.gz \
  --outdir runs/epidis
```

The frozen Kp scripts enforce the accepted 3,279-locus, 2,609,249-pair universe. A general-input implementation must be versioned separately rather than bypassing those guards silently.

## Rebuild repository figures

```bash
python scripts/generate_figures.py
```

The figures use only aggregate tables under `docs/source-results`. Repeated generation is deterministic.

## Refresh a deliberate source freeze

```bash
python scripts/update_freeze_manifest.py
cinch-dev2 verify
```

Do this only after an intentional code change and review. Updating hashes is not a substitute for numerical equivalence testing.

