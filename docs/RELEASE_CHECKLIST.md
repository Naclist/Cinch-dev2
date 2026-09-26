# Release checklist

1. `python -m pip install -e ".[test]"`
2. `pytest`
3. `cinch-dev2 doctor`
4. `python scripts/update_freeze_manifest.py`
5. `cinch-dev2 verify`
6. `python scripts/generate_figures.py`
7. confirm `git diff --exit-code -- docs/figures site/assets` after a second figure-generation run;
8. inspect every PNG and the responsive static site;
9. confirm no raw genomes, profiles, pair matrices, credentials, raw conversations, caches or virtual environments are tracked;
10. record the commit hash and create the remote repository as private until visibility and licensing are chosen explicitly.

