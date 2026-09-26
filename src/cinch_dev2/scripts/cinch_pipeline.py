#!/usr/bin/env python3
"""Run the complete frozen Cinch analysis and visualization pipeline."""

from __future__ import annotations

import argparse
import json
import os
import shlex
import subprocess
import sys
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Callable


PACKAGE = Path(__file__).resolve().parents[1]
REPO = Path.cwd()
SNAPSHOTS = PACKAGE / "scripts" / "snapshots"
DEFAULT_PYTHON = PACKAGE / ".venv" / "bin" / "python"
PIPELINE_PYTHON = str(DEFAULT_PYTHON if DEFAULT_PYTHON.exists() else Path(sys.executable))


@dataclass(frozen=True)
class Context:
    repo: Path
    run: Path
    profile: Path
    pairs: Path
    positions: Path
    tree: Path
    metadata: Path
    gff: Path
    cazyme: Path
    effectors: Path
    mobile: Path
    transporters: Path


@dataclass(frozen=True)
class Stage:
    id: str
    track: str
    description: str
    command: Callable[[Context], list[str]]
    inputs: Callable[[Context], list[Path]]
    outputs: Callable[[Context], list[Path]]


def py(script: str, *args: object) -> list[str]:
    return [PIPELINE_PYTHON, str(SNAPSHOTS / script), *(str(x) for x in args)]


def out(c: Context, name: str) -> Path:
    return c.run / name


STAGES = [
    Stage(
        "01-significant-ordination",
        "historical",
        "PCA and UMAP of significant r2 profiles with strength and degree views.",
        lambda c: py("plot_r2_significant_pca_umap.py", c.pairs, "--outdir", out(c, "01_significant_ordination")),
        lambda c: [c.pairs],
        lambda c: [out(c, "01_significant_ordination/summary.tsv")],
    ),
    Stage(
        "02-distance-decay",
        "historical",
        "Order, BP, and phylogenetic distance decay, elbow thresholds, and filtered candidates.",
        lambda c: py("secondary_filter_ld_decay.py", c.pairs, "--outdir", out(c, "02_distance_decay")),
        lambda c: [c.pairs],
        lambda c: [out(c, "02_distance_decay/summary.tsv"), out(c, "02_distance_decay/elbow_thresholds.tsv")],
    ),
    Stage(
        "03-qq",
        "historical",
        "QQ plots for raw, adjusted, and distance-filtered p-values.",
        lambda c: py(
            "qq_plots_ldlike.py",
            "--all-pairs", c.pairs,
            "--filtered-pairs", out(c, "02_distance_decay/epistasis_candidate_pairs_after_elbow_filter.tsv.gz"),
            "--outdir", out(c, "03_qq"),
        ),
        lambda c: [c.pairs, out(c, "02_distance_decay/epistasis_candidate_pairs_after_elbow_filter.tsv.gz")],
        lambda c: [out(c, "03_qq/qq_summary.tsv")],
    ),
    Stage(
        "04-sign-aware-phylogeny",
        "historical",
        "Diagnostic comparison of sign-specific phylogenetic filtering rules.",
        lambda c: py("sign_specific_phylo_filter_check.py", "--pairs", c.pairs, "--outdir", out(c, "04_sign_aware_phylogeny")),
        lambda c: [c.pairs],
        lambda c: [out(c, "04_sign_aware_phylogeny/sign_specific_filter_counts.tsv")],
    ),
    Stage(
        "05-order-bp-pca",
        "historical",
        "Order-versus-BP filtered edge sets, Venn diagram, PCA, K-means, and cluster evaluation.",
        lambda c: py(
            "compare_order_bp_phylo_pca_venn.py", c.pairs,
            "--thresholds", out(c, "02_distance_decay/elbow_thresholds.tsv"),
            "--outdir", out(c, "05_order_bp_pca"),
        ),
        lambda c: [c.pairs, out(c, "02_distance_decay/elbow_thresholds.tsv")],
        lambda c: [out(c, "05_order_bp_pca/summary.tsv"), out(c, "05_order_bp_pca/venn_order_phylo_vs_bp_phylo.png")],
    ),
    Stage(
        "06-pcoa-metadata",
        "historical",
        "PCoA plus gene-by-metadata association PCA/PCoA.",
        lambda c: py(
            "pcoa_and_gene_metadata_assoc.py",
            "--order-matrix", out(c, "05_order_bp_pca/order_phylo_gene_by_gene_r2_matrix.tsv"),
            "--bp-matrix", out(c, "05_order_bp_pca/bp_phylo_gene_by_gene_r2_matrix.tsv"),
            "--profile", c.profile, "--metadata", c.metadata,
            "--outdir", out(c, "06_pcoa_metadata"),
        ),
        lambda c: [out(c, "05_order_bp_pca/order_phylo_gene_by_gene_r2_matrix.tsv"), out(c, "05_order_bp_pca/bp_phylo_gene_by_gene_r2_matrix.tsv"), c.profile, c.metadata],
        lambda c: [out(c, "06_pcoa_metadata/summary.tsv"), out(c, "06_pcoa_metadata/gene_by_metadata_nmi_matrix.tsv")],
    ),
    Stage(
        "07-rda-arrows",
        "historical",
        "Metadata arrow projections on PCA and PCoA results.",
        lambda c: py(
            "add_metadata_rda_arrows.py",
            "--assoc", out(c, "06_pcoa_metadata/gene_by_metadata_nmi_matrix.tsv"),
            "--order-pca", out(c, "05_order_bp_pca/order_phylo_pca_coords.tsv"),
            "--bp-pca", out(c, "05_order_bp_pca/bp_phylo_pca_coords.tsv"),
            "--order-pcoa", out(c, "06_pcoa_metadata/order_phylo_r2_profile_pcoa_coords.tsv"),
            "--bp-pcoa", out(c, "06_pcoa_metadata/bp_phylo_r2_profile_pcoa_coords.tsv"),
            "--outdir", out(c, "07_rda_arrows"),
        ),
        lambda c: [out(c, "06_pcoa_metadata/gene_by_metadata_nmi_matrix.tsv"), out(c, "05_order_bp_pca/order_phylo_pca_coords.tsv")],
        lambda c: [out(c, "07_rda_arrows/summary.tsv")],
    ),
    Stage(
        "08-nmds-constrained-rda",
        "historical",
        "NMDS, genuine constrained RDA, K-means views, vectors, and diagnostics.",
        lambda c: py(
            "nmds_and_constrained_rda.py",
            "--order-matrix", out(c, "05_order_bp_pca/order_phylo_gene_by_gene_r2_matrix.tsv"),
            "--bp-matrix", out(c, "05_order_bp_pca/bp_phylo_gene_by_gene_r2_matrix.tsv"),
            "--assoc", out(c, "06_pcoa_metadata/gene_by_metadata_expanded_assoc_matrix.tsv"),
            "--outdir", out(c, "08_nmds_constrained_rda"),
        ),
        lambda c: [out(c, "05_order_bp_pca/order_phylo_gene_by_gene_r2_matrix.tsv"), out(c, "06_pcoa_metadata/gene_by_metadata_expanded_assoc_matrix.tsv")],
        lambda c: [out(c, "08_nmds_constrained_rda/summary.tsv")],
    ),
    Stage(
        "09-population-pca",
        "historical",
        "Genome-by-gene presence PCA, variance curves, and selected metadata association.",
        lambda c: py("profile_population_pca.py", "--profile", c.profile, "--metadata", c.metadata, "--outdir", out(c, "09_population_pca")),
        lambda c: [c.profile, c.metadata],
        lambda c: [out(c, "09_population_pca/summary.tsv"), out(c, "09_population_pca/genome_presence_pca_scree.png")],
    ),
    Stage(
        "10-cluster-annotation",
        "historical",
        "Selected PCA-cluster annotation and Phandango-like tree/presence heatmaps.",
        lambda c: py(
            "cluster3_annotation_and_tree_heatmap.py",
            "--outdir", out(c, "10_cluster_annotation"),
            "--gff", c.gff, "--profile", c.profile, "--tree", c.tree,
            "--order-coords", out(c, "05_order_bp_pca/order_phylo_pca_coords.tsv"),
            "--bp-coords", out(c, "05_order_bp_pca/bp_phylo_pca_coords.tsv"),
            "--cazyme", c.cazyme, "--effectors", c.effectors,
            "--mobile-elements", c.mobile, "--transporters", c.transporters,
        ),
        lambda c: [c.gff, c.profile, c.tree, out(c, "05_order_bp_pca/order_phylo_pca_coords.tsv"), c.cazyme, c.effectors, c.mobile, c.transporters],
        lambda c: [out(c, "10_cluster_annotation/selected_cluster_annotation_summary.tsv")],
    ),
    Stage(
        "10b-cluster-geography",
        "historical",
        "India/Taiwan enrichment and NMI for the selected PCA cluster, with a comparison visual.",
        lambda c: py(
            "selected_cluster_geography_association.py",
            "--cluster-table", out(c, "10_cluster_annotation/order_phylo_pca_cluster3_gene_annotations_with_prevalence.tsv"),
            "--profile", c.profile, "--metadata", c.metadata,
            "--nmi", out(c, "06_pcoa_metadata/gene_by_metadata_nmi_matrix.tsv"),
            "--levels", "India", "Taiwan", "--outdir", out(c, "10b_cluster_geography"),
        ),
        lambda c: [out(c, "10_cluster_annotation/order_phylo_pca_cluster3_gene_annotations_with_prevalence.tsv"), c.profile, c.metadata, out(c, "06_pcoa_metadata/gene_by_metadata_nmi_matrix.tsv")],
        lambda c: [out(c, "10b_cluster_geography/selected_cluster_india_association.tsv"), out(c, "10b_cluster_geography/selected_cluster_taiwan_association.tsv"), out(c, "10b_cluster_geography/selected_cluster_geography_summary.png")],
    ),
    Stage(
        "11-cgmlst-allelic-scan",
        "historical",
        "cgMLST allele-state gene-pair NMI, chi-square approximation, and allele-pair detail.",
        lambda c: py(
            "cgmlst_allelic_scan.py", "--profile", c.profile, "--metadata", c.metadata,
            "--tree", c.tree, "--outdir", out(c, "11_cgmlst_allelic_scan"),
        ),
        lambda c: [c.profile, c.metadata, c.tree],
        lambda c: [out(c, "11_cgmlst_allelic_scan/summary.tsv")],
    ),
    Stage(
        "12-cgmlst-gene-nmi",
        "historical",
        "Unfiltered cgMLST gene-by-gene NMI, approximate chi-square, distance decay, and top links.",
        lambda c: py(
            "cgmlst_gene_gene_nmi_distance.py", "--profile", c.profile,
            "--positions", c.positions, "--outdir", out(c, "12_cgmlst_gene_nmi"),
        ),
        lambda c: [c.profile, c.positions],
        lambda c: [out(c, "12_cgmlst_gene_nmi/summary.tsv"), out(c, "12_cgmlst_gene_nmi/top20_flattened_genome_links.png")],
    ),
    Stage(
        "13-gene-phylo-mixing",
        "historical",
        "Endpoint plot: per-gene Simpson allele diversity versus same-allele mean tree distance.",
        lambda c: py("gene_phylo_mixing_metrics.py", "--profile", c.profile, "--tree", c.tree, "--outdir", out(c, "13_gene_phylo_mixing")),
        lambda c: [c.profile, c.tree],
        lambda c: [out(c, "13_gene_phylo_mixing/summary.tsv"), out(c, "13_gene_phylo_mixing/allele_diversity_vs_same_allele_tree_distance.png"), out(c, "13_gene_phylo_mixing/lineage_strong_genes.tsv")],
    ),
    Stage(
        "14-cgmlst-phylo-filtered",
        "historical",
        "Repeat cgMLST gene-pair NMI after excluding every pair containing a lineage-strong gene.",
        lambda c: py(
            "cgmlst_gene_gene_nmi_distance.py", "--profile", c.profile,
            "--positions", c.positions, "--exclude-genes", out(c, "13_gene_phylo_mixing/lineage_strong_genes.tsv"),
            "--outdir", out(c, "14_cgmlst_phylo_filtered"),
        ),
        lambda c: [c.profile, c.positions, out(c, "13_gene_phylo_mixing/lineage_strong_genes.tsv")],
        lambda c: [out(c, "14_cgmlst_phylo_filtered/summary.tsv"), out(c, "14_cgmlst_phylo_filtered/top20_flattened_genome_links.png")],
    ),
    Stage(
        "15-profile-similarity",
        "extension",
        "wg/cg isolate similarities and sample-reweighting tables.",
        lambda c: py("profile_wg_cg_similarity.py", "--profile", c.profile, "--outdir", out(c, "15_profile_similarity")),
        lambda c: [c.profile],
        lambda c: [out(c, "15_profile_similarity/spydrpick_style_sample_weights.tsv"), out(c, "15_profile_similarity/wg_cg_similarity_matrices.npz")],
    ),
    Stage(
        "16-old-weighted-mi",
        "extension",
        "Existing pseudocount-weighted MI/NMI comparator and ARACNE candidates.",
        lambda c: py(
            "spydrpick_locus_wg_nmi.py", "--profile", c.profile, "--pairs", c.pairs,
            "--weights", out(c, "15_profile_similarity/spydrpick_style_sample_weights.tsv"),
            "--gff", c.gff, "--outdir", out(c, "16_old_weighted_mi"),
        ),
        lambda c: [c.profile, c.pairs, out(c, "15_profile_similarity/spydrpick_style_sample_weights.tsv"), c.gff],
        lambda c: [out(c, "16_old_weighted_mi/wgmlst_all_weighted_mi_nmi_pairs.tsv.gz")],
    ),
    Stage(
        "17-mi-distance-visuals",
        "extension",
        "MI/NMI distance comparisons, elbows, outliers, and direct/indirect views.",
        lambda c: py(
            "plot_wgmlst_mi_nmi_distance.py", "--pairs", out(c, "16_old_weighted_mi/wgmlst_all_weighted_mi_nmi_pairs.tsv.gz"),
            "--gff", c.gff, "--outdir", out(c, "17_mi_distance_visuals"),
        ),
        lambda c: [out(c, "16_old_weighted_mi/wgmlst_all_weighted_mi_nmi_pairs.tsv.gz"), c.gff],
        lambda c: [out(c, "17_mi_distance_visuals/mi_nmi_distance_summary.tsv")],
    ),
    Stage(
        "18-weighted-mi-v2",
        "extension",
        "Thesis-style Hamming-neighbour sample-reweighted MI v2.",
        lambda c: py(
            "weighted_mi_v2.py", "--profile", c.profile, "--pairs", c.pairs,
            "--weights", out(c, "15_profile_similarity/spydrpick_style_sample_weights.tsv"),
            "--old-mi", out(c, "16_old_weighted_mi/wgmlst_all_weighted_mi_nmi_pairs.tsv.gz"),
            "--outdir", out(c, "18_weighted_mi_v2"),
        ),
        lambda c: [c.profile, c.pairs, out(c, "15_profile_similarity/spydrpick_style_sample_weights.tsv"), out(c, "16_old_weighted_mi/wgmlst_all_weighted_mi_nmi_pairs.tsv.gz")],
        lambda c: [out(c, "18_weighted_mi_v2/validation_summary.json"), out(c, "18_weighted_mi_v2/kp_weighted_mi_v2_all_pairs.tsv.gz")],
    ),
    Stage(
        "19-pairwise-epidis",
        "extension",
        "Thesis-defined pairwise EpiDis and direct MI-equivalence audit.",
        lambda c: py(
            "epidis_pairwise_kp.py", "--profile", c.profile, "--pairs", c.pairs,
            "--weights", out(c, "15_profile_similarity/spydrpick_style_sample_weights.tsv"),
            "--mi-v2", out(c, "18_weighted_mi_v2/kp_weighted_mi_v2_all_pairs.tsv.gz"),
            "--old-mi", out(c, "16_old_weighted_mi/wgmlst_all_weighted_mi_nmi_pairs.tsv.gz"),
            "--outdir", out(c, "19_pairwise_epidis"),
        ),
        lambda c: [c.profile, c.pairs, out(c, "15_profile_similarity/spydrpick_style_sample_weights.tsv"), out(c, "18_weighted_mi_v2/kp_weighted_mi_v2_all_pairs.tsv.gz")],
        lambda c: [out(c, "19_pairwise_epidis/summary.tsv"), out(c, "19_pairwise_epidis/kp_pairwise_epidis_all_pairs.tsv.gz")],
    ),
    Stage(
        "20-epidis-visuals",
        "extension",
        "EpiDis distribution, MI comparison, and order/BP distance visualizations.",
        lambda c: py(
            "plot_epidis_pairwise_kp.py", "--input", out(c, "19_pairwise_epidis/kp_pairwise_epidis_all_pairs.tsv.gz"),
            "--threshold", out(c, "19_pairwise_epidis/epidis_iqr_threshold.tsv"),
            "--outdir", out(c, "20_epidis_visuals"),
        ),
        lambda c: [out(c, "19_pairwise_epidis/kp_pairwise_epidis_all_pairs.tsv.gz"), out(c, "19_pairwise_epidis/epidis_iqr_threshold.tsv")],
        lambda c: [out(c, "20_epidis_visuals/kp_epidis_overview.png"), out(c, "20_epidis_visuals/kp_epidis_vs_previous_mi.png")],
    ),
    Stage(
        "21-shc-correction",
        "extension",
        "Stable hierarchical cluster construction and within-SHC permutation filtering.",
        lambda c: py(
            "hiercc_shc_phylo_correct_mi.py", "--profile", c.profile,
            "--similarity", out(c, "15_profile_similarity/wg_cg_similarity_matrices.npz"),
            "--weights", out(c, "15_profile_similarity/spydrpick_style_sample_weights.tsv"),
            "--candidate-dir", out(c, "17_mi_distance_visuals"),
            "--outdir", out(c, "21_shc_correction"),
        ),
        lambda c: [c.profile, out(c, "15_profile_similarity/wg_cg_similarity_matrices.npz"), out(c, "17_mi_distance_visuals")],
        lambda c: [out(c, "21_shc_correction/genome_shc_assignments.tsv"), out(c, "21_shc_correction/shc_filter_summary.tsv")],
    ),
]


def default_context(args: argparse.Namespace) -> Context:
    run = args.run_dir.resolve()
    return Context(
        repo=REPO,
        run=run,
        profile=args.profile.resolve(),
        pairs=args.pairs.resolve(),
        positions=args.positions.resolve(),
        tree=args.tree.resolve(),
        metadata=args.metadata.resolve(),
        gff=args.gff.resolve(),
        cazyme=args.cazyme.resolve(),
        effectors=args.effectors.resolve(),
        mobile=args.mobile_elements.resolve(),
        transporters=args.transporters.resolve(),
    )


def select_stages(args: argparse.Namespace) -> list[Stage]:
    selected = [s for s in STAGES if args.track == "all" or s.track == args.track]
    if args.only:
        requested = set(args.only)
        selected = [s for s in selected if s.id in requested]
        unknown = requested - {s.id for s in STAGES}
        if unknown:
            raise SystemExit(f"Unknown stage(s): {', '.join(sorted(unknown))}")
    ids = [s.id for s in selected]
    if args.from_stage:
        if args.from_stage not in ids:
            raise SystemExit(f"--from stage is outside selected track: {args.from_stage}")
        selected = selected[ids.index(args.from_stage):]
        ids = [s.id for s in selected]
    if args.through:
        if args.through not in ids:
            raise SystemExit(f"--through stage is outside selected range: {args.through}")
        selected = selected[: ids.index(args.through) + 1]
    return selected


def write_manifest(path: Path, payload: dict) -> None:
    path.write_text(json.dumps(payload, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")


def main() -> None:
    global PIPELINE_PYTHON
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--track", choices=["historical", "extension", "all"], default="all")
    parser.add_argument("--run-dir", type=Path, default=PACKAGE / "runs" / datetime.now().strftime("%Y%m%d_%H%M%S"))
    parser.add_argument("--from", dest="from_stage")
    parser.add_argument("--through")
    parser.add_argument("--only", nargs="+")
    parser.add_argument("--dry-run", action="store_true")
    parser.add_argument("--force", action="store_true")
    parser.add_argument("--profile", type=Path, default=REPO / "Leca/source/profile/0427_1.profile")
    parser.add_argument("--pairs", type=Path, default=REPO / "Leca/source/profile/0427_1.ldlike_pairs.tsv.gz")
    parser.add_argument("--positions", type=Path, default=REPO / "Leca/source/positions/test0427_1_GCA_905232925.fna.positions.tsv")
    parser.add_argument("--tree", type=Path, default=REPO / "Leca/source/tree/Kp.labelled.nwk")
    parser.add_argument("--metadata", type=Path, default=REPO / "Leca/source/metadata/Kp.mp.merge")
    parser.add_argument("--gff", type=Path, default=REPO / "Leca/source/annotation/klebsiella_reference.gff")
    parser.add_argument("--cazyme", type=Path, default=REPO / "Leca/source/annotation/klebsiella_reference.cazyme.txt")
    parser.add_argument("--effectors", type=Path, default=REPO / "Leca/source/annotation/klebsiella_reference.effectors.txt")
    parser.add_argument("--mobile-elements", type=Path, default=REPO / "Leca/source/annotation/klebsiella_reference.mobile_elements.txt")
    parser.add_argument("--transporters", type=Path, default=REPO / "Leca/source/annotation/klebsiella_reference.transporters.txt")
    parser.add_argument("--list-stages", action="store_true")
    parser.add_argument("--python", type=Path, default=DEFAULT_PYTHON if DEFAULT_PYTHON.exists() else Path(sys.executable))
    args = parser.parse_args()
    # Preserve a virtual-environment interpreter symlink; resolving it would
    # bypass the venv's site-packages and execute the base interpreter instead.
    PIPELINE_PYTHON = str(args.python.absolute())

    if args.list_stages:
        for stage in STAGES:
            print(f"{stage.id:30s} {stage.track:10s} {stage.description}")
        return

    context = default_context(args)
    stages = select_stages(args)
    if not args.python.exists():
        raise SystemExit(f"Python interpreter not found: {args.python}; run scripts/bootstrap_env.py")
    if args.dry_run:
        for stage in stages:
            print(f"[{stage.id}] {shlex.join(stage.command(context))}")
        return

    context.run.mkdir(parents=True, exist_ok=True)
    log_dir = context.run / "logs"
    log_dir.mkdir(exist_ok=True)
    manifest = {
        "schema": "cinch.pipeline-run",
        "schema_version": "1.0.0",
        "started_at": datetime.now(timezone.utc).isoformat(),
        "track": args.track,
        "run_dir": str(context.run),
        "stages": [],
    }
    manifest_path = context.run / "run_manifest.json"

    for stage in stages:
        outputs = stage.outputs(context)
        if outputs and all(path.exists() for path in outputs) and not args.force:
            print(f"SKIP {stage.id}: expected outputs already exist")
            manifest["stages"].append({"id": stage.id, "status": "skipped_existing"})
            write_manifest(manifest_path, manifest)
            continue
        missing = [path for path in stage.inputs(context) if not path.exists()]
        if missing:
            manifest["stages"].append({"id": stage.id, "status": "blocked_missing_input", "missing": [str(p) for p in missing]})
            write_manifest(manifest_path, manifest)
            raise SystemExit(f"{stage.id} missing input(s): {', '.join(str(p) for p in missing)}")

        command = stage.command(context)
        print(f"RUN  {stage.id}")
        started = datetime.now(timezone.utc).isoformat()
        env = os.environ.copy()
        env["MPLBACKEND"] = "Agg"
        log_path = log_dir / f"{stage.id}.log"
        with log_path.open("w", encoding="utf-8") as log:
            result = subprocess.run(command, cwd=context.repo, env=env, stdout=log, stderr=subprocess.STDOUT, check=False)
        record = {
            "id": stage.id,
            "track": stage.track,
            "status": "completed" if result.returncode == 0 else "failed",
            "started_at": started,
            "completed_at": datetime.now(timezone.utc).isoformat(),
            "returncode": result.returncode,
            "command": command,
            "log": str(log_path),
            "outputs": [str(p) for p in outputs],
        }
        manifest["stages"].append(record)
        write_manifest(manifest_path, manifest)
        if result.returncode != 0:
            raise SystemExit(f"{stage.id} failed; inspect {log_path}")
        missing_outputs = [path for path in outputs if not path.exists()]
        if missing_outputs:
            raise SystemExit(f"{stage.id} returned success but missed outputs: {', '.join(str(p) for p in missing_outputs)}")

    manifest["completed_at"] = datetime.now(timezone.utc).isoformat()
    write_manifest(manifest_path, manifest)
    print(f"COMPLETE {context.run}")


if __name__ == "__main__":
    main()
