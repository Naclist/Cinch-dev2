#!/usr/bin/env python3
import argparse
import re
from dataclasses import dataclass, field
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib.colors import ListedColormap


@dataclass
class Node:
    name: str = ""
    length: float = 0.0
    children: list["Node"] = field(default_factory=list)
    x: float = 0.0
    y: float = 0.0

    @property
    def is_leaf(self) -> bool:
        return not self.children


def gene_number(gene: str) -> int:
    match = re.search(r"ref_gene_(\d+)$", str(gene))
    return int(match.group(1)) if match else 10**12


def clean_tip_name(name: str) -> str:
    return str(name).removesuffix(".fna")


def parse_newick(text: str) -> Node:
    text = text.strip().rstrip(";")
    idx = 0

    def parse_label_and_length() -> tuple[str, float]:
        nonlocal idx
        start = idx
        while idx < len(text) and text[idx] not in ",()":
            idx += 1
        token = text[start:idx].strip()
        if not token:
            return "", 0.0
        if ":" in token:
            name, length = token.rsplit(":", 1)
            try:
                return name, float(length)
            except ValueError:
                return name, 0.0
        return token, 0.0

    def parse_subtree() -> Node:
        nonlocal idx
        if text[idx] == "(":
            idx += 1
            children = []
            while True:
                children.append(parse_subtree())
                if idx >= len(text):
                    break
                if text[idx] == ",":
                    idx += 1
                    continue
                if text[idx] == ")":
                    idx += 1
                    break
            name, length = parse_label_and_length()
            return Node(name=name, length=length, children=children)
        name, length = parse_label_and_length()
        return Node(name=name, length=length)

    return parse_subtree()


def prune_tree(node: Node, keep: set[str]) -> Node | None:
    if node.is_leaf:
        return Node(name=clean_tip_name(node.name), length=node.length) if clean_tip_name(node.name) in keep else None
    children = [child for child in (prune_tree(child, keep) for child in node.children) if child is not None]
    if not children:
        return None
    return Node(name=node.name, length=node.length, children=children)


def assign_tree_layout(node: Node, leaf_order: list[str], parent_x: float = 0.0) -> None:
    node.x = parent_x + max(node.length, 0.0)
    if node.is_leaf:
        node.y = float(leaf_order.index(clean_tip_name(node.name)))
        return
    for child in node.children:
        assign_tree_layout(child, leaf_order, node.x)
    node.y = float(np.mean([child.y for child in node.children]))


def draw_tree(ax: plt.Axes, node: Node) -> None:
    if node.is_leaf:
        return
    child_ys = [child.y for child in node.children]
    ax.plot([node.x, node.x], [min(child_ys), max(child_ys)], color="#303030", lw=0.35)
    for child in node.children:
        ax.plot([node.x, child.x], [child.y, child.y], color="#303030", lw=0.35)
        draw_tree(ax, child)


def read_ref_gene_md5(gff: Path, prefix: str) -> pd.DataFrame:
    rows = []
    pattern = re.compile(r"ref_gene_(\d+)")
    with gff.open() as handle:
        for line in handle:
            if not line.startswith("##sequence-region"):
                continue
            gene_match = pattern.search(line)
            md5_match = re.search(r"md5=([0-9a-fA-F]+)", line)
            if gene_match and md5_match:
                rows.append(
                    {
                        "gene": f"{prefix}:ref_gene_{gene_match.group(1)}",
                        "md5": md5_match.group(1).lower(),
                    }
                )
    return pd.DataFrame(rows).drop_duplicates()


def read_annotation(path: Path, source: str) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t")
    df = df.rename(columns={df.columns[0]: "locus_tag"})
    df["md5"] = df["contig"].astype(str).str.extract(r"md5=([0-9a-fA-F]+)", expand=False).str.lower()
    df["source"] = source
    return df.dropna(subset=["md5"])


def collapse_annotations(df: pd.DataFrame) -> pd.DataFrame:
    if df.empty:
        return pd.DataFrame(columns=["md5", "annotation_sources", "annotation_categories", "annotation_subcategories", "annotation_genes", "annotation_functions"])

    def uniq_join(values: pd.Series, limit: int = 8) -> str:
        vals = [str(v) for v in values.dropna().unique() if str(v) and str(v) != "-"]
        return "; ".join(vals[:limit])

    grouped = (
        df.groupby("md5", as_index=False)
        .agg(
            annotation_sources=("source", uniq_join),
            annotation_categories=("category", uniq_join),
            annotation_subcategories=("subcategory", uniq_join),
            annotation_genes=("gene", uniq_join),
            annotation_functions=("function", uniq_join),
        )
    )
    return grouped


def make_cluster_annotation_table(coords_path: Path, cluster: int, md5_map: pd.DataFrame, annotations: pd.DataFrame) -> pd.DataFrame:
    coords = pd.read_csv(coords_path, sep="\t")
    cluster_df = coords.loc[coords["kmeans_cluster"] == cluster].copy()
    cluster_df = cluster_df.sort_values("gene", key=lambda s: s.map(gene_number))
    table = cluster_df.merge(md5_map, on="gene", how="left").merge(annotations, on="md5", how="left")
    fill_cols = ["annotation_sources", "annotation_categories", "annotation_subcategories", "annotation_genes", "annotation_functions"]
    table[fill_cols] = table[fill_cols].fillna("")
    return table


def read_presence(profile: Path, genes: list[str]) -> pd.DataFrame:
    usecols = ["genome", *genes]
    df = pd.read_csv(profile, sep="\t", usecols=lambda col: col in set(usecols))
    present = df.set_index("genome").apply(pd.to_numeric, errors="coerce").fillna(0)
    present.index = [clean_tip_name(idx) for idx in present.index]
    return (present[genes] != 0).astype(int)


def plot_tree_heatmap(tree_path: Path, presence: pd.DataFrame, title: str, out_base: Path) -> pd.DataFrame:
    tree = parse_newick(tree_path.read_text())
    tree_tips = []

    def collect(node: Node) -> None:
        if node.is_leaf:
            tree_tips.append(clean_tip_name(node.name))
        for child in node.children:
            collect(child)

    collect(tree)
    keep = set(presence.index).intersection(tree_tips)
    pruned = prune_tree(tree, keep)
    if pruned is None:
        raise ValueError("No overlap between tree tips and profile genomes")

    leaf_order = [tip for tip in tree_tips if tip in keep]
    assign_tree_layout(pruned, leaf_order)
    ordered_presence = presence.loc[leaf_order]

    sns.set_theme(style="white", context="paper")
    height = max(9, len(leaf_order) * 0.010)
    width = max(12, len(ordered_presence.columns) * 0.13 + 4.5)
    fig = plt.figure(figsize=(width, height))
    grid = fig.add_gridspec(1, 2, width_ratios=[1.4, 4.8], wspace=0.015)
    tree_ax = fig.add_subplot(grid[0, 0])
    heat_ax = fig.add_subplot(grid[0, 1], sharey=tree_ax)

    draw_tree(tree_ax, pruned)
    tree_ax.set_ylim(len(leaf_order) - 0.5, -0.5)
    tree_ax.set_axis_off()

    cmap = ListedColormap(["#f4f4f4", "#1f77b4"])
    heat_ax.imshow(ordered_presence.to_numpy(), aspect="auto", interpolation="nearest", cmap=cmap, vmin=0, vmax=1)
    heat_ax.set_xticks(np.arange(len(ordered_presence.columns)))
    heat_ax.set_xticklabels([g.split(":")[-1] for g in ordered_presence.columns], rotation=90, fontsize=6)
    heat_ax.set_yticks([])
    heat_ax.set_title(title, fontsize=12)
    heat_ax.set_xlabel("cluster genes")
    for spine in heat_ax.spines.values():
        spine.set_visible(False)
    fig.tight_layout()
    fig.savefig(out_base.with_suffix(".png"), dpi=240)
    fig.savefig(out_base.with_suffix(".pdf"))
    plt.close(fig)

    return pd.DataFrame({"genome": leaf_order})


def summarize_table(table: pd.DataFrame, label: str) -> pd.DataFrame:
    rows = [{"dataset": label, "metric": "n_genes", "value": len(table)}]
    annotated = table["annotation_sources"].astype(bool)
    rows.append({"dataset": label, "metric": "n_annotated_any", "value": int(annotated.sum())})
    rows.append({"dataset": label, "metric": "n_unannotated", "value": int((~annotated).sum())})
    for source in ["Cazyme", "VIRULENCE", "Mobile_Element", "Transporter"]:
        rows.append(
            {
                "dataset": label,
                "metric": f"n_{source}",
                "value": int(table["annotation_categories"].str.contains(source, regex=False, na=False).sum()),
            }
        )
    return pd.DataFrame(rows)


def main() -> None:
    parser = argparse.ArgumentParser(description="Annotate selected PCA clusters and draw Phandango-like tree heatmaps.")
    parser.add_argument("--outdir", type=Path, default=Path("Leca/cluster3_annotations_0427_1"))
    parser.add_argument("--order-cluster", type=int, default=3)
    parser.add_argument("--bp-cluster", type=int, default=2)
    parser.add_argument("--gene-prefix", default="0427_1")
    parser.add_argument("--gff", type=Path, default=Path("Leca/klebsiella_reference.gff"))
    parser.add_argument("--profile", type=Path, default=Path("Leca/0427_1.profile"))
    parser.add_argument("--tree", type=Path, default=Path("Leca/Kp.labelled.nwk"))
    parser.add_argument("--order-coords", type=Path, default=Path("Leca/order_bp_phylo_compare_0427_1/order_phylo_pca_coords.tsv"))
    parser.add_argument("--bp-coords", type=Path, default=Path("Leca/order_bp_phylo_compare_0427_1/bp_phylo_pca_coords.tsv"))
    parser.add_argument("--cazyme", type=Path, default=Path("Leca/klebsiella_reference.cazyme.txt"))
    parser.add_argument("--effectors", type=Path, default=Path("Leca/klebsiella_reference.effectors.txt"))
    parser.add_argument("--mobile-elements", type=Path, default=Path("Leca/klebsiella_reference.mobile_elements.txt"))
    parser.add_argument("--transporters", type=Path, default=Path("Leca/klebsiella_reference.transporters.txt"))
    args = parser.parse_args()

    args.outdir.mkdir(parents=True, exist_ok=True)
    md5_map = read_ref_gene_md5(args.gff, args.gene_prefix)
    ann = pd.concat(
        [
            read_annotation(args.cazyme, "cazyme"),
            read_annotation(args.effectors, "effectors"),
            read_annotation(args.mobile_elements, "mobile_elements"),
            read_annotation(args.transporters, "transporters"),
        ],
        ignore_index=True,
    )
    annotations = collapse_annotations(ann)

    summary_tables = []
    for label, coords_path, cluster in [
        ("order_phylo_pca", args.order_coords, args.order_cluster),
        ("bp_phylo_pca", args.bp_coords, args.bp_cluster),
    ]:
        table = make_cluster_annotation_table(coords_path, cluster, md5_map, annotations)
        table.to_csv(args.outdir / f"{label}_cluster{cluster}_gene_annotations.tsv", sep="\t", index=False)
        summary_tables.append(summarize_table(table, label))

        presence = read_presence(args.profile, table["gene"].tolist())
        prevalence = presence.mean(axis=0).rename("presence_fraction").reset_index().rename(columns={"index": "gene"})
        table.merge(prevalence, on="gene", how="left").to_csv(
            args.outdir / f"{label}_cluster{cluster}_gene_annotations_with_prevalence.tsv",
            sep="\t",
            index=False,
        )
        ordered = plot_tree_heatmap(
            args.tree,
            presence,
            f"{label} cluster {cluster} gene presence",
            args.outdir / f"{label}_cluster{cluster}_tree_presence_heatmap",
        )
        ordered.to_csv(args.outdir / f"{label}_cluster{cluster}_tree_order.tsv", sep="\t", index=False)

    summary = pd.concat(summary_tables, ignore_index=True)
    summary.to_csv(args.outdir / "selected_cluster_annotation_summary.tsv", sep="\t", index=False)
    print(summary.to_string(index=False))


if __name__ == "__main__":
    main()
