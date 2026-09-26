#!/usr/bin/env python3
"""Draw tree-aligned presence/absence heatmaps for real SHC seeds and toy triples."""

from pathlib import Path
from dataclasses import dataclass, field
import re

import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[4]
WF = ROOT / "Leca/workflows/shc_diff_gwes"


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


def clean_tip_name(name: str) -> str:
    return str(name).removesuffix(".fna")


def parse_newick(text: str) -> Node:
    text = text.strip().rstrip(";")
    index = 0

    def label_and_length() -> tuple[str, float]:
        nonlocal index
        start = index
        while index < len(text) and text[index] not in ",()":
            index += 1
        token = text[start:index].strip()
        if not token:
            return "", 0.0
        if ":" not in token:
            return token, 0.0
        name, length = token.rsplit(":", 1)
        try:
            return name, float(length)
        except ValueError:
            return name, 0.0

    def subtree() -> Node:
        nonlocal index
        if text[index] == "(":
            index += 1
            children = []
            while True:
                children.append(subtree())
                if index >= len(text):
                    break
                if text[index] == ",":
                    index += 1
                    continue
                if text[index] == ")":
                    index += 1
                    break
            name, length = label_and_length()
            return Node(name=name, length=length, children=children)
        name, length = label_and_length()
        return Node(name=name, length=length)

    return subtree()


def prune_tree(node: Node, keep: set[str]) -> Node | None:
    if node.is_leaf:
        name = clean_tip_name(node.name)
        return Node(name=name, length=node.length) if name in keep else None
    children = [child for child in (prune_tree(child, keep) for child in node.children) if child]
    return Node(name=node.name, length=node.length, children=children) if children else None


def assign_layout(node: Node, y_lookup: dict[str, int], parent_x: float = 0.0) -> None:
    node.x = parent_x + max(node.length, 0.0)
    if node.is_leaf:
        node.y = float(y_lookup[clean_tip_name(node.name)])
        return
    for child in node.children:
        assign_layout(child, y_lookup, node.x)
    node.y = float(np.mean([child.y for child in node.children]))


def draw_tree(axis: plt.Axes, node: Node) -> None:
    if node.is_leaf:
        return
    child_y = [child.y for child in node.children]
    axis.plot([node.x, node.x], [min(child_y), max(child_y)], color="#252A30", lw=0.35)
    for child in node.children:
        axis.plot([node.x, child.x], [child.y, child.y], color="#252A30", lw=0.35)
        draw_tree(axis, child)


def plot_tree_heatmap(
    tree_path: Path,
    presence: pd.DataFrame,
    title: str,
    out_base: Path,
    present_color: str,
    column_groups: list[int] | None = None,
) -> pd.DataFrame:
    tree = parse_newick(tree_path.read_text())
    tips = []

    def collect(node: Node) -> None:
        if node.is_leaf:
            tips.append(clean_tip_name(node.name))
        for child in node.children:
            collect(child)

    collect(tree)
    leaf_order = [tip for tip in tips if tip in presence.index]
    pruned = prune_tree(tree, set(leaf_order))
    if pruned is None:
        raise ValueError("No overlap between tree tips and profile genomes")
    assign_layout(pruned, {tip: index for index, tip in enumerate(leaf_order)})
    ordered = presence.loc[leaf_order]

    height = max(9.0, len(leaf_order) * 0.010)
    width = max(12.0, len(ordered.columns) * 0.13 + 4.8)
    fig = plt.figure(figsize=(width, height), facecolor="white")
    grid = fig.add_gridspec(1, 2, width_ratios=[1.35, 4.8], wspace=0.012)
    tree_axis = fig.add_subplot(grid[0, 0])
    heat_axis = fig.add_subplot(grid[0, 1], sharey=tree_axis)
    draw_tree(tree_axis, pruned)
    tree_axis.set_ylim(len(leaf_order) - 0.5, -0.5)
    tree_axis.set_axis_off()
    heat_axis.imshow(
        ordered.to_numpy(), aspect="auto", interpolation="nearest",
        cmap=ListedColormap(["#F1F3F4", present_color]), vmin=0, vmax=1,
    )
    heat_axis.set_xticks(np.arange(len(ordered.columns)))
    heat_axis.set_xticklabels(ordered.columns, rotation=90, fontsize=6)
    heat_axis.set_yticks([])
    heat_axis.set_title(title, fontsize=12, loc="left")
    heat_axis.set_xlabel("Locus (light grey = absent; color = present)")
    if column_groups is not None:
        if len(column_groups) != len(ordered.columns):
            raise ValueError("column_groups must have one value per heatmap column")
        starts = [0]
        for index in range(1, len(column_groups)):
            if column_groups[index] != column_groups[index - 1]:
                heat_axis.axvline(index - 0.5, color="#FFFFFF", linewidth=2.2)
                starts.append(index)
        ends = starts[1:] + [len(column_groups)]
        top_axis = heat_axis.secondary_xaxis("top")
        top_axis.set_xticks([(start + end - 1) / 2 for start, end in zip(starts, ends)])
        top_axis.set_xticklabels([f"Community {column_groups[start]}" for start in starts], fontsize=8)
        top_axis.tick_params(length=0, pad=4)
        top_axis.spines["top"].set_visible(False)
    for spine in heat_axis.spines.values():
        spine.set_visible(False)
    fig.subplots_adjust(left=0.02, right=0.995, top=0.95, bottom=0.18)
    fig.savefig(out_base.with_suffix(".png"), dpi=240)
    fig.savefig(out_base.with_suffix(".pdf"))
    plt.close(fig)
    return pd.DataFrame({"genome": leaf_order})


def read_annotations(path: Path, genes: set[str]) -> dict[str, str]:
    annotations = {}
    pattern = re.compile(r"ref_gene_(\d+)")
    with path.open() as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 9 or fields[2] != "CDS":
                continue
            match = pattern.search(fields[0])
            if not match:
                continue
            gene = f"0427_1:ref_gene_{match.group(1)}"
            if gene not in genes:
                continue
            attributes = dict(item.split("=", 1) for item in fields[8].split(";") if "=" in item)
            annotations[gene] = attributes.get("gene", attributes.get("Name", gene.split(":")[-1]))
    return annotations


def read_presence(profile: Path, genes: list[str]) -> pd.DataFrame:
    wanted = {"genome", *genes}
    table = pd.read_csv(profile, sep="\t", usecols=lambda column: column in wanted)
    missing = sorted(set(genes) - set(table.columns))
    if missing:
        raise ValueError(f"Profile is missing {len(missing)} requested loci: {missing[:5]}")
    values = table.set_index("genome").apply(pd.to_numeric, errors="coerce").fillna(0)
    values.index = [clean_tip_name(index) for index in values.index]
    return (values[genes] > 0).astype(int)


def real_locus_order(seeds: pd.DataFrame) -> list[str]:
    def number(gene: str) -> int:
        match = re.search(r"ref_gene_(\d+)$", gene)
        return int(match.group(1)) if match else 10**12

    return sorted(set(seeds.gene1) | set(seeds.gene2), key=number)


def main() -> None:
    figure_dir = WF / "output/figures"
    table_dir = WF / "output/phandango"
    figure_dir.mkdir(parents=True, exist_ok=True)
    table_dir.mkdir(parents=True, exist_ok=True)
    tree = ROOT / "Leca/source/tree/Kp.labelled.nwk"

    seeds = pd.read_csv(WF / "input/seed_pairs.tsv", sep="\t")
    real_genes = real_locus_order(seeds)
    annotations = read_annotations(
        ROOT / "Leca/source/annotation/klebsiella_reference.gff", set(real_genes)
    )
    labels = {
        gene: (
            f"{annotations[gene]} | {gene.split(':')[-1]}"
            if annotations.get(gene, "") not in ("", gene.split(":")[-1])
            else gene.split(":")[-1]
        )
        for gene in real_genes
    }
    real_presence = read_presence(ROOT / "Leca/source/profile/0427_1.profile", real_genes)
    real_order = plot_tree_heatmap(
        tree,
        real_presence.rename(columns=labels),
        "True Kp: 12 loci from 25 SHC-strict seed pairs",
        figure_dir / "true_shc_seed_tree_presence_heatmap",
        "#3976A8",
    )
    real_presence.loc[real_order.genome].to_csv(
        table_dir / "true_shc_seed_phandango_presence.tsv", sep="\t", index_label="genome"
    )
    real_order.to_csv(table_dir / "true_shc_seed_tree_tip_order.tsv", sep="\t", index=False)
    pd.DataFrame(
        {
            "gene": real_genes,
            "display_label": [labels[gene] for gene in real_genes],
            "annotation": [annotations.get(gene, "") for gene in real_genes],
            "presence_count": [int(real_presence[gene].sum()) for gene in real_genes],
            "presence_fraction": [float(real_presence[gene].mean()) for gene in real_genes],
        }
    ).to_csv(table_dir / "true_shc_seed_locus_summary.tsv", sep="\t", index=False)

    truth = pd.read_csv(WF / "input/toy_truth.tsv", sep="\t")
    toy_genes = [gene for row in truth.itertuples(index=False) for gene in (row.gene1, row.gene2, row.background)]
    toy_labels = {
        gene: gene.replace("toy_epi_", "E").replace("_", "-") for gene in toy_genes
    }
    toy_presence = read_presence(WF / "input/toy_profile.tsv.gz", toy_genes)
    toy_order = plot_tree_heatmap(
        tree,
        toy_presence.rename(columns=toy_labels),
        "Toy Kp: 20 significant injected A-B-C triples",
        figure_dir / "toy_significant_triples_tree_presence_heatmap",
        "#D46A3A",
    )
    toy_presence.loc[toy_order.genome].to_csv(
        table_dir / "toy_significant_triples_phandango_presence.tsv",
        sep="\t",
        index_label="genome",
    )
    toy_order.to_csv(
        table_dir / "toy_significant_triples_tree_tip_order.tsv", sep="\t", index=False
    )
    pd.DataFrame(
        {
            "gene": toy_genes,
            "display_label": [toy_labels[gene] for gene in toy_genes],
            "triple": [index for index in range(1, 21) for _ in range(3)],
            "role": [role for _ in range(20) for role in ("A", "B", "C")],
            "presence_count": [int(toy_presence[gene].sum()) for gene in toy_genes],
            "presence_fraction": [float(toy_presence[gene].mean()) for gene in toy_genes],
        }
    ).to_csv(table_dir / "toy_significant_triple_locus_summary.tsv", sep="\t", index=False)

    pd.DataFrame(
        [
            {"dataset": "true", "tree_tips": len(real_order), "loci": len(real_genes), "pairs_or_triples": len(seeds)},
            {"dataset": "toy", "tree_tips": len(toy_order), "loci": len(toy_genes), "pairs_or_triples": len(truth)},
        ]
    ).to_csv(table_dir / "summary.tsv", sep="\t", index=False)


if __name__ == "__main__":
    main()
