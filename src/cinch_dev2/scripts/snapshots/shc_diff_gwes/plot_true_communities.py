#!/usr/bin/env python3
"""Community views of the true SHC-strict seed-pair network."""

from itertools import combinations
from pathlib import Path
import re

import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.lines import Line2D
import numpy as np
import pandas as pd

from plot_phandango_presence import plot_tree_heatmap, read_annotations, read_presence


ROOT = Path(__file__).resolve().parents[4]
WF = ROOT / "Leca/workflows/shc_diff_gwes"
COLORS = ["#3976A8", "#D46A3A", "#3B8C6E", "#8A65A5", "#C49A2C", "#64727D"]
INK, MUTED, GRID = "#20252C", "#68717D", "#E2E7EB"


def gene_number(gene: str) -> int:
    match = re.search(r"ref_gene_(\d+)$", gene)
    return int(match.group(1)) if match else 10**12


def modularity(communities: list[set[str]], weights: dict[tuple[str, str], float]) -> float:
    nodes = sorted({node for edge in weights for node in edge})
    degrees = {node: sum(value for (left, _), value in weights.items() if left == node) for node in nodes}
    total_edge_weight = sum(degrees.values()) / 2.0
    if total_edge_weight <= 0:
        return 0.0
    score = 0.0
    for community in communities:
        internal = sum(weights.get((left, right), 0.0) for left, right in combinations(community, 2))
        degree_sum = sum(degrees[node] for node in community)
        score += internal / total_edge_weight - (degree_sum / (2.0 * total_edge_weight)) ** 2
    return float(score)


def greedy_modularity(nodes: list[str], weights: dict[tuple[str, str], float]) -> tuple[list[set[str]], float]:
    communities = [{node} for node in nodes]
    while True:
        current = modularity(communities, weights)
        best_delta = 0.0
        best = None
        for left, right in combinations(range(len(communities)), 2):
            candidate = [community for index, community in enumerate(communities) if index not in (left, right)]
            candidate.append(communities[left] | communities[right])
            delta = modularity(candidate, weights) - current
            tie_key = tuple(sorted(communities[left] | communities[right], key=gene_number))
            if delta > best_delta + 1e-15 or (
                abs(delta - best_delta) <= 1e-15 and best is not None and tie_key < best[0]
            ):
                best_delta = delta
                best = (tie_key, candidate)
        if best is None or best_delta <= 1e-15:
            break
        communities = best[1]
    communities = sorted(communities, key=lambda group: min(gene_number(node) for node in group))
    return communities, modularity(communities, weights)


def community_layout(communities: list[set[str]]) -> dict[str, tuple[float, float]]:
    centers = [(-1.5, 1.0), (1.5, 1.0), (-1.5, -1.0), (1.5, -1.0), (0.0, 0.0), (0.0, -2.0)]
    positions = {}
    for community_index, community in enumerate(communities):
        ordered = sorted(community, key=gene_number)
        center_x, center_y = centers[community_index]
        if len(ordered) == 1:
            positions[ordered[0]] = (center_x, center_y)
            continue
        angles = np.linspace(0, 2 * np.pi, len(ordered), endpoint=False) + np.pi / 4
        radius = 0.50 if len(ordered) <= 3 else 0.62
        for node, angle in zip(ordered, angles):
            positions[node] = (center_x + radius * np.cos(angle), center_y + radius * np.sin(angle))
    return positions


def main() -> None:
    figure_dir = WF / "output/figures"
    table_dir = WF / "output/true"
    phandango_dir = WF / "output/phandango"
    for directory in (figure_dir, table_dir, phandango_dir):
        directory.mkdir(parents=True, exist_ok=True)

    edges = pd.read_csv(WF / "input/seed_pairs.tsv", sep="\t")
    nodes = sorted(set(edges.gene1) | set(edges.gene2), key=gene_number)
    weights = {}
    for row in edges.itertuples(index=False):
        weights[(row.gene1, row.gene2)] = float(row.conditional_mi)
        weights[(row.gene2, row.gene1)] = float(row.conditional_mi)
    communities, score = greedy_modularity(nodes, weights)
    assignment = {node: index + 1 for index, community in enumerate(communities) for node in community}
    weighted_degree = {node: sum(weights.get((node, other), 0.0) for other in nodes) for node in nodes}
    ordered_nodes = [
        node
        for community in communities
        for node in sorted(community, key=lambda item: (-weighted_degree[item], gene_number(item)))
    ]

    annotations = read_annotations(ROOT / "Leca/source/annotation/klebsiella_reference.gff", set(nodes))
    assignments = pd.DataFrame(
        {
            "gene": ordered_nodes,
            "community": [assignment[node] for node in ordered_nodes],
            "weighted_degree": [weighted_degree[node] for node in ordered_nodes],
            "annotation": [annotations.get(node, "") for node in ordered_nodes],
        }
    )
    assignments.to_csv(table_dir / "true_shc_community_assignments.tsv", sep="\t", index=False)
    edge_output = edges.copy()
    edge_output["gene1_community"] = edge_output.gene1.map(assignment)
    edge_output["gene2_community"] = edge_output.gene2.map(assignment)
    edge_output["within_community"] = edge_output.gene1_community == edge_output.gene2_community
    edge_output.to_csv(table_dir / "true_shc_community_edges.tsv", sep="\t", index=False)
    pd.DataFrame(
        [
            {
                "community": index + 1,
                "n_loci": len(community),
                "n_internal_edges": int(
                    sum(1 for row in edge_output.itertuples(index=False) if row.within_community and row.gene1 in community)
                ),
                "loci": ";".join(sorted(community, key=gene_number)),
                "weighted_modularity_total": score,
            }
            for index, community in enumerate(communities)
        ]
    ).to_csv(table_dir / "true_shc_community_summary.tsv", sep="\t", index=False)

    positions = community_layout(communities)
    fig, axes = plt.subplots(1, 2, figsize=(14, 7), facecolor="white", gridspec_kw={"width_ratios": [1.05, 1.0]})
    fig.subplots_adjust(left=0.055, right=0.98, bottom=0.13, top=0.82, wspace=0.20)
    fig.text(0.055, 0.965, "True Kp SHC-strict association communities", ha="left", va="top", fontsize=18, weight="bold", color=INK)
    fig.text(0.055, 0.915, f"Weighted greedy-modularity partition of 12 loci and 25 conditional-MI edges; four communities, Q={score:.3f}. Communities summarize network topology and are not additional significance tests.", ha="left", va="top", fontsize=9, color=MUTED)

    network_axis = axes[0]
    maximum_weight = edges.conditional_mi.max()
    for row in edges.itertuples(index=False):
        x1, y1 = positions[row.gene1]
        x2, y2 = positions[row.gene2]
        within = assignment[row.gene1] == assignment[row.gene2]
        network_axis.plot(
            [x1, x2], [y1, y2],
            color=COLORS[assignment[row.gene1] - 1] if within else "#AAB2BA",
            linewidth=0.6 + 3.4 * row.conditional_mi / maximum_weight,
            alpha=0.72 if within else 0.42,
            linestyle="-" if "order" in row.source else "--",
            zorder=1,
        )
    for node in nodes:
        x, y = positions[node]
        community = assignment[node]
        network_axis.scatter(x, y, s=260, color=COLORS[community - 1], edgecolors="white", linewidths=1.2, zorder=3)
        network_axis.text(x, y - 0.16, node.split(":")[-1], ha="center", va="top", fontsize=7.5, color=INK, zorder=4)
    network_axis.set_title("A  Conditional-MI network", loc="left", fontsize=12, weight="bold")
    network_axis.set_xlim(-2.4, 2.4); network_axis.set_ylim(-1.85, 1.85); network_axis.set_axis_off()
    network_axis.legend(
        handles=[
            Line2D([0], [0], color="#68717D", lw=2.2, linestyle="-", label="Order + BP selected"),
            Line2D([0], [0], color="#68717D", lw=2.2, linestyle="--", label="BP-only selected"),
            Line2D([0], [0], color="#AAB2BA", lw=2.2, label="Between communities"),
        ],
        frameon=False,
        fontsize=7,
        loc="lower center",
        ncol=3,
    )

    matrix = np.full((len(nodes), len(nodes)), np.nan)
    for index in range(len(nodes)):
        matrix[index, index] = 0.0
    lookup = {node: index for index, node in enumerate(ordered_nodes)}
    for row in edges.itertuples(index=False):
        left, right = lookup[row.gene1], lookup[row.gene2]
        matrix[left, right] = matrix[right, left] = row.conditional_mi
    matrix_axis = axes[1]
    cmap = LinearSegmentedColormap.from_list("conditional_mi", ["#F4F5F6", "#9CC8BC", "#176B58"])
    masked = np.ma.masked_invalid(matrix)
    image = matrix_axis.imshow(masked, cmap=cmap, vmin=0, vmax=maximum_weight, interpolation="nearest")
    matrix_axis.set_xticks(np.arange(len(nodes)), [node.split(":")[-1] for node in ordered_nodes], rotation=90, fontsize=7)
    matrix_axis.set_yticks(np.arange(len(nodes)), [node.split(":")[-1] for node in ordered_nodes], fontsize=7)
    boundary = np.cumsum([len(community) for community in communities])[:-1] - 0.5
    for value in boundary:
        matrix_axis.axhline(value, color="white", linewidth=2.2)
        matrix_axis.axvline(value, color="white", linewidth=2.2)
    matrix_axis.set_title("B  Community-sorted edge matrix", loc="left", fontsize=12, weight="bold")
    colorbar = fig.colorbar(image, ax=matrix_axis, fraction=0.046, pad=0.04)
    colorbar_ticks = np.linspace(0, maximum_weight, 5)
    colorbar.set_ticks(colorbar_ticks)
    colorbar.set_ticklabels([f"{value:.3f}" for value in colorbar_ticks])
    colorbar.set_label("Conditional MI (nats)", fontsize=8)
    for spine in matrix_axis.spines.values():
        spine.set_visible(False)
    fig.savefig(figure_dir / "true_shc_community_network.png", dpi=240)
    fig.savefig(figure_dir / "true_shc_community_network.svg")
    plt.close(fig)

    presence = read_presence(ROOT / "Leca/source/profile/0427_1.profile", ordered_nodes)
    display = {
        node: f"C{assignment[node]} | {node.split(':')[-1]}" for node in ordered_nodes
    }
    tree_order = plot_tree_heatmap(
        ROOT / "Leca/source/tree/Kp.labelled.nwk",
        presence.rename(columns=display),
        "True Kp: SHC-strict loci ordered by weighted network community",
        figure_dir / "true_shc_community_tree_presence_heatmap",
        "#3976A8",
        [assignment[node] for node in ordered_nodes],
    )
    presence.loc[tree_order.genome, ordered_nodes].to_csv(
        phandango_dir / "true_shc_community_phandango_presence.tsv",
        sep="\t",
        index_label="genome",
    )
    tree_order.to_csv(
        phandango_dir / "true_shc_community_tree_tip_order.tsv", sep="\t", index=False
    )


if __name__ == "__main__":
    main()
