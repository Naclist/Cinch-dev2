#!/usr/bin/env python3
import argparse
import re
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from numba import njit, prange
from scipy.stats import chi2


def canonical_genome(x: str) -> str:
    match = re.search(r"(GC[AF]_\d+)", str(x))
    if match:
        return match.group(1)
    s = str(x).strip()
    for suf in (".fna", ".fa", ".fasta"):
        if s.endswith(suf):
            s = s[: -len(suf)]
    return s


UNKNOWN_VALUES = {
    "",
    "unknown",
    "not_available",
    "not available",
    "missing",
    "na",
    "nan",
    "none",
    "null",
    "unspecified",
}


def ref_order(gene: str) -> int:
    m = re.search(r"ref_gene_(\d+)$", str(gene))
    return int(m.group(1)) if m else -1


class NewickNode:
    def __init__(self) -> None:
        self.name = ""
        self.length = 0.0
        self.children: list["NewickNode"] = []
        self.parent: "NewickNode | None" = None


def parse_newick(text: str) -> NewickNode:
    text = text.strip().rstrip(";")

    def parse_name_length(node: NewickNode, pos: int) -> int:
        start = pos
        while pos < len(text) and text[pos] not in ",():":
            pos += 1
        node.name = text[start:pos].strip()
        if pos < len(text) and text[pos] == ":":
            pos += 1
            start = pos
            while pos < len(text) and text[pos] not in ",()":
                pos += 1
            try:
                node.length = float(text[start:pos])
            except ValueError:
                node.length = 0.0
        return pos

    def parse_node(pos: int) -> tuple[NewickNode, int]:
        node = NewickNode()
        if text[pos] == "(":
            pos += 1
            while True:
                child, pos = parse_node(pos)
                child.parent = node
                node.children.append(child)
                if pos >= len(text):
                    break
                if text[pos] == ",":
                    pos += 1
                    continue
                if text[pos] == ")":
                    pos += 1
                    break
        pos = parse_name_length(node, pos)
        return node, pos

    root, _ = parse_node(0)
    return root


def tree_distance_matrix(tree_path: Path, genome_order: pd.Index) -> np.ndarray:
    root = parse_newick(tree_path.read_text())
    tips: dict[str, NewickNode] = {}
    depths: dict[NewickNode, float] = {}
    ancestors: dict[NewickNode, list[NewickNode]] = {}

    def walk(node: NewickNode, depth: float, path: list[NewickNode]) -> None:
        depths[node] = depth
        ancestors[node] = path + [node]
        if not node.children and node.name:
            tips[canonical_genome(node.name)] = node
        for child in node.children:
            walk(child, depth + child.length, ancestors[node])

    walk(root, 0.0, [])
    nodes = []
    missing = []
    for genome in genome_order:
        key = canonical_genome(genome)
        if key in tips:
            nodes.append(tips[key])
        else:
            missing.append(str(genome))
    if missing:
        raise ValueError(f"Tree is missing {len(missing)} profile genomes, first missing: {missing[:5]}")

    n = len(nodes)
    dist = np.zeros((n, n), dtype=np.float32)
    ancestor_sets = [set(ancestors[node]) for node in nodes]
    for i in range(n - 1):
        ai = ancestors[nodes[i]]
        for j in range(i + 1, n):
            common = ancestor_sets[j]
            lca = root
            for a in ai:
                if a in common:
                    lca = a
            d = depths[nodes[i]] + depths[nodes[j]] - 2.0 * depths[lca]
            dist[i, j] = dist[j, i] = d
    return dist


def classical_pcoa(distance_matrix: np.ndarray, n_components: int) -> np.ndarray:
    n = distance_matrix.shape[0]
    d2 = np.square(distance_matrix.astype(np.float64))
    j = np.eye(n) - np.ones((n, n)) / n
    b = -0.5 * j @ d2 @ j
    eigvals, eigvecs = np.linalg.eigh(b)
    order = np.argsort(eigvals)[::-1]
    eigvals = eigvals[order]
    eigvecs = eigvecs[:, order]
    return eigvecs[:, :n_components] * np.sqrt(np.maximum(eigvals[:n_components], 0))


def tree_pc_strata(
    tree_path: Path,
    genome_order: pd.Index,
    n_pcs: int,
    bins: int,
    min_stratum_size: int,
    outdir: Path,
) -> np.ndarray:
    dist = tree_distance_matrix(tree_path, genome_order)
    pcs = classical_pcoa(dist, n_pcs)
    pc_df = pd.DataFrame(pcs, columns=[f"treePC{i + 1}" for i in range(n_pcs)])
    pc_df.insert(0, "genome", genome_order.to_list())
    pc_df.to_csv(outdir / "tree_pc_coordinates.tsv", sep="\t", index=False)
    labels = pd.Series([""] * len(genome_order), dtype=object)
    for i in range(n_pcs):
        q = pd.qcut(pcs[:, i], q=bins, labels=False, duplicates="drop")
        labels = labels + ("|" if i else "") + pd.Series(q).fillna(-1).astype(int).astype(str)
    counts = labels.value_counts()
    labels = labels.where(labels.map(counts) >= min_stratum_size, "SMALL_TREE_PC_STRATA")
    counts = labels.value_counts()
    if "SMALL_TREE_PC_STRATA" in counts and counts["SMALL_TREE_PC_STRATA"] < min_stratum_size and len(counts) > 1:
        largest = counts.drop(index="SMALL_TREE_PC_STRATA").idxmax()
        labels = labels.where(labels != "SMALL_TREE_PC_STRATA", largest)
    return pd.factorize(labels)[0].astype(np.int32)


def clean_metadata(path: Path) -> pd.DataFrame:
    meta = pd.read_csv(path, sep="\t", dtype=str, on_bad_lines="skip", index_col=False)
    meta["genome_key"] = meta["Query"].map(canonical_genome)
    meta = meta[meta["genome_key"].str.match(r"GC[AF]_\d+", na=False)].copy()
    meta = meta.drop_duplicates("genome_key", keep="first")
    for col in meta.columns:
        if col == "genome_key":
            continue
        vals = meta[col].astype(str).str.strip()
        meta[col] = vals.mask(vals.str.lower().isin(UNKNOWN_VALUES), np.nan)
    return meta


def build_strata(profile_index: pd.Index, metadata: Path | None, covariates: list[str], min_stratum_size: int) -> np.ndarray:
    if metadata is None or not covariates:
        return np.zeros(len(profile_index), dtype=np.int32)
    meta = clean_metadata(metadata)
    keys = pd.DataFrame({"genome_key": [canonical_genome(x) for x in profile_index]})
    merged = keys.merge(meta, on="genome_key", how="left")
    available = [c for c in covariates if c in merged.columns]
    if not available:
        return np.zeros(len(profile_index), dtype=np.int32)
    parts = []
    for c in available:
        v = merged[c].fillna("Unknown").astype(str).str.strip()
        v = v.mask(v.eq("") | v.str.lower().isin(["nan", "unknown", "not_available"]), "Unknown")
        parts.append(v)
    labels = parts[0]
    for p in parts[1:]:
        labels = labels + "|" + p
    counts = labels.value_counts()
    labels = labels.where(labels.map(counts) >= min_stratum_size, "SMALL_OR_UNKNOWN")
    return pd.factorize(labels)[0].astype(np.int32)


def encode_core_alleles(profile: pd.DataFrame, core_min: float, allele_min_count: int) -> tuple[np.ndarray, list[str], np.ndarray, pd.DataFrame]:
    X = profile.to_numpy()
    present = X > 0
    prev = present.mean(axis=0)
    keep = prev >= core_min
    genes = profile.columns[keep].tolist()
    X = X[:, keep]
    encoded = np.zeros_like(X, dtype=np.int16)
    k_eff = np.zeros(X.shape[1], dtype=np.int16)
    rows = []
    for j, g in enumerate(genes):
        vals = X[:, j]
        pos = vals[vals > 0]
        alleles, counts = np.unique(pos, return_counts=True)
        keep_alleles = alleles[counts >= allele_min_count]
        mapping = {int(a): i + 1 for i, a in enumerate(keep_alleles)}
        out = np.zeros(vals.shape[0], dtype=np.int16)
        for a, code in mapping.items():
            out[vals == a] = code
        encoded[:, j] = out
        k_eff[j] = len(keep_alleles) + 1
        rows.append(
            {
                "gene": g,
                "prevalence": float((vals > 0).mean()),
                "raw_n_alleles": int(len(alleles)),
                "kept_alleles": int(len(keep_alleles)),
                "encoded_categories_including_rare_or_missing": int(k_eff[j]),
                "singletons": int((counts == 1).sum()),
            }
        )
    return encoded, genes, k_eff, pd.DataFrame(rows)


@njit
def entropy_from_counts(counts, n):
    h = 0.0
    for c in counts:
        if c > 0:
            p = c / n
            h -= p * np.log(p)
    return h


@njit
def pair_nmi_gtest(x, y, kx, ky):
    n = x.shape[0]
    cx = np.zeros(kx, dtype=np.int32)
    cy = np.zeros(ky, dtype=np.int32)
    cxy = np.zeros((kx, ky), dtype=np.int32)
    for i in range(n):
        a = x[i]
        b = y[i]
        cx[a] += 1
        cy[b] += 1
        cxy[a, b] += 1
    hx = entropy_from_counts(cx, n)
    hy = entropy_from_counts(cy, n)
    hxy = 0.0
    for a in range(kx):
        for b in range(ky):
            c = cxy[a, b]
            if c > 0:
                p = c / n
                hxy -= p * np.log(p)
    mi = hx + hy - hxy
    g_stat = 2.0 * n * max(0.0, mi)
    used_x = 0
    used_y = 0
    for a in range(kx):
        if cx[a] > 0:
            used_x += 1
    for b in range(ky):
        if cy[b] > 0:
            used_y += 1
    df = max(1, (used_x - 1) * (used_y - 1))
    if hx + hy <= 0:
        return 0.0, g_stat, df
    return max(0.0, mi / ((hx + hy) / 2.0)), g_stat, df


@njit
def pair_nmi(x, y, kx, ky):
    nmi, _, _ = pair_nmi_gtest(x, y, kx, ky)
    return nmi


@njit
def pair_stratified_gtest(x, y, kx, ky, strata):
    max_s = 0
    for i in range(strata.shape[0]):
        if strata[i] > max_s:
            max_s = strata[i]
    total_g = 0.0
    total_df = 0
    for s in range(max_s + 1):
        n = 0
        for i in range(strata.shape[0]):
            if strata[i] == s:
                n += 1
        if n <= 1:
            continue
        xs = np.empty(n, dtype=np.int16)
        ys = np.empty(n, dtype=np.int16)
        pos = 0
        for i in range(strata.shape[0]):
            if strata[i] == s:
                xs[pos] = x[i]
                ys[pos] = y[i]
                pos += 1
        _, g, df = pair_nmi_gtest(xs, ys, kx, ky)
        total_g += g
        total_df += df
    if total_df < 1:
        total_df = 1
    return total_g, total_df


@njit(parallel=True)
def all_pair_nmi_gtest(encoded, k_eff, strata):
    n_genes = encoded.shape[1]
    n_pairs = n_genes * (n_genes - 1) // 2
    g1 = np.empty(n_pairs, dtype=np.int32)
    g2 = np.empty(n_pairs, dtype=np.int32)
    scores = np.empty(n_pairs, dtype=np.float32)
    g_stats = np.empty(n_pairs, dtype=np.float32)
    dfs = np.empty(n_pairs, dtype=np.int32)
    idx = 0
    for i in range(n_genes - 1):
        for j in range(i + 1, n_genes):
            g1[idx] = i
            g2[idx] = j
            idx += 1
    for p in prange(n_pairs):
        i = g1[p]
        j = g2[p]
        nmi, _, _ = pair_nmi_gtest(encoded[:, i], encoded[:, j], int(k_eff[i]), int(k_eff[j]))
        g, df = pair_stratified_gtest(encoded[:, i], encoded[:, j], int(k_eff[i]), int(k_eff[j]), strata)
        scores[p] = nmi
        g_stats[p] = g
        dfs[p] = df
    return g1, g2, scores, g_stats, dfs


@njit
def permute_within_strata(y, strata, seed):
    np.random.seed(seed)
    out = y.copy()
    max_s = 0
    for i in range(strata.shape[0]):
        if strata[i] > max_s:
            max_s = strata[i]
    for s in range(max_s + 1):
        idx = np.where(strata == s)[0]
        n = idx.shape[0]
        for i in range(n - 1, 0, -1):
            j = np.random.randint(0, i + 1)
            a = idx[i]
            b = idx[j]
            tmp = out[a]
            out[a] = out[b]
            out[b] = tmp
    return out


def empirical_for_candidates(encoded, k_eff, candidates, strata, n_perm: int, seed: int) -> pd.DataFrame:
    rng = np.random.default_rng(seed)
    rows = []
    for row in candidates.itertuples(index=False):
        i = int(row.gene1_idx)
        j = int(row.gene2_idx)
        obs = float(row.nmi)
        null = np.zeros(n_perm, dtype=np.float32)
        y = encoded[:, j]
        for p in range(n_perm):
            yp = permute_within_strata(y, strata, int(rng.integers(1, 2**31 - 1)))
            null[p] = pair_nmi(encoded[:, i], yp, int(k_eff[i]), int(k_eff[j]))
        emp_p = (1.0 + np.sum(null >= obs)) / (n_perm + 1.0)
        rows.append(
            {
                "gene1_idx": i,
                "gene2_idx": j,
                "empirical_p": emp_p,
                "null_mean": float(null.mean()),
                "null_sd": float(null.std(ddof=1)) if n_perm > 1 else np.nan,
                "nmi_z": float((obs - null.mean()) / null.std(ddof=1)) if n_perm > 1 and null.std(ddof=1) > 0 else np.nan,
            }
        )
    return pd.DataFrame(rows)


def top_allele_pairs(encoded, genes, gene_i: int, gene_j: int, min_count: int, max_rows: int) -> list[dict]:
    x = encoded[:, gene_i]
    y = encoded[:, gene_j]
    rows = []
    n = len(x)
    ux, cx = np.unique(x[x > 0], return_counts=True)
    uy, cy = np.unique(y[y > 0], return_counts=True)
    x_counts = dict(zip(ux.tolist(), cx.tolist()))
    y_counts = dict(zip(uy.tolist(), cy.tolist()))
    for a, ca in x_counts.items():
        if ca < min_count:
            continue
        xa = x == a
        for b, cb in y_counts.items():
            if cb < min_count:
                continue
            yb = y == b
            joint = int(np.sum(xa & yb))
            exp = ca * cb / n
            delta = joint - exp
            pa = ca / n
            pb = cb / n
            pab = joint / n
            mi = 0.0
            for aa, bb, c in [
                (1, 1, joint),
                (1, 0, ca - joint),
                (0, 1, cb - joint),
                (0, 0, n - ca - cb + joint),
            ]:
                if c <= 0:
                    continue
                pcell = c / n
                px = pa if aa else 1 - pa
                py = pb if bb else 1 - pb
                mi += pcell * np.log(pcell / (px * py))
            hx = -(pa * np.log(pa) + (1 - pa) * np.log(1 - pa)) if 0 < pa < 1 else 0.0
            hy = -(pb * np.log(pb) + (1 - pb) * np.log(1 - pb)) if 0 < pb < 1 else 0.0
            nmi = mi / ((hx + hy) / 2.0) if hx + hy > 0 else 0.0
            rows.append(
                {
                    "gene1": genes[gene_i],
                    "allele1_code": int(a),
                    "gene2": genes[gene_j],
                    "allele2_code": int(b),
                    "count1": int(ca),
                    "count2": int(cb),
                    "joint": joint,
                    "expected_joint": float(exp),
                    "delta_joint": float(delta),
                    "direction": "cooccur" if delta > 0 else "mutual_exclusion" if delta < 0 else "neutral",
                    "allele_pair_nmi": float(max(0.0, nmi)),
                    "weighted_nmi": float(max(0.0, nmi) * ((ca + cb) / (2 * n))),
                }
            )
    rows.sort(key=lambda r: abs(r["delta_joint"]), reverse=True)
    return rows[:max_rows]


def qvalues_bh(p: np.ndarray) -> np.ndarray:
    order = np.argsort(p)
    ranked = p[order]
    q = ranked * len(p) / (np.arange(len(p)) + 1)
    q = np.minimum.accumulate(q[::-1])[::-1]
    out = np.empty_like(q)
    out[order] = np.minimum(q, 1.0)
    return out


def plot_qq(p: np.ndarray, out_base: Path, title: str) -> None:
    p = np.asarray(p, dtype=float)
    p = p[np.isfinite(p) & (p >= 0) & (p <= 1)]
    if len(p) == 0:
        return
    positive = p[p > 0]
    floor = positive.min() / 10.0 if len(positive) else np.nextafter(0, 1)
    p[p <= 0] = max(floor, np.nextafter(0, 1))
    p = np.sort(p)
    n = len(p)
    expected = -np.log10((np.arange(1, n + 1) - 0.5) / n)
    observed = -np.log10(p)
    xlim = np.ceil((-np.log10(0.5 / n)) * 1.05 * 2) / 2
    ylim = np.ceil(observed.max() * 1.05)
    sns.set_theme(style="white", context="talk")
    fig, ax = plt.subplots(figsize=(7.2, 7.0))
    ax.scatter(expected, observed, s=13, alpha=0.35, linewidths=0, color="black")
    ax.plot([0, xlim], [0, xlim], color="red", linewidth=2)
    ax.set_xlim(0, xlim)
    ax.set_ylim(0, ylim)
    ax.set_xlabel(r"Expected $-\log_{10}(pvalue)$")
    ax.set_ylabel(r"Observed $-\log_{10}(pvalue)$")
    ax.set_title(title)
    sns.despine(fig=fig, ax=ax)
    fig.tight_layout()
    fig.savefig(out_base.with_suffix(".png"), dpi=260)
    fig.savefig(out_base.with_suffix(".pdf"))
    plt.close(fig)


def plot_elbow(df: pd.DataFrame, x: str, y: str, out_base: Path) -> float:
    sub = df[[x, y]].dropna().copy()
    sub = sub[np.isfinite(sub[x]) & np.isfinite(sub[y])]
    edges = np.unique(sub[x].quantile(np.linspace(0, 1, 81)).to_numpy())
    sub["bin"] = pd.cut(sub[x], edges, include_lowest=True, duplicates="drop")
    trend = sub.groupby("bin", observed=True).agg(x_mid=(x, "median"), y_q90=(y, lambda s: s.quantile(0.90))).dropna()
    pts = trend[["x_mid", "y_q90"]].to_numpy()
    if len(pts) < 3:
        elbow = float(sub[x].median())
    else:
        a = pts[0]
        b = pts[-1]
        denom = np.linalg.norm(b - a)
        dist = np.abs(np.cross(b - a, pts - a) / denom) if denom > 0 else np.zeros(len(pts))
        elbow = float(pts[int(np.argmax(dist)), 0])
    sns.set_theme(style="white", context="talk")
    fig, ax = plt.subplots(figsize=(8.6, 6.2))
    sample = sub.sample(min(len(sub), 200_000), random_state=42)
    ax.scatter(sample[x], sample[y], s=3, alpha=0.08, linewidths=0, color="#222222", rasterized=True)
    ax.plot(trend["x_mid"], trend["y_q90"], color="#d62728", linewidth=2.2)
    ax.axvline(elbow, color="#1f77b4", linestyle="--", linewidth=1.6)
    ax.set_xlabel(x)
    ax.set_ylabel(y)
    ax.set_ylim(bottom=0)
    sns.despine(fig=fig, ax=ax)
    fig.tight_layout()
    fig.savefig(out_base.with_suffix(".png"), dpi=260)
    fig.savefig(out_base.with_suffix(".pdf"))
    plt.close(fig)
    return elbow


def main() -> None:
    parser = argparse.ArgumentParser(description="cgMLST allelic association prototype with metadata-stratified permutation.")
    parser.add_argument("--profile", type=Path, default=Path("Leca/0427_1.profile"))
    parser.add_argument("--metadata", type=Path, default=Path("Leca/Kp.mp.merge"))
    parser.add_argument("--covariates", nargs="*", default=["Geography", "Region", "Sub_Species"])
    parser.add_argument("--tree", type=Path, default=Path("Leca/Kp.labelled.nwk"))
    parser.add_argument("--tree-pc-strata", type=int, default=0)
    parser.add_argument("--tree-pc-bins", type=int, default=4)
    parser.add_argument("--outdir", type=Path, default=Path("Leca/cgmlst_allelic_scan_0427_1"))
    parser.add_argument("--core-min", type=float, default=0.95)
    parser.add_argument("--allele-min-count", type=int, default=20)
    parser.add_argument("--min-stratum-size", type=int, default=10)
    parser.add_argument("--top-gene-pairs", type=int, default=20000)
    parser.add_argument("--n-perm", type=int, default=0)
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--allele-detail-top-pairs", type=int, default=500)
    parser.add_argument("--allele-detail-per-gene-pair", type=int, default=20)
    args = parser.parse_args()

    args.outdir.mkdir(parents=True, exist_ok=True)
    profile = pd.read_csv(args.profile, sep="\t", index_col=0).fillna(0)
    if args.tree_pc_strata > 0:
        strata = tree_pc_strata(
            args.tree,
            profile.index,
            args.tree_pc_strata,
            args.tree_pc_bins,
            args.min_stratum_size,
            args.outdir,
        )
        strata_mode = f"treePC{args.tree_pc_strata}_bins{args.tree_pc_bins}"
    else:
        strata = build_strata(profile.index, args.metadata, args.covariates, args.min_stratum_size)
        strata_mode = "metadata" if args.covariates else "none"
    encoded, genes, k_eff, allele_summary = encode_core_alleles(profile, args.core_min, args.allele_min_count)
    allele_summary.to_csv(args.outdir / "core_gene_allele_summary.tsv", sep="\t", index=False)

    g1, g2, scores, g_stats, dfs = all_pair_nmi_gtest(encoded, k_eff, strata)
    p_gtest = chi2.sf(g_stats.astype(np.float64), dfs)
    order_vals = np.array([ref_order(g) for g in genes], dtype=np.int32)
    order_dist = np.abs(order_vals[g1] - order_vals[g2]).astype(np.float32)
    pairs = pd.DataFrame(
        {
            "gene1_idx": g1,
            "gene2_idx": g2,
            "gene1": np.array(genes, dtype=object)[g1],
            "gene2": np.array(genes, dtype=object)[g2],
            "nmi": scores,
            "stratified_g": g_stats,
            "stratified_df": dfs,
            "stratified_p": p_gtest,
            "order_dist": order_dist,
            "gene1_kept_alleles": k_eff[g1] - 1,
            "gene2_kept_alleles": k_eff[g2] - 1,
        }
    )
    elbow = plot_elbow(pairs, "order_dist", "nmi", args.outdir / "elbow_nmi_vs_order_dist")
    pairs["passes_order_elbow"] = pairs["order_dist"] >= elbow
    pairs["stratified_q"] = qvalues_bh(pairs["stratified_p"].to_numpy())
    top = pairs[pairs["passes_order_elbow"]].nsmallest(args.top_gene_pairs, "stratified_p").copy()
    top["stratified_q_top_subset"] = qvalues_bh(top["stratified_p"].to_numpy())
    if args.n_perm > 0:
        emp = empirical_for_candidates(encoded, k_eff, top, strata, args.n_perm, args.seed)
        top = top.merge(emp, on=["gene1_idx", "gene2_idx"], how="left")
        top["empirical_q"] = qvalues_bh(top["empirical_p"].to_numpy())
    top.to_csv(args.outdir / "top_gene_pair_conditioned_nmi.tsv", sep="\t", index=False)
    plot_qq(pairs["stratified_p"].to_numpy(), args.outdir / "qq_stratified_gtest_all_gene_pairs", "All cgMLST gene-pair stratified G-test p QQ")
    plot_qq(top["stratified_p"].to_numpy(), args.outdir / "qq_stratified_gtest_top_gene_pairs", "Top cgMLST gene-pair stratified G-test p QQ")
    if args.n_perm > 0:
        plot_qq(top["empirical_p"].to_numpy(), args.outdir / "qq_empirical_p_top_gene_pairs", "Top cgMLST gene-pair empirical p QQ")

    details = []
    detail_rank_col = "empirical_p" if "empirical_p" in top.columns else "stratified_p"
    for row in top.nsmallest(args.allele_detail_top_pairs, detail_rank_col).itertuples(index=False):
        details.extend(top_allele_pairs(encoded, genes, int(row.gene1_idx), int(row.gene2_idx), args.allele_min_count, args.allele_detail_per_gene_pair))
    detail_df = pd.DataFrame(details)
    if not detail_df.empty:
        detail_df.to_csv(args.outdir / "top_allele_pair_details.tsv", sep="\t", index=False)

    summary = pd.DataFrame(
        [
            {
                "n_genomes": encoded.shape[0],
                "n_core_genes": encoded.shape[1],
                "allele_min_count": args.allele_min_count,
                "sum_kept_alleles": int((k_eff - 1).sum()),
                "gene_pairs_total": len(pairs),
                "order_elbow": elbow,
                "gene_pairs_after_order_elbow": int(pairs["passes_order_elbow"].sum()),
                "top_gene_pairs_permuted": len(top),
                "n_perm": args.n_perm,
                "n_strata": int(len(np.unique(strata))),
                "strata_mode": strata_mode,
                "covariates": ",".join(args.covariates),
                "min_stratified_p": float(top["stratified_p"].min()),
                "min_stratified_q_all_pairs": float(top["stratified_q"].min()),
                "min_empirical_p": float(top["empirical_p"].min()) if "empirical_p" in top else np.nan,
                "min_empirical_q": float(top["empirical_q"].min()) if "empirical_q" in top else np.nan,
            }
        ]
    )
    summary.to_csv(args.outdir / "summary.tsv", sep="\t", index=False)
    print(summary.to_string(index=False))


if __name__ == "__main__":
    main()
