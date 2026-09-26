#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Cinch v8 HPC Pipeline:
Phase 1: High-speed Reference Mapping & Profile Generation
Phase 2: Numba-Accelerated Phylogeny-aware LD Scan (wgMLST & cgMLST Stratified Z-score)

v8 change:
- Deduplicate identical reference CDS sequences during preprocessing so one
  biological reference sequence is not assigned multiple ref_gene IDs.
"""

import argparse
import glob
import hashlib
import logging
import os
import subprocess
import sys
import warnings
import shutil
from collections import OrderedDict, defaultdict
from concurrent.futures import ThreadPoolExecutor, as_completed
from itertools import combinations
from typing import Dict, List, Tuple, Optional, Iterable, Set
import re

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from Bio import Phylo, SeqIO
from sklearn.metrics import normalized_mutual_info_score
from sklearn.decomposition import PCA
from sklearn.cluster import KMeans
from numba import njit, prange
from scipy.stats import t as t_dist

try:
    from configure import scheme_info
    from uberBlast import uberBlast, readFastq, rc
except ImportError:
    logging.warning("Modules 'configure' or 'uberBlast' not found. Ensure they are in the PYTHONPATH.")

warnings.simplefilter(action='ignore')
logging.basicConfig(level=logging.INFO, format='%(asctime)s %(levelname)s: %(message)s')


# =========================================================
# =================  PHASE 1: FRONTEND  ===================
# =========================================================

def get_md5(value, dtype=str):
    m = hashlib.md5(str(value).encode()).hexdigest()
    return m if dtype == str else int(m, 16)

def normalize_interval(s, e):
    start, end = min(int(s), int(e)), max(int(s), int(e))
    strand = '+' if int(s) <= int(e) else '-'
    return start, end, strand

def canonical_genome_name(name: str) -> str:
    s = os.path.basename(str(name).strip())
    m = re.search(r'(GC[AF]_\d+)', s)
    return m.group(1) if m else s

def to_canonical_set(values: Optional[Iterable[str]]) -> Set[str]:
    return {canonical_genome_name(v) for v in values} if values else set()

def preprocess(prefix, ref_cds, output_path):
    ref_cds_seq = SeqIO.parse(ref_cds, 'fasta')
    ref_base = os.path.basename(ref_cds)
    modi_ref_cds = os.path.join(output_path, f"{prefix}.{ref_base}.modi.cds")
    
    gene_count = 1
    total_count = 0
    skipped_duplicates = 0
    seen_seq_md5 = {}
    with open(modi_ref_cds, 'w') as fout:
        for record in ref_cds_seq:
            total_count += 1
            seq = str(record.seq).upper()
            seq_md5 = get_md5(seq)
            if seq_md5 in seen_seq_md5:
                skipped_duplicates += 1
                continue
            seen_seq_md5[seq_md5] = f"{prefix}:ref_gene_{gene_count}"
            fout.write(f">{prefix}:ref_gene_{gene_count}_1\n{seq}\n")
            gene_count += 1
    logging.info(
        "Reference CDS deduplicated: %s unique sequences kept from %s records; %s duplicates skipped.",
        gene_count - 1, total_count, skipped_duplicates
    )

    genes_file, md5_out, pos_out = nomenclature(prefix, modi_ref_cds, modi_ref_cds, output_path)
    return modi_ref_cds, gene_count - 1, genes_file, md5_out, pos_out

def nomenclature(prefix, query, refseqs, output_path, n_thread=1):
    blastab = uberBlast(
        f"-r {query} -q {refseqs} -f --blastn --diamond --min_id {scheme_info['min_iden'] - 0.05} "
        f"--min_ratio 0.05 -t {n_thread} -p -s 2 -e 21,21 -m --merge_gap 300".split()
    )

    merges = {}
    for b in blastab.T[16]:
        if len(b) > 4:
            merges[tuple(b[3:])] = b[:3]

    bsn = {b[15]: b[:15] for b in blastab}
    for ids, (score, iden, size) in merges.items():
        bs = np.array([bsn[i] for i in ids], dtype=object)
        if np.unique(bs.T[1]).size == 1:
            bs[0][2], bs[0][3], bs[0][7], bs[0][9], bs[0][11], bs[0][14] = (
                iden, size, bs[-1][7], bs[-1][9], score, 'COMPLEX'
            )
        for i in ids[1:]:
            bsn.pop(i)
        bsn[ids[0]] = bs[0]

    blastab = np.array(list(bsn.values()), dtype=object)
    blastab = blastab[
        (blastab.T[2] >= scheme_info['min_iden']) &
        ((blastab.T[7] - blastab.T[6] + 1) >= scheme_info['min_frag'] * blastab.T[12])
    ]

    q_base = os.path.basename(query)
    cinch_md5 = os.path.join(output_path, f"{prefix}_{q_base}.md5")
    genes_file = os.path.join(output_path, f"{prefix}_{q_base}.genes")
    pos_file = os.path.join(output_path, f"{prefix}_{q_base}.positions.tsv")

    if blastab.size == 0:
        open(cinch_md5, 'w').close()
        open(genes_file, 'w').close()
        with open(pos_file, 'w') as fout:
            fout.write("genome\tgene\tcontig\tstart\tend\tstrand\tgene_length\tcontig_length\tmd5\tflag\tidentity\tcoordinates\n")
        return genes_file, cinch_md5, pos_file

    blastab.T[11] = blastab.T[11] * (blastab.T[7] - blastab.T[6] + 1) / blastab.T[12]
    blastab = blastab[np.lexsort((blastab[:, 8], blastab[:, 1]))]

    for i, b0 in enumerate(blastab[:-1]):
        if b0[0] == '': continue
        s0, e0 = sorted(b0[8:10])
        todel = []
        for b1 in blastab[i + 1:]:
            s1, e1 = sorted(b1[8:10])
            if b0[1] != b1[1] or e0 < s1: break
            ovl = min(e0, e1) - max(s0, s1) + 1
            if ovl >= 0.5 * (e0 - s0 + 1) or ovl >= 0.5 * (e1 - s1 + 1):
                sc0, sc1 = abs(b0[11]), abs(b1[11])
                g0, g1 = b0[0].rsplit('_', 1)[0], b1[0].rsplit('_', 1)[0]
                if b0[2] < b1[2] * scheme_info['max_iden'] or (b1[2] >= b0[2] * scheme_info['max_iden'] and (sc0 < sc1 or (sc0 == sc1 and b0[0] > b1[0]))):
                    b0[11] = -sc0
                    if g0 == g1 or sc0 < sc1 * scheme_info['max_iden'] or b0[2] < b1[2] * scheme_info['max_iden']:
                        b0[0] = ''
                        break
                else:
                    b1[11] = -sc1
                    if g0 == g1 or sc1 < sc0 * scheme_info['max_iden'] or b1[2] < b0[2] * scheme_info['max_iden']:
                        todel.append(b1)
        if b0[0] and todel:
            for b1 in todel: b1[0] = ''

    blastab = blastab[blastab.T[0] != '']
    if blastab.size == 0:
        open(cinch_md5, 'w').close()
        open(genes_file, 'w').close()
        with open(pos_file, 'w') as fout:
            fout.write("genome\tgene\tcontig\tstart\tend\tstrand\tgene_length\tcontig_length\tmd5\tflag\tidentity\tcoordinates\n")
        return genes_file, cinch_md5, pos_file

    blastab = blastab[np.lexsort((-blastab.T[11], [b.rsplit('_', 1)[0] for b in blastab.T[0]]))]
    alleles = OrderedDict()
    cinch_fa = os.path.join(output_path, f"{prefix}_{q_base}.cinch.fasta")

    with open(cinch_fa, 'w') as fasta_out, open(cinch_md5, 'w') as md5_out, open(pos_file, 'w') as pos_out:
        pos_out.write("genome\tgene\tcontig\tstart\tend\tstrand\tgene_length\tcontig_length\tmd5\tflag\tidentity\tcoordinates\n")
        for bsn in blastab:
            gene = bsn[0].rsplit('_', 1)[0]
            if gene in alleles:
                if alleles[gene]['score'] * scheme_info['max_iden'] > bsn[11]: continue
                alleles[gene]['coordinates'].append((bsn[1], int(bsn[8]), int(bsn[9])))
                alleles[gene]['flag'] |= 32
                continue

            flag = 0
            if bsn[6] > 1 or bsn[7] < bsn[12]:
                flag = 64
                if bsn[14] != 'COMPLEX':
                    if bsn[6] > 1: bsn[14] = f"{bsn[6] - 1}D{bsn[14]}"
                    if bsn[7] < bsn[12]: bsn[14] = f"{bsn[14]}{bsn[12] - bsn[7]}D"

            alleles[gene] = {
                'gene_name': gene, 'CIGAR': f"{bsn[0]}:{bsn[14]}",
                'reference': q_base, 'identity': float(bsn[2]),
                'coordinates': [(bsn[1], int(bsn[8]), int(bsn[9]))],
                'flag': flag, 'score': bsn[11]
            }

        seq, qual = readFastq(query)
        contig_lengths = {k: len(v) for k, v in seq.items()}
        remove_keys = []

        for gene, allele in sorted(alleles.items()):
            if allele['flag'] & 96 == 96:
                remove_keys.append(gene)
                continue

            if allele['flag'] & 32 > 0:
                allele['sequence'] = 'DUPLICATED'
            else:
                c, s, e = allele['coordinates'][0]
                ss = seq[c][s - 1:e] if s < e else rc(seq[c][e - 1:s])
                qs = min(qual[c][s - 1:e] if s < e else qual[c][e - 1:s]) if qual else 0
                allele['sequence'] = ss
                if qs < 10: allele['flag'] |= 2
                if allele['flag'] == 0: allele['flag'] = 1

            allele['coordinates_str'] = ','.join([f"{c[0]}:{c[1]}..{c[2]}" for c in allele['coordinates']])
            md5_value = ('' if allele['flag'] < 16 else '-') + get_md5(allele['sequence'])
            
            fasta_out.write(f">{gene}_in_{query}|{allele['coordinates_str']}|md5={md5_value}\n{allele['sequence']}\n")
            md5_out.write(f"{gene}_in_{query}\t{md5_value}\n")

            gene_len = len(allele['sequence']) if allele['sequence'] != 'DUPLICATED' else None
            for contig, s, e in allele['coordinates']:
                start, end, strand = normalize_interval(s, e)
                local_gene_len = gene_len if gene_len is not None else (end - start + 1)
                contig_len = contig_lengths.get(contig, np.nan)
                pos_out.write(f"{q_base}\t{gene}\t{contig}\t{start}\t{end}\t{strand}\t{local_gene_len}\t{contig_len}\t"
                              f"{md5_value}\t{allele['flag']}\t{allele['identity']}\t{allele['coordinates_str']}\n")

        for k in remove_keys: alleles.pop(k, None)
    
    with open(genes_file, 'w') as f:
        f.write('\n'.join(alleles.keys()) + '\n')

    return genes_file, cinch_md5, pos_file

# --- Fast File Merging with shutil.copyfileobj ---
def merge_files_optimized(prefix, work_dir, output_suffix, file_suffix, skip_header=False):
    files = glob.glob(os.path.join(work_dir, f"{prefix}*{file_suffix}"))
    out_name = os.path.join(work_dir, f"{prefix}{output_suffix}")
    
    with open(out_name, 'w') as fout:
        header_written = False
        for fname in files:
            with open(fname, 'r') as fin:
                if skip_header:
                    header = fin.readline()
                    if not header_written:
                        fout.write(header)
                        header_written = True
                shutil.copyfileobj(fin, fout)
    return out_name

def gene_conservative(merged_genes, total_num, ratio, work_dir, prefix):
    genes = defaultdict(int)
    with open(merged_genes, 'r') as gin:
        for line in gin: genes[line.rstrip()] += 1

    g_sum = os.path.join(work_dir, f"{prefix}.conservation.summary")
    cgs = os.path.join(work_dir, f"{prefix}.core_genes")
    num_cg = 0
    with open(g_sum, 'w') as gout, open(cgs, 'w') as cgout:
        for g, count in genes.items():
            if count >= ratio * total_num:
                cgout.write(g + '\n')
                num_cg += 1
            gout.write(f"{g}\t{round(count / total_num, 4)}\n")
    return g_sum, cgs, num_cg

def md5_numerate(gene_title_with_md5, core_genes, ref_md5_file):
    mapping_dict = defaultdict(dict)
    with open(ref_md5_file, 'r') as f:
        for line in f:
            ref_gene, md5 = line.rstrip().split('\t')
            title = ref_gene.split('_in_')[0]
            if title in core_genes:
                mapping_dict[title][md5] = 1

    with open(gene_title_with_md5, 'r') as f:
        for line in f:
            part = line.rstrip().split('\t')
            gene = part[0].split('_in_')[0]
            if gene in core_genes:
                md5 = part[-1].replace('-', '') if '-' in part[-1] else part[-1]
                if md5 not in mapping_dict[gene]:
                    mapping_dict[gene][md5] = len(mapping_dict[gene]) + 1
    return mapping_dict

def generate_profile(mapping_dict, gene_title_with_md5, prefix, output_path, exclude_genomes=None):
    exclude_genomes = to_canonical_set(exclude_genomes)
    genome_dict = defaultdict(dict)

    with open(gene_title_with_md5, 'r') as f:
        for line in f:
            part = line.rstrip().split('\t')
            md5_value = part[-1]
            gene, genome = part[0].split('_in_')[0], part[0].split('_in_')[1]
            
            if canonical_genome_name(genome) in exclude_genomes or gene not in mapping_dict: continue
            
            cleaned_md5 = md5_value.replace('-', '') if '-' in md5_value else md5_value
            val = mapping_dict[gene].get(cleaned_md5, 0)
            genome_dict[genome][gene] = -val if '-' in md5_value else val

    df = pd.DataFrame.from_dict(genome_dict, orient='index').fillna(0).astype(int)
    df.index.name = 'genome'
    df.reset_index(inplace=True)
    ordered_cols = ['genome'] + list(mapping_dict.keys())
    for col in ordered_cols:
        if col not in df.columns: df[col] = 0
    df = df[ordered_cols]
    
    out_profile = os.path.join(output_path, f"{prefix}.profile")
    df.to_csv(out_profile, sep='\t', index=False)
    return out_profile


# =========================================================
# =================  PHASE 2: BACKEND =====================
# Numba-Accelerated HPC LD Scan
# =========================================================

def get_interval_dist_backend(s1, e1, s2, e2, c_len=None, mode='linear'):
    st1, en1 = min(s1, e1), max(s1, e1)
    st2, en2 = min(s2, e2), max(s2, e2)
    if not (en1 < st2 or en2 < st1): return 0
    dist = st2 - en1 - 1 if en1 < st2 else st1 - en2 - 1
    if mode == 'circular' and c_len:
        wrap_dist = (min(st1, st2) + c_len) - max(en1, en2) - 1
        return min(dist, wrap_dist)
    return dist

def calculate_p_values(r_matrix: np.ndarray, n_samples: int, n_covariates: int = 0) -> np.ndarray:
    df = n_samples - 2 - n_covariates
    if df <= 0: return np.full_like(r_matrix, np.nan)
    with np.errstate(divide='ignore', invalid='ignore'):
        t_stat = r_matrix * np.sqrt(df / (1.0 - np.square(r_matrix) + 1e-9))
        p_matrix = 2 * t_dist.sf(np.abs(t_stat), df=df)
    np.fill_diagonal(p_matrix, 0)
    return p_matrix.astype(np.float32)

@njit(fastmath=True)
def _entropy(arr):
    n = len(arr)
    if n == 0: return 0.0
    counts = np.bincount(arr)
    probs = counts[counts > 0] / n
    return -np.sum(probs * np.log(probs + 1e-12))

@njit(parallel=True, fastmath=True)
def compute_stratified_nmi_z(X_int, strata, n_perm, min_joint):
    n_samples, n_genes = X_int.shape
    nmi_obs = np.full((n_genes, n_genes), np.nan, dtype=np.float32)
    z_scores = np.full((n_genes, n_genes), np.nan, dtype=np.float32)

    for i in prange(n_genes):
        for j in range(i + 1, n_genes):
            mask = (X_int[:, i] >= 0) & (X_int[:, j] >= 0)
            v1 = X_int[mask, i]
            v2 = X_int[mask, j]
            if len(v1) < min_joint: continue

            h1, h2 = _entropy(v1), _entropy(v2)
            if h1 + h2 == 0: continue
            
            v2_max = v2.max() + 1
            combined = v1 * v2_max + v2
            mi_obs = max(0.0, h1 + h2 - _entropy(combined))
            obs_nmi = mi_obs / ((h1 + h2) / 2.0)
            nmi_obs[i, j] = obs_nmi

            perm_values = np.zeros(n_perm, dtype=np.float32)
            cur_strata = strata[mask]
            
            unique_strata = []
            for s in cur_strata:
                if s not in unique_strata: unique_strata.append(s)

            for p in range(n_perm):
                v2_perm = v2.copy()
                for s_id in unique_strata:
                    idx = np.where(cur_strata == s_id)[0]
                    if len(idx) > 1:
                        perm_idx = np.random.permutation(len(idx))
                        v2_perm[idx] = v2[idx[perm_idx]]
                
                c_perm = v1 * v2_max + v2_perm
                mi_p = max(0.0, h1 + h2 - _entropy(c_perm))
                perm_values[p] = mi_p / ((h1 + h2) / 2.0)
            
            m_null, s_null = np.mean(perm_values), np.std(perm_values)
            z_scores[i, j] = (obs_nmi - m_null) / s_null if s_null > 0 else 0.0

    return nmi_obs, z_scores

def compute_adjusted_r_matrix(X: np.ndarray, n_pcs: int = 5) -> np.ndarray:
    logging.info(f"Extracting {n_pcs} Phylo-PCs and projecting residuals...")
    X_centered = X - np.mean(X, axis=0)
    pca = PCA(n_components=n_pcs)
    U = pca.fit_transform(X_centered)  
    Q, _ = np.linalg.qr(U)
    X_proj = Q @ (Q.T @ X_centered)
    X_res = X_centered - X_proj
    with np.errstate(divide='ignore', invalid='ignore'):
        r_adj = np.corrcoef(X_res.T).astype(np.float32)
    return r_adj

def compute_vectorized_phylo_dist(X: np.ndarray, D: np.ndarray, joint_counts: np.ndarray) -> np.ndarray:
    S = (X.T @ D) @ X
    pair_denom = (joint_counts * (joint_counts - 1)) / 2.0
    with np.errstate(divide='ignore', invalid='ignore'):
        return (S / (2.0 * pair_denom)).astype(np.float32)

def get_physical_stats(positions_file: str, genes_to_keep: List[str], dist_mode='linear'):
    logging.info(f"Reconstructing GLOBAL Order and distances (Mode: {dist_mode})...")
    full_df = pd.read_csv(positions_file, sep='\t')
    full_df['genome'] = full_df['genome'].map(canonical_genome_name)
    full_df = full_df.sort_values(['genome', 'contig', 'start'])
    full_df['locus_start'] = full_df[['start', 'end']].min(axis=1)
    full_df['locus_end'] = full_df[['start', 'end']].max(axis=1)

    gene_set = set(genes_to_keep)
    pos_df = full_df[full_df['gene'].isin(gene_set)].copy()
    
    n = len(genes_to_keep)
    g2idx = {g: i for i, g in enumerate(genes_to_keep)}
    bp_sum = np.zeros((n, n), dtype=np.float64)
    order_sum = np.zeros((n, n), dtype=np.float64)
    same_replicon_count = np.zeros((n, n), dtype=np.int32)

    for (genome, contig), gdf in pos_df.groupby(['genome', 'contig']):
        gdf = gdf.sort_values(['locus_start', 'locus_end']).copy()
        locus_orders = []
        locus_order = -1
        last_start = None
        last_end = None
        for row in gdf.itertuples():
            start = int(row.locus_start)
            end = int(row.locus_end)
            if last_start is None or start != last_start or end != last_end:
                locus_order += 1
                last_start, last_end = start, end
            locus_orders.append(locus_order)
        gdf['locus_order'] = locus_orders
        total_loci_on_contig = locus_order + 1

        records = gdf.to_dict('records')
        if len(records) < 2: continue
        c_len = records[0].get('contig_length', None) if dist_mode == 'circular' else None
        local_min_bp = {}; local_min_order = {}
        
        for r1, r2 in combinations(records, 2):
            i, j = sorted([g2idx[r1['gene']], g2idx[r2['gene']]])
            pair = (i, j)
            bp = get_interval_dist_backend(r1['start'], r1['end'], r2['start'], r2['end'], c_len, dist_mode)
            od_lin = abs(r1['locus_order'] - r2['locus_order'])
            od = min(od_lin, total_loci_on_contig - od_lin) if dist_mode == 'circular' else od_lin
            
            if pair not in local_min_bp or bp < local_min_bp[pair]: local_min_bp[pair] = bp
            if pair not in local_min_order or od < local_min_order[pair]: local_min_order[pair] = od
            
        for (i, j), bp in local_min_bp.items():
            bp_sum[i, j] += bp
            order_sum[i, j] += local_min_order[(i, j)]
            same_replicon_count[i, j] += 1
            
    return bp_sum, order_sum, same_replicon_count

def run_hpc_ld_scan(profile_file, positions_file, tree_file, output_prefix, 
                    assoc_mode, n_pcs, min_prev, max_prev, min_joint_count, distance_mode, make_plots):
    logging.info(f"Phase 2 HPC Scan Mode: {assoc_mode} | Phylo-PCs: {n_pcs}")
    
    profile_df = pd.read_csv(profile_file, sep='\t', index_col=0).fillna(0)
    profile_df.index = [canonical_genome_name(x) for x in profile_df.index]
    n_samples = profile_df.shape[0]
    
    bin_df = (profile_df > 0).astype(np.float32)
    prev_series = bin_df.mean(axis=0)
    
    keep_mask = (prev_series >= min_prev) if assoc_mode == 'cgMLST' else (prev_series >= min_prev) & (prev_series <= max_prev)
    genes = prev_series[keep_mask].index.tolist()
    gene_prevs = prev_series[keep_mask].values
    X_bin = bin_df[genes].values
    genomes = bin_df.index.tolist()
    n = len(genes)

    # 1. 相关性计算
    joint_counts = (X_bin.T @ X_bin).astype(np.int32)
    r_raw = np.corrcoef(X_bin.T).astype(np.float32)
    r2_raw = np.square(r_raw)
    p_raw = calculate_p_values(r_raw, n_samples) if assoc_mode == 'wgMLST' else None
    
    r_adj = compute_adjusted_r_matrix(X_bin, n_pcs=n_pcs)
    r2_adj = np.square(r_adj)
    p_adj = calculate_p_values(r_adj, n_samples, n_covariates=n_pcs) if assoc_mode == 'wgMLST' else None

    # 2. cgMLST 专项 NMI & Z-score
    nmi_obs, nmi_z = None, None
    if assoc_mode == 'cgMLST':
        logging.info("Encoding alleles & clustering for stratified permutation...")
        X_int = np.zeros((n_samples, n), dtype=np.int32)
        for i, g in enumerate(genes):
            X_int[:, i] = pd.factorize(profile_df[g])[0]
        
        pca_cluster = PCA(n_components=n_pcs)
        u_cluster = pca_cluster.fit_transform(X_bin)
        kmeans = KMeans(n_clusters=max(2, min(8, n_samples//10)), random_state=42, n_init=10)
        strata = kmeans.fit_predict(u_cluster).astype(np.int32)
        
        logging.info("Computing NMI and Z-score (Numba Accelerated)...")
        nmi_obs, nmi_z = compute_stratified_nmi_z(X_int, strata, n_perm=100, min_joint=min_joint_count)

    # 3. 系统发育距离与物理距离
    logging.info("Calculating mean patristic distances...")
    tree = Phylo.read(tree_file, 'newick')
    tip_map = {canonical_genome_name(t.name): t.name for t in tree.get_terminals() if t.name}
    D = np.zeros((n_samples, n_samples), dtype=np.float32)
    for i, g1 in enumerate(genomes):
        for j in range(i + 1, n_samples):
            if g1 in tip_map and genomes[j] in tip_map:
                D[i, j] = D[j, i] = tree.distance(tip_map[g1], tip_map[genomes[j]])
    
    phylo_matrix = compute_vectorized_phylo_dist(X_bin, D, joint_counts)
    bp_sum_mtx, order_sum_mtx, same_cnt_mtx = get_physical_stats(positions_file, genes, distance_mode)

    # 4. 导出结果
    out_file = f"{output_prefix}.ldlike_pairs.tsv"
    logging.info(f"Saving combined results to {out_file}...")
    with open(out_file, 'w') as f:
        f.write("gene1\tgene2\tprev1\tprev2\tr2_raw\tp_raw\tr2_adj\tp_adj\tnmi\tnmi_zscore\tn_joint\tmean_phylo_dist\tmean_bp_dist\tmean_order_dist\tsame_replicon_rate\n")
        for i in range(n):
            for j in range(i + 1, n):
                nj = joint_counts[i, j]
                if nj < min_joint_count: continue
                pr = p_raw[i, j] if p_raw is not None else np.nan
                pa = p_adj[i, j] if p_adj is not None else np.nan
                no = nmi_obs[i, j] if nmi_obs is not None else np.nan
                nz = nmi_z[i, j] if nmi_z is not None else np.nan
                s_cnt = same_cnt_mtx[i, j]
                m_bp = bp_sum_mtx[i, j] / s_cnt if s_cnt > 0 else np.nan
                m_order = order_sum_mtx[i, j] / s_cnt if s_cnt > 0 else np.nan
                f.write(f"{genes[i]}\t{genes[j]}\t{gene_prevs[i]:.4f}\t{gene_prevs[j]:.4f}\t"
                        f"{r2_raw[i,j]:.4f}\t{pr:.2e}\t{r2_adj[i,j]:.4f}\t{pa:.2e}\t"
                        f"{no:.4f}\t{nz:.2f}\t{nj}\t{phylo_matrix[i,j]:.6f}\t"
                        f"{m_bp:.1f}\t{m_order:.2f}\t{s_cnt/nj:.4f}\n")

    if make_plots:
        visualize_results(out_file, output_prefix, assoc_mode)


def visualize_results(tsv, prefix, mode):
    logging.info("Generating comparison plots...")
    df = pd.read_csv(tsv, sep='\t')
    if len(df) < 10: return
    if mode == 'wgMLST':
        sub = df[['mean_order_dist', 'r2_raw', 'r2_adj']].dropna()
        fig, axes = plt.subplots(1, 2, figsize=(12, 5), sharey=True)
        axes[0].hexbin(sub['mean_order_dist'], sub['r2_raw'], gridsize=40, cmap='Blues', bins='log', mincnt=1)
        axes[0].set_title("Raw R2")
        axes[1].hexbin(sub['mean_order_dist'], sub['r2_adj'], gridsize=40, cmap='Reds', bins='log', mincnt=1)
        axes[1].set_title("Phylo-Adjusted R2")
        plt.savefig(f"{prefix}.r2_correction_comparison.png", dpi=300); plt.close()
    else:
        sub = df[['mean_order_dist', 'nmi_zscore']].dropna()
        plt.figure(figsize=(8, 6))
        plt.hexbin(sub['mean_order_dist'], sub['nmi_zscore'], gridsize=40, cmap='Spectral_r', bins='log', mincnt=1)
        plt.xlabel("Order Distance"); plt.ylabel("NMI Z-score")
        plt.title("cgMLST Allele Linkage (Phylo-Corrected Z-score)")
        plt.savefig(f"{prefix}.nmi_z_decay.png", dpi=300); plt.close()

# =========================================================
# =================      MAIN LOOP      ===================
# =========================================================

def main():
    parser = argparse.ArgumentParser(description='Cinch HPC Pipeline: Fast Reference Mapping + Numba LD Scan')
    parser.add_argument('-r', '--reference', required=True, help='Reference CDS file')
    parser.add_argument('-t', '--threads', type=int, default=8, help='Number of threads')
    parser.add_argument('-p', '--prefix', default='Cinch', help='Prefix for output files')
    parser.add_argument('-C', '--conservation_ratio', type=float, default=0.95, help='Conservation ratio for core-gene (0-1)')
    
    # Phase 2 Args
    parser.add_argument('--pairwise_assoc', action='store_true', help='Run phylogeny-aware association analysis')
    parser.add_argument('--assoc_mode', choices=['cgMLST', 'wgMLST'], default='cgMLST', help='Association mode')
    parser.add_argument('--tree', required=False, help='Newick tree for phylogeny-aware association analysis')
    parser.add_argument('--distance_mode', choices=['linear', 'circular'], default='linear', help='Physical distance mode')
    parser.add_argument('--n_pcs', type=int, default=10, help='Number of phylogenetic PCs used as covariates')
    parser.add_argument('--min_joint_count', type=int, default=20, help='Minimum joint non-zero sample size for test')
    parser.add_argument('--wg_min_presence', type=float, default=0.03, help='Minimum prevalence in wgMLST mode')
    parser.add_argument('--wg_max_presence', type=float, default=0.97, help='Maximum prevalence in wgMLST mode')
    parser.add_argument('--make_plots', action='store_true', help='Generate LD decay plots')
    parser.add_argument('genomes', nargs='+', help='List of genome files')
    
    args = parser.parse_args()

    if args.pairwise_assoc and not args.tree:
        parser.error('--tree is required when --pairwise_assoc is enabled')

    work_dir = os.path.join(os.getcwd(), args.prefix)
    if not os.path.exists(work_dir):
        os.makedirs(work_dir)

    # ---------------------------------------------------------
    # PHASE 1: Reference Mapping & Profile Generation
    # ---------------------------------------------------------
    logging.info('PHASE 1: Parsing reference CDS & Mapping...')
    modified_cds, cg_count, ref_genes_file, ref_md5_out, ref_pos_out = preprocess(args.prefix, args.reference, work_dir)
    ref_basename = os.path.basename(modified_cds)
    
    with ThreadPoolExecutor(max_workers=args.threads) as executor:
        futures = {executor.submit(nomenclature, args.prefix, genome, modified_cds, work_dir): genome for genome in args.genomes}
        for future in as_completed(futures):
            future.result()

    logging.info('Mapping finished. Merging results...')
    genes_presence = merge_files_optimized(args.prefix, work_dir, '.genes.presence', '.genes')
    gene_title_md5 = merge_files_optimized(args.prefix, work_dir, '.md5_summary', '.md5')
    positions_summary = merge_files_optimized(args.prefix, work_dir, '.positions.summary.tsv', '.positions.tsv', skip_header=True)

    genes_sum, core_genes, num_cg = gene_conservative(genes_presence, len(args.genomes), args.conservation_ratio, work_dir, args.prefix)
    core_set = [line.strip() for line in open(core_genes)]
    
    mapping_dict = md5_numerate(gene_title_md5, core_set, ref_md5_out)
    output_profile = generate_profile(mapping_dict, gene_title_md5, args.prefix, work_dir, exclude_genomes={ref_basename})
    logging.info(f'Profile generation complete: {output_profile}')

    # ---------------------------------------------------------
    # PHASE 2: HPC LD-Like Scan
    # ---------------------------------------------------------
    if args.pairwise_assoc:
        logging.info('PHASE 2: Activating High-Performance Association Analysis...')
        pairwise_prefix = os.path.join(work_dir, args.prefix)
        run_hpc_ld_scan(
            profile_file=output_profile,
            positions_file=positions_summary,
            tree_file=args.tree,
            output_prefix=pairwise_prefix,
            assoc_mode=args.assoc_mode,
            n_pcs=args.n_pcs,
            min_prev=args.wg_min_presence,
            max_prev=args.wg_max_presence,
            min_joint_count=args.min_joint_count,
            distance_mode=args.distance_mode,
            make_plots=args.make_plots
        )

if __name__ == '__main__':
    main()
