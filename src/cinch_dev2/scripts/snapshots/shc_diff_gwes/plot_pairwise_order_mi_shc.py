#!/usr/bin/env python3
"""Plot distance versus MI for true and toy data, marking selection-specific SHC pairs."""

from pathlib import Path
import textwrap
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm, LinearSegmentedColormap

ROOT = Path(__file__).resolve().parents[4]
WF = ROOT / "Leca/workflows/shc_diff_gwes"
INK, MUTED, GRID = "#20252C", "#68717D", "#E2E7EB"
BLUE, ORANGE, GREY = "#3976A8", "#D46A3A", "#B8C1C9"
CMAP = LinearSegmentedColormap.from_list("density", ["#FFFFFF", "#E3EDF3", "#80ACC5", "#285D7A"])


def mi(x, y, w):
    total = w.sum()
    cells = np.array([w[(x == a) & (y == b)].sum() for a, b in [(0, 0), (0, 1), (1, 0), (1, 1)]]) / total
    px = np.array([cells[:2].sum(), cells[2:].sum()]); py = np.array([cells[[0, 2]].sum(), cells[[1, 3]].sum()])
    exp = np.array([px[0]*py[0], px[0]*py[1], px[1]*py[0], px[1]*py[1]])
    ok = (cells > 0) & (exp > 0)
    return float(np.sum(cells[ok] * np.log(cells[ok] / exp[ok])))


def informative(x, y, shc, minimum=5):
    keep, positive = [], 0
    for h in np.unique(shc):
        m = shc == h; n = int(m.sum())
        if n < minimum or not (0 < x[m].sum() < n and 0 < y[m].sum() < n): continue
        keep.append(h)
        n11=np.sum(x[m]&y[m]); n10=np.sum(x[m]&~y[m]); n01=np.sum(~x[m]&y[m]); n00=np.sum(~x[m]&~y[m])
        positive += int(n11*n00 > n10*n01)
    return keep, positive


def conditional(x, y, w, shc, keep):
    mask=np.isin(shc,keep); total=w[mask].sum(); value=0.0
    for h in keep:
        m=shc==h; value += w[m].sum()/total * mi(x[m],y[m],w[m])
    return value


def bh(p):
    order=np.argsort(p); ranked=p[order]; q=np.minimum.accumulate((ranked*len(p)/np.arange(1,len(p)+1))[::-1])[::-1]
    out=np.empty_like(q); out[order]=np.minimum(q,1); return out


def synthetic_pair_tests():
    profile=pd.read_csv(WF/'input/toy_profile.tsv.gz',sep='\t',compression='gzip',index_col=0,usecols=lambda c:c=='genome' or c.startswith('toy_epi_')).fillna(0)>0
    shc=pd.read_csv(WF/'input/sample_split.tsv',sep='\t').set_index('genome').loc[profile.index].shc.to_numpy()
    wt=pd.read_csv(ROOT/'Leca/results/profile_similarity/wg_cg_profile_0427_1/spydrpick_style_sample_weights.tsv',sep='\t').set_index('genome').loc[profile.index].wg_content_hamming_weight.to_numpy(float)
    rng=np.random.default_rng(20260713); rows=[]
    for i in range(1,21):
        a=f'toy_epi_{i:02d}_A'; b=f'toy_epi_{i:02d}_B'; x=profile[a].to_numpy(); y=profile[b].to_numpy()
        keep,positive=informative(x,y,shc); obs=conditional(x,y,wt,shc,keep)
        null=np.empty(999)
        for p in range(999):
            yp=y.copy()
            for h in keep:
                idx=np.flatnonzero(shc==h); yp[idx]=y[idx][rng.permutation(len(idx))]
            null[p]=conditional(x,yp,wt,shc,keep)
        sd=null.std(ddof=1)
        rows.append({'pair_id':'|'.join(sorted((a,b))),'gene1':a,'gene2':b,'weighted_mi':mi(x,y,wt),'conditional_mi':obs,'conditional_null_mean':null.mean(),'z':(obs-null.mean())/sd,'p':(1+np.sum(null>=obs))/1000,'informative_shc':len(keep),'positive_shc':positive,'mean_order_dist':300.0,'mean_bp_dist':300000.0})
    out=pd.DataFrame(rows); out['q']=bh(out.p.to_numpy()); out['passes_shc']=(out.q<.05)&(out.z>0)&(out.positive_shc>=2)
    out.to_csv(WF/'output/toy/synthetic_pairwise_shc_tests.tsv',sep='\t',index=False)
    return out


def style(ax):
    ax.spines[['top','right']].set_visible(False); ax.spines[['left','bottom']].set_color('#B8C0C8')
    ax.grid(axis='y',color=GRID,linewidth=.7); ax.set_axisbelow(True); ax.tick_params(colors=MUTED,labelsize=8)


def plot_distance(base, strict, toy, distance, elbow, xlabel, output_stem):
    outdir=WF/'output/figures'
    strict_ids={'|'.join(sorted((a,b))) for a,b in zip(strict.gene1,strict.gene2)}
    strict_points=base[base.pair_id.isin(strict_ids)]
    fig,axes=plt.subplots(1,2,figsize=(13,5.6),facecolor='white'); fig.subplots_adjust(left=.07,right=.98,bottom=.13,top=.78,wspace=.23)
    distance_name='gene-order distance' if distance=='mean_order_dist' else 'physical distance'
    simulated='300 order units' if distance=='mean_order_dist' else '300,000 bp'
    fig.text(.07,.965,f'Pairwise MI versus {distance_name} after SHC filtering',ha='left',va='top',fontsize=17,weight='bold',color=INK)
    subtitle=(f'Grey density: all Kp locus pairs with callable {distance_name}. Blue: only real pairs passing the {distance_name}-specific '
              f'99.5th-percentile elbow, SHC permutation, and ARACNE filters. Orange: synthetic A-B pairs passing an independent '
              f'999-permutation SHC test; simulated distance is {simulated}. Dotted line: 99.5th-percentile elbow at {elbow:,.2f}.')
    fig.text(.07,.91,textwrap.fill(subtitle,145),ha='left',va='top',fontsize=9,color=MUTED)
    valid=base[distance].notna()&base.weighted_mi.notna()
    for ax,title,show_toy in [(axes[0],'A  True Kp data',False),(axes[1],'B  Toy profile with 20 injected pairs',True)]:
        ax.hexbin(base.loc[valid,distance],base.loc[valid,'weighted_mi'],gridsize=(105,60),mincnt=1,cmap=CMAP,norm=LogNorm(),linewidths=0,rasterized=True)
        ax.axvline(elbow,color=ORANGE,linestyle=':',linewidth=1.2,label=f'99.5% elbow ({elbow:,.2f})')
        ax.scatter(strict_points[distance],strict_points.weighted_mi,s=30,color=BLUE,edgecolors=INK,linewidths=.35,label=f'Real {distance_name} SHC strict (n={len(strict_points)})',zorder=3)
        if show_toy:
            passed=toy[toy.passes_shc]
            ax.scatter(passed[distance],passed.weighted_mi,s=38,color=ORANGE,edgecolors=INK,linewidths=.4,label=f'Injected SHC pass (n={len(passed)})',zorder=4)
        ax.set(xlabel=xlabel,ylabel='Weighted MI v2 (nats)',title=title); ax.legend(frameon=False,fontsize=7,loc='upper right'); style(ax)
    fig.savefig(outdir/f'{output_stem}.png',dpi=220)
    fig.savefig(outdir/f'{output_stem}.svg')
    plt.close(fig)


def main():
    outdir=WF/'output/figures'; outdir.mkdir(parents=True,exist_ok=True)
    base=pd.read_csv(ROOT/'Leca/results/gene_gene_wgmlst/weighted_mi_v2/kp_weighted_mi_v2_all_pairs.tsv.gz',sep='\t',compression='gzip',usecols=['gene1','gene2','mi','mean_order_dist','mean_bp_dist']).rename(columns={'mi':'weighted_mi'})
    base['pair_id']=['|'.join(sorted(x)) for x in zip(base.gene1,base.gene2)]
    strict_dir=ROOT/'Leca/results/gene_gene_wgmlst/spydrpick_locus_nmi/hiercc_shc_correction'
    elbow_dir=ROOT/'Leca/results/gene_gene_wgmlst/spydrpick_locus_nmi/mi_nmi_distance_comparison'
    strict_order=pd.read_csv(strict_dir/'order_dist_shc_strict_direct_edges.tsv',sep='\t')
    strict_bp=pd.read_csv(strict_dir/'bp_dist_shc_strict_direct_edges.tsv',sep='\t')
    order_elbow=float(pd.read_csv(elbow_dir/'order_dist_mi_aracne_elbow_summary.tsv',sep='\t').iloc[0].upper_q995_segmented_elbow)
    bp_elbow=float(pd.read_csv(elbow_dir/'bp_dist_mi_aracne_elbow_summary.tsv',sep='\t').iloc[0].upper_q995_segmented_elbow)
    toy=synthetic_pair_tests()
    plot_distance(base,strict_order,toy,'mean_order_dist',order_elbow,'Mean order distance','pairwise_order_distance_mi_shc_true_toy')
    plot_distance(base,strict_bp,toy,'mean_bp_dist',bp_elbow,'Mean physical distance (bp)','pairwise_bp_distance_mi_shc_true_toy')
    pd.DataFrame([
        {'selection':'order_dist','upper_q995_elbow':order_elbow,'strict_pairs':len(strict_order),'minimum_strict_distance':strict_order.mean_order_dist.min(),'strict_pairs_below_elbow':int((strict_order.mean_order_dist<order_elbow).sum())},
        {'selection':'bp_dist','upper_q995_elbow':bp_elbow,'strict_pairs':len(strict_bp),'minimum_strict_distance':strict_bp.mean_bp_dist.min(),'strict_pairs_below_elbow':int((strict_bp.mean_bp_dist<bp_elbow).sum())},
    ]).to_csv(WF/'output/pairwise_selection_qa.tsv',sep='\t',index=False)


if __name__=='__main__': main()
