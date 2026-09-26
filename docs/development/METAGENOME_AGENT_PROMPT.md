# Cinch 全量方法学交接与拔高任务 Prompt

你是一位同时熟悉统计遗传学、细菌群体基因组学、宏基因组学、微生物生态学、组成数据分析、信息论、因果推断、置换检验和计算方法开发的资深方法学家。请对下面的 Cinch 项目进行独立、尖锐、可证伪的方法学评价，并尝试把它提升为一套真正成立的“环境宏基因组多尺度 epistasis”理论与算法。不要迎合作者，不要因为框架复杂、数据量大或概念漂亮就默认它成立。必须区分数学恒等式、已实现算法、已验证结果、探索性发现、尚未实现构思和纯推测。凡是需要外部文献支持的内容，请检索并引用原始论文或权威方法文献，不得编造引用。

这是一个完全自包含的文字交接。你不会收到、也不需要索取原始脚本、结果文件、聊天记录或本地目录；以下内容就是本轮评价的完整项目事实底稿。文件名和组件名如果偶尔出现，只用于说明实现来源与可复现性缺口，不构成需要访问的附件。若某项细节在本 prompt 中没有给出，请明确标记为未知，不要假设可从附件补齐。

## 一、你必须完成的任务

1. 判断 Cinch 当前哪些部分具有真正的方法学原创性，哪些部分是已有 MI/NMI、JSD、SpydrPick、ARACNE、群体结构分层、条件置换、生态共现分析或 GWAS 思路的重新组合。
2. 判断“epistasis”这一术语在以下层次是否科学成立：物种级、菌株级、基因级、突变级、同尺度组合以及所有跨尺度组合。若不成立，请给出更准确的术语、可证伪定义和升级为严格 epistasis 所需的额外证据。
3. 把现有 Cinch 从细菌 isolate wg/cgMLST 流程形式化为环境宏基因组框架，明确 observation model、estimand、null model、score、筛选、推断、网络、实验验证和解释边界。
4. 重点解决：组成效应、测序深度/检出过程、零值与 non-call、物种归属、菌株归属、基因归属、突变归属、物理连锁、水平转移、谱系结构、环境异质性、批次效应、研究间异质性和多重检验。
5. 设计一个能处理超过 20,000 个公共环境宏基因组的可扩展算法，不允许直接穷举所有跨尺度 pair/triple 后再依赖经验阈值。
6. 设计模拟、benchmark、外部验证和实验室竞争实验，使结果能够区分“普通共现”“共同生境偏好”“分类或物理归属”“谱系标记”“潜在统计 epistasis”和“功能/因果 epistasis”。
7. 给出论文级理论叙事、核心图、必须证明的定理或性质、主要反例、最致命审稿意见、补救路线和合理投稿档次。
8. 最后给出一个比当前方案更强但仍可实施的 Cinch v2/v3 蓝图，不要只提出模糊方向。

## 二、证据状态标签

阅读以下材料时请严格使用这些标签：

- `[FROZEN-VALIDATED]`：代码、输入、结果和数值已经冻结并做过一致性验证。
- `[IMPLEMENTED]`：已有可运行代码和输出，但不一定已经纳入主流水线或完成外部验证。
- `[EXPLORATORY]`：已做探索分析；阈值或解释不能视为正式推断。
- `[PROPOSED]`：当前讨论形成的宏基因组扩展构思，尚未实现。
- `[UNRESOLVED]`：关键方法决定、可识别性或实现缺口。
- `[FORBIDDEN-CLAIM]`：现有证据不允许作出的结论。

## 三、项目定位与解释边界

`[FROZEN-VALIDATED]` Cinch 当前是一套面向细菌 wg/cgMLST profile 的分层统计候选筛选流程。它把原始关联、样本去冗余权重、MI/NMI、EpiDis、物理距离、系统发育/群体结构条件化、置换推断、网络剪枝、高阶背景差异和可视化拆成独立层。

`[FORBIDDEN-CLAIM]` 任意 retained edge、IQR outlier、SHC-passing edge、网络 community 或 Diff-GWES triple 都不是分子 epistasis、适合度互作或因果关系的直接证明。

`[FROZEN-VALIDATED]` 命名必须保持：论文公式的直接实现叫 `pairwise EpiDis`；无 pseudocount、Hamming-neighbour 加权的信息量叫 `weighted MI v2`；完整分层过程叫 `Cinch`。三者不能互换。

## 四、当前数据、输入契约与 pair universe

`[FROZEN-VALIDATED]` 主要 Kp 数据为 1,000 个基因组、41,312 个 profile loci。当前 profile 以 `value > 0` 表示 presence 或正 allele call，`value <= 0` 表示 absence/non-call。这个契约把真实缺失、不可调用、组装失败和生物学 absence 合并，是核心限制。

`[FROZEN-VALIDATED]` 冻结的 wgMLST pair universe 有 3,279 个 loci、2,609,249 个唯一有效 pair。cgMLST-like core 集合使用 prevalence ≥0.95，得到 3,632 loci；gene-level allele 编码把出现至少 20 次的 allele 保留为独立类别，其余稀有 allele 与 missing/non-call 进入 0 类。

`[FROZEN-VALIDATED]` 上游 `Cinch_v8.py` 对参考 CDS 去重、将 genome 映射到参考、构建 allele profile 和 pairwise 字段。wgMLST 分支使用二元 presence、raw correlation、PCA residualized correlation、tree distance 和 physical distance；cgMLST 分支计算 allele-state NMI 和分层 permutation/G-test 诊断。PCA residualization 只是候选宇宙构建，不是最终系统发育检验。

## 五、历史探索和可视化流程

主 runner 注册了 22 个 stage，其中 15 个 historical stage 和 7 个 MI/EpiDis/SHC extension stage：

1. `01-significant-ordination`：显著 r2 profile 的 PCA/UMAP，按 strength 和 degree 着色。
2. `02-distance-decay`：order、bp、phylogenetic distance 衰减、分箱趋势、elbow 和过滤候选。
3. `03-qq`：raw、adjusted 和 distance-filtered p-value QQ 图；历史结果显示明显 genomic inflation，提示结构未被简单调整消除。
4. `04-sign-aware-phylogeny`：正负关联分别进行 phylogeny filter 的诊断。
5. `05-order-bp-pca`：order-filter 与 bp-filter 的 gene association-profile PCA、K-means、Venn 和一致性比较。
6. `06-pcoa-metadata`：gene–gene profile PCoA；gene presence 与 metadata level 的 NMI 矩阵；gene–metadata PCA/PCoA。
7. `07-rda-arrows`：把 gene–metadata association 作为 RDA/envfit 风格方向投影到已有 PCA/PCoA。
8. `08-nmds-constrained-rda`：NMDS、真正的 constrained RDA、K-means、向量和诊断；必须区分事后投影和 constrained ordination。
9. `09-population-pca`：genome × gene presence PCA、scree/cumulative variance 和 metadata association。
10. `10-cluster-annotation`：PCA cluster 注释，CAZyme、virulence/effectors、mobile element、transporter 注释，以及 Phandango-like tree/presence heatmap。
11. `10b-cluster-geography`：选定 70-gene cluster 的 India/Taiwan presence enrichment、Fisher exact、BH-FDR、gene–geography NMI 和比较图。旧表逐字段复现成功；India 58/70、Taiwan 61/70 达到 q<0.05，但 Taiwan NMI 明显强于 India，因此不能称 India-driven。
12. `11-cgmlst-allelic-scan`：cgMLST allele-state pair association 原型、NMI、G-test/近似推断和 allele-pair detail。
13. `12-cgmlst-gene-nmi`：全 core gene–gene multistate NMI、近似 chi-square、距离衰减和 top flattened-genome links。
14. `13-gene-phylo-mixing`：每个 gene 的 allele diversity 与 same-allele carrier phylogenetic mixing endpoint。
15. `14-cgmlst-phylo-filtered`：排除任何涉及 lineage-strong gene 的 pair 后重新进行 cgMLST NMI。
16. `15-profile-similarity`：wg/cg similarity 和 sample weights。
17. `16-old-weighted-mi`：旧 pseudocount weighted MI/NMI comparator 和 ARACNE candidates。
18. `17-mi-distance-visuals`：MI/NMI 与 order/bp distance、elbows、outliers、direct/indirect 图。
19. `18-weighted-mi-v2`：thesis-style Hamming-neighbour weighted MI v2。
20. `19-pairwise-epidis`：论文公式 pairwise EpiDis。
21. `20-epidis-visuals`：EpiDis distribution、MI comparison、distance relationship 和 outlier 图。
22. `21-shc-correction`：HierCC-like SHC 构建、SHC conditional MI、within-SHC permutation、BH 和 ARACNE。

`[UNRESOLVED]` 当前主 DAG 并非完全自包含：stage 15 所声明的 snapshot 当前缺失，只有主工程之外的原始版本。因此不能把 22-stage registry 等同于已经再次端到端运行成功的完整产品。

## 六、样本相似度与 inverse-neighbour weighting

`[IMPLEMENTED]` 旧 comparator 计算：wg feature 集 prevalence ≥1%（8,651 loci），cg 集 prevalence ≥95%（3,632 loci）。wg gene-content Jaccard 中位数约 0.8194；cg allele concordance 中位数约 0.1403。CGAV similarity 使用 `(same+0.5)/(jointly callable+1)`，distance 为 `-log(similarity)`。

`[IMPLEMENTED]` 旧正式 NMI 使用 wg-content Jaccard ≥0.90 的邻居数 `m_i`，样本权重 `w_i=1/m_i`。当时 binary Hamming 在 0.9 similarity 阈值下的有效权重和仅约 7.0，Jaccard 权重约 101.2，因此旧流程选了 Jaccard。这是 Kp 实用决策，不是 EpiDis 论文公式。

`[FROZEN-VALIDATED]` weighted MI v2 和 pairwise EpiDis 使用 thesis-style genome-wide binary Hamming difference threshold `tau=0.10` 的 inverse-neighbour 权重：

```text
w_k_raw = 1 / count{l : mean_r Hamming(G_kr,G_lr) < tau}
w_k = w_k_raw / sum_s w_s_raw
```

self-neighbour 被包含，归一化后 `sum w_k=1`。需要在未来实现中把 distance definition、tau、self inclusion、normalization 和 effective sample size 全部写入 provenance。

## 七、weighted MI、NMI 与旧 comparator

`[IMPLEMENTED]` 旧 weighted MI 以二元 presence/absence 形成加权 2×2 表，每格加入 0.5 Jeffreys pseudocount。NMI 为 `2MI/(H(X)+H(Y))`。旧流程在 distance elbow 外按“每 locus 最大 non-LD NMI”的分布使用 `Q3+1.5 IQR` 和 `Q3+3 IQR`，再作 ARACNE。这些阈值不是 pairwise EpiDis 的全 pair IQR 阈值。

`[FROZEN-VALIDATED]` weighted MI v2 不加 pseudocount，使用加权联合概率：

```text
I_w(X;Y) = sum_xy p_w(x,y) ln[p_w(x,y)/(p_w(x)p_w(y))]
NMI_w = 2 I_w / [H_w(X)+H_w(Y)]
```

alpha/beta marginal proportions已经进入联合与边际概率，不能再乘一个 `4f(1-f)` 或所谓 frequency scalar，否则会重复加权。

`[FROZEN-VALIDATED]` 2,609,249 个 pair 全部唯一且字段完整；MI v2 范围约 `3.70e-14–0.693138` nats，NMI 范围约 `6.29e-14–1.0`，state frequency sum 最大误差约 `2.22e-16`。

`[IMPLEMENTED]` 旧 Jaccard/pseudocount NMI 在 3,279 loci、2,609,249 pairs 上曾得到 1,303 条 non-LD Tukey outlier、ARACNE 后 1,177 direct、126 indirect、9 extreme；weighted NMI 与旧 r2_adj Pearson≈0.802、Spearman≈0.201，说明大效应趋同但总体排序差异显著。

## 八、pairwise EpiDis 的公式、审计结论与阈值

`[FROZEN-VALIDATED]` EpiDis 按论文 equations 2-8 至 2-14 直接计算。对于 conditioning variable `X_i`：

```text
d_p = sum_{k:X_i=1} w_k
d_q = sum_{k:X_i=0} w_k
p1 = sum_{k:X_i=1} X_jk w_k / d_p
q1 = sum_{k:X_i=0} X_jk w_k / d_q
p=(p1,1-p1)+epsilon
q=(q1,1-q1)+epsilon
alpha=d_p/(d_p+d_q), beta=1-alpha
u=alpha p+beta q
JSD_w = alpha sum_s p_s log2(p_s/u_s)+beta sum_s q_s log2(q_s/u_s)
EpiDis=sqrt(JSD_w)
epsilon=1e-27
```

`[UNRESOLVED]` 论文没有说明加 epsilon 后是否重新归一化；冻结实现采用 literal no-renormalization，同时双方向计算并审计方向差异。单态 conditioning locus 会导致分母为零，必须过滤或标记 undefined，而不是设成 0。

`[FROZEN-VALIDATED]` 对有效归一化经验概率，核心恒等式是：

```text
EpiDis^2 = MI_bits = MI_nats/ln(2)
```

因此相同编码、pair universe 和 sample weights 下，pairwise EpiDis 与 weighted MI v2 排序完全相同；EpiDis 是 MI 的开方尺度而非新的依赖信息量。Kp 全量结果 Spearman≈1；`EpiDis^2` 对 `MI_bits` Pearson/Spearman≈1，最大误差约 `2.87e-15`。EpiDis 对 MI v2 Pearson≈0.912491；对旧 weighted MI/NMI 的 Spearman 仅约 0.360846/0.340807，因为旧权重和 pseudocount 改变了排序。

`[FROZEN-VALIDATED]` EpiDis 范围约 `2.31e-7–0.999993`，median≈0.117832，mean≈0.152389，最大双方向差约 `1.82e-10`。论文 pairwise outlier threshold 是在全 pair 分布上使用 `Q3+k×1.5IQR`，默认 `k=2`，即 `Q3+3IQR`；Kp threshold≈0.684741，得到 15,888 outliers，占约 0.6089%。非线性开方会改变 IQR fence，因此即使排序相同，MI 与 EpiDis 的 tail count 也可不同。

`[UNRESOLVED]` 公开 EpiDive repository 的 audited commit 没有实现 pairwise EpiDis、inverse-neighbour weights、IQR、CopSCC/Tarjan 或 differential-background workflow；论文公式是 score authority，不能声称复现了公开 repository 的完整核心算法。

## 九、物理距离、局部连锁与空间层

`[FROZEN-VALIDATED]` raw EpiDis/MI 高分尾部强烈富集于物理相邻 loci。例如 `ref_gene_692–693` order distance=1、bp distance=14；`692–694` order=2、bp=477。全局 IQR threshold 不是 linkage correction。

`[IMPLEMENTED]` 项目同时分析 order distance 与 flattened single-reference bp distance，使用分箱高分位趋势和 elbow。不同阶段曾产生不同操作阈值：早期 r2/order≈50.09、bp≈105,586；后期 MI 上尾 99.5% elbow 为 order≈102.85、bp≈246,474.5。阈值属于数据和统计量特定的操作点，不是生物常数；order distance 也不等于 kb。

`[UNRESOLVED]` accessory loci 只部分映射到统一参考：3,279 loci 中仅 952 有 single-reference flattened position。对未放置 loci，局部 LD/block preservation 不完整。宏基因组版必须区分同一 contig/replicon/mobile element 的物理连锁、taxonomic linkage 和跨宿主共现。

## 十、系统发育/群体结构校正：SHC

`[FROZEN-VALIDATED]` inverse-neighbour weighting 只降低过采样克隆的贡献，不消除 phylogenetic covariance。Cinch 另外用 3,632 cgMLST loci 的 allele concordance distance 建 single-linkage hierarchy，在连续层级之间用 partition NMI≥0.9 找 stable blocks，排除 singleton 主导层级，并在 stable block 内用 silhouette 选 SHC。

`[IMPLEMENTED]` Kp 选择 stable block HC129–HC430，最终 SHC=HC426，142 clusters，silhouette≈0.698。对每条 binary gene pair，在每个 SHC 中要求 cluster size≥5 且 X/Y 均多态；记录 informative SHC 和正向 odds-pattern SHC。至少两个 informative 且至少两个 positive SHC 才进入 permutation。

`[FROZEN-VALIDATED]` 条件 MI：

```text
I(X;Y|SHC)=sum_h P(h) I(X;Y|h)
```

在每个 eligible SHC 内独立 permute Y，999 permutations；`p=(1+#null>=observed)/(R+1)`，BH 后要求 q<0.05 且 z>0。之后才执行 ARACNE。

`[IMPLEMENTED]` 99.5% distance candidates：order 3,941、bp 3,793；至少两个 SHC 内重复正关联后为 108/90；SHC permutation q<0.05 后为 16/28；ARACNE strict direct 为 16/25。order 候选减少约 99.59%，bp 候选减少约 99.34%。多个看似远距离的强 pair 因没有任何 SHC 内变异被完全删除，说明原网络多数是 population background。

`[IMPLEMENTED]` 剩余主要涉及 `ref_gene_5282/5283/5285/6405/36159` 与 `ref_gene_992/993/995/996/999` 等两组；在约 7–11 个 SHC 重复，conditional MI≈0.08–0.135。它们仍可能是反复共同获得、共同选择或未解析模块，不能直接称分子 epistasis。

## 十一、ARACNE、网络与 PCA

`[IMPLEMENTED]` ARACNE 对任意闭合三角形删除 MI 唯一最弱边；若最小值并列则保留。它只是 data-processing/topology heuristic，不是显著性检验，也不保证恢复真实因果图。

`[EXPLORATORY]` SHC 后 25-edge 网络按 community 重画时 weighted modularity 仅 Q≈0.108；25 edges 中 10 intra-community、15 inter-community。四个 community 更像交叉关联中心而非四套独立 epistasis modules。community assignment 只用于阅读，不作为额外证据。

`[EXPLORATORY]` 把 3,279 loci 的 MI-neighbourhood 组成 locus×locus MI matrix 后做 PCA：raw MI PC1/PC2≈22.7%/16.7%，column-z-scored≈7.6%/5.3%。strict SHC 相关 loci 约 12 个；top MI profiles 常由局部相邻 gene block 驱动。PCA 几何更像强关联块/谱系块，不能包装成干净 epistasis modules。另有按 mean MI=`mi_strength/degree` 着色的 order/bp PCA；高 mean MI 并不与极端 PCA outlier 完全重合。

## 十二、cgMLST multistate NMI 与 layer3 探索

`[IMPLEMENTED]` cgMLST 把每个 gene 的 allele ID 当离散多分类变量。NMI 避免把 allele ID 当有序数值，但 allele 类别数、频率不均、稀疏 contingency table 和稀有 allele 合并会造成有限样本偏差。allele ID 只表示型别不同，不表达相差 1 SNP 还是 50 SNP。

`[EXPLORATORY]` allele-diversity图中曾用局部 y-density valley、trend residual K-means 及 12 种 clustering 方法识别三条纵向 layer。HDBSCAN 会把中层折叠为空，因此不适合预期的三轨迹目标。layer3 取 135 genes，共 9,045 pairs；旧 accessory-gene phylogenetic MPD threshold 对 core allele pairs 完全不适用，因为共同携带者近乎所有样本。

`[EXPLORATORY]` layer3 后改用 phylogeny-restricted allele permutation：按 tree-tip order 切成 200 blocks、每 block 5 genomes，在 block 内每 locus 独立 permute，500 次，保留局部 allele counts 和 missing rate。9,045 pairs 中 8,754 通过单侧 permutation-z BH q≤0.05；再加旧 order≥50.09 剩 6,510，网络过密。top 5% corrected-NMI backbone 为 326 edges、73 connected nodes。这里的 edge weight 是 observed NMI minus permutation-null mean NMI。

`[EXPLORATORY]` layer3 专属 order-decay elbow：raw NMI≈524、corrected NMI≈406；旧 50.09 位于明显衰减区，不能沿用。order≥406 后 raw NMI 分布单峰，峰约 0.40–0.45，无 0.5 自然谷；median≈0.4147，95%≈0.6028，99%≈0.7313。`NMI>0.5` 只是保留约上尾 17.3% 的人为强关联阈值。

## 十三、allele diversity–phylogenetic mixing 与 turnover baseline

`[FROZEN-VALIDATED]` 经典 endpoint 对每 locus 计算：

```text
Y = Simpson diversity = 1-sum_a p_a^2
X = pair-count-weighted mean patristic tree distance among genome pairs sharing an allele
```

只让 carrier count≥20 的 allele 贡献 X；core prevalence≥0.95。Kp 评估 3,632 loci；全 genome-pair mean tree distance≈0.0069329。早期按 X 下 10%/上 10% 标记 lineage-strong/mixed-high，各 364 loci；lineage threshold≈0.0003715，mixed threshold≈0.0047316。这些只是经验 quantile labels。随后 stage 14 排除所有涉及 lineage-strong locus 的 pair 重跑 cgMLST NMI。

`[IMPLEMENTED]` 后续提出理论 exponential allelic-turnover baseline，而不是只用经验 LOESS：给定某物种 pairwise tree-distance distribution D 和 turnover rate λ：

```text
H(lambda)=E[exp(-lambda D)]
SDI(lambda)=1-H(lambda)
MPD(lambda)=E[D exp(-lambda D)]/H(lambda)
```

对 observed locus 计算 pairwise SDI=`1-P(same allele)` 和 same-allele MPD，将其按相同 SDI 投影到理论曲线，定义：

```text
normalized residual = (observed MPD-model expected MPD)/mean pairwise D
```

positive residual 表示相同 allele 跨越的系统发育距离超过垂直 turnover baseline 预期。探索 threshold 为 species/profile 内 `median + 3×1.4826×MAD`；这是 robust outlier threshold，不是已校准显著性检验。

`[IMPLEMENTED]` original `0427_1.profile`：1,000 genomes、3,632 core loci、172 phylogeny-external candidates，占 4.74%，residual median≈-0.0425，threshold≈0.1987。

`[EXPLORATORY]` 在 1,663-genome Mycobacterium 数据上比较原始 CGAV、95% usearch 和 usearch-rm allele profiles，并按 top species/complex 分层。不同 allele-calling pipeline 的候选比例差异远大于 species 间差异：例如 M. abscessus 为约 0.05%、3.80%、21.93%；usearch-rm 下 439/2,002 loci 为候选。这个差异说明 allele definition 是主要 measurement model，21.93% 不能视为物种固有参数。

`[EXPLORATORY]` 用外部 HGT frequency 给 usearch-rm 三物种上色，HGT 与 normalized residual 的 Spearman 接近 0，候选与背景 HGT median 差异不显著；HGT score 在部分物种近饱和，当前指标不能解释上方离群层。

## 十四、SHC-conditioned Diff-GWES 高阶扩展

`[FROZEN-VALIDATED]` 高阶不是穷举任意 tuple，而是测试 background C 是否改变 seed pair A–B 的 SHC-conditioned association：

```text
I_c(A;B|SHC)=sum_h P(h|c, informative) I(A;B|h,c)
E_c=sqrt(I_c/ln(2))
Delta(A,B;C)=E_1-E_0
```

Delta>0 表示 C=1 背景增强，Delta<0 表示减弱。这是 Cinch/Leca 的明确扩展，不是论文高阶公式的逐字复刻。

`[FROZEN-VALIDATED]` 样本在每个 SHC 内 deterministic 50:50 split 为 discovery/validation。background prevalence 0.05–0.95；每个 validation background state≥25 samples；每个 state 至少 2 informative SHCs；每 SHC-state cell≥3 且 A/B 均变异。discovery 保留 `abs(Delta)≥0.15`，每 seed 最多 25 backgrounds。若 A/B/C 都有位置，排除与 A/B 共 100 kb reference block 的 C。

`[FROZEN-VALIDATED]` validation 中 C 在 SHC 内置换；同一 reference block 的 C 使用一致 permutation streams，A/B 固定。双侧 p=`(1+#|Delta_null|≥|Delta_obs|)/(999+1)`；对独立选择且 validation-eligible hypotheses 做 BH。

`[FROZEN-VALIDATED]` true Kp：63,644 discovery tests，351 validation-eligible，0 个 q<0.05，minimum q≈0.5383。toy：20/20 synthetic truth triples 显著，20/20 predeclared negative controls 均不显著，non-truth false positives=0；20 个 synthetic A–B 独立 SHC permutation 都 p=q=0.001。该 toy 只证明对其构造的强效应有 sensitivity/specificity，不证明对真实中弱效应有效。

## 十五、当前软件打包与可复现性

`[IMPLEMENTED]` 当前打包包含：命令行入口、22-stage workflow registry、主 runner、方法学说明、结果纪录、EpiDis 对照说明、历史覆盖审计、算法规格、聊天 JSON、独立的 SHC-conditioned Diff-GWES workflow，以及独立的 allele-diversity/phylogenetic-mixing workflow。这里仅报告组件状态，不要求你访问这些文件。

`[FROZEN-VALIDATED]` runner 支持 historical/extension/all tracks、`--from`、`--through`、`--only`、`--dry-run`、`--force`、独立 stage logs、declared outputs 和 `run_manifest.json`；缺 input、非零 exit 或缺 declared output 会停止。冻结文件使用 SHA-256 验证。

`[UNRESOLVED]` 生成 `0427_1.ldlike_pairs.tsv.gz` 的最初精确命令没有与 pair file 同地保存；`Cinch_v8.py` 还依赖非 pip 的 `configure` 与 `uberBlast`；Diff-GWES 是单独 frozen workflow，未完整并入 22-stage runner；turnover baseline、HGT、layer3 和 MI-PCA 也是另行实现，未全部进入 `workflow.json`。

## 十六、宏基因组多尺度 Cinch 构想

目标是把 epistasis 扩展到四个分辨率，并允许所有同尺度与跨尺度类型进入候选检测：

1. species level：species–species；
2. strain level：strain–strain；
3. gene level：gene–gene；
4. mutation level：variant/allele–variant/allele；
5. cross-scale：species–strain、species–gene、species–mutation、strain–gene、strain–mutation、gene–mutation，以及必要的 higher-order background modification。

`[PROPOSED]` “全部检测”应解释为所有 interaction type 都有资格进入，不是穷举所有组合。总变量数 p 时 pair 为 O(p²)、triple 为 O(p³)；20k 样本并不能消除组合爆炸和多重检验。

### 16.1 两级 epistasis 定义

`[PROPOSED]` 在没有 phenotype Y 的环境宏基因组中，先定义 `structure-corrected potential epistasis`：

```text
potential epistasis = observed dependence - dependence expected under a structure-aware null
```

structure-aware null 至少应考虑 taxonomic hierarchy、physical linkage、population/strain structure、study/batch、habitat、sampling design、composition 和 callability。它仍是统计候选。

`[PROPOSED]` 有明确生态/功能结果 Y 时定义 `functional epistasis`：联合 perturbation 对 Y 的作用偏离预先声明的 additive、multiplicative 或 mechanistic null。只有 factorial perturbation、方向预注册、阴性组合、必要时 rescue/reconstruction 才能把候选升级为功能 epistasis。

`[UNRESOLVED]` 如果拒绝两级定义，就必须提出一个不依赖 Y、且能把 epistasis 与 ordinary conditional co-occurrence 区分开的可证伪标准。

### 16.2 多尺度 observation model

同一批 metagenomic reads 同时被用来估计 species、strain、gene 和 mutation，所以四层不是独立数据源，而是嵌套且共享 measurement error。每层应输出：

```text
Z_sj = latent/estimated biological state
M_sj = callable/measurement confidence
U_sj = uncertainty or posterior distribution
```

不能直接把四层 raw relative abundance 拼成一个 MI matrix。

建议分别定义分母和状态：

- species：相对全 microbial community，理想情况下有 absolute cell/genome abundance anchor；
- strain：相对所属 species 的 strain mixture；
- gene：相对可信 host species、genome equivalent 或 single-copy marker 的 gene copy/presence；
- mutation：相对对应 callable locus coverage 的 alternate/reference counts 或 allele posterior。

gene–species、mutation–gene 等确定性的层级归属不能被当作 interaction。可靠 linkage 工作可以降低问题，但 mobile/shared genes、多宿主 gene、strain mixture、short-read ambiguity 和低覆盖区应保留 probabilistic assignment。

### 16.3 测序深度与抽样效应

测序深度可视为生态抽样/检出过程，但不应默认校正后无影响。若 sample depth 为 N、真实相对丰度 q，则近似 detection probability：

```text
P(detect)=1-(1-q)^N
```

两个稀有 feature 在浅测序样本中共同为 0、深测序样本中共同为 1，可制造正 MI。mutation 的 zero alternate read 也可能是 reference state 或 coverage 不足。必须区分 present/absent/non-callable，设置 feature-specific callable masks、minimum coverage、effective sample size，并用 downsampling/simulation 检查 score 方向和排序稳定性。

### 16.4 组成效应与两个 estimand

相对丰度满足 closure：

```text
q_sj=x_sj/sum_k x_sk,  sum_j q_sj=1
```

一个 feature 绝对增加会机械性压低其他 feature 的相对比例。因此必须区分：

- `absolute-Cinch`：推断 cell/genome/gene/molecule 的绝对数量 interaction；通常需要 spike-in、qPCR、flow cytometry、总 microbial load 或强识别假设。
- `compositional-Cinch`：推断 log-ratio/balance 空间中的 interaction；结果只能解释为相对 balance，不可重新解释为两个对象的绝对协同增长。

`[PROPOSED]` 最科学的方案可能是两个不可互换分支，并研究何时结论一致、何时相反。需要评价 CLR/ALR/ILR、reference frames、multinomial/Dirichlet-multinomial/count models、absolute abundance anchors、zero handling 与 MI/conditional MI 的兼容性。

### 16.5 habitat 与跨生境 interaction

单纯比较“加入/去除 habitat covariate 前后”只是 sensitivity analysis，不能自动识别 habitat-specific epistasis。建议直接估计：

```text
I(X;Y | H=h)
Delta_h1,h2 = I(X;Y|H=h1)-I(X;Y|H=h2)
```

habitat-specific candidate 应要求该 habitat 内 effect、跨 habitat heterogeneity test 和 structure-preserving null 同时通过；cross-habitat candidate 应要求方向一致、效应稳定，并通过 leave-one-habitat/study validation。若环境变量既是 confounder 又是 effect modifier，应明确 estimand，而不是机械“校正掉”。

### 16.6 计划中的大规模研究设计

`[PROPOSED]` 使用 >20,000 个公共环境宏基因组进行 discovery，寻找 habitat-related multiscale epi combinations，然后用实验室生态竞争实验验证。

不能 random sample split；应至少包括：study-level deduplication、protocol/batch harmonization、leave-one-study-out、leave-one-habitat-out、完全冻结的 external test、negative-control habitats、negative-control feature pairs、study-aware permutation 或 hierarchical model。公共样本数量不等于有效独立样本量。

实验必须在数据分析完成前冻结候选、效应方向和 null。最低设计是 factorial：neither/A only/B only/A+B；物种/菌株层做 pairwise 或 synthetic-community competition，gene/mutation 层做 isogenic knockout/knock-in/allele swap，跨尺度最好形成 mutation→gene function→strain fitness→species coexistence→community outcome 的可验证链。单次 co-culture 中“竞争成功”只验证一个生态现象，不足以证明多尺度理论。

## 十七、下一位 agent 必须重点攻击的问题

1. 没有 phenotype Y 时，把 residual conditional dependence 称为 epistasis 是否只是改名？
2. `potential epistasis` 的 null 应是 independence、conditional independence、maximum entropy、log-linear no-interaction、生态 neutral model，还是生成式 measurement+ecology model？
3. 二元 EpiDis 已等价于 sqrt(MI_bits)，多尺度理论的新信息量究竟在哪里？
4. 多分类、多变量和连续 abundance 情况下应使用 MI、conditional MI、interaction information、co-information、partial information decomposition、total correlation、log-linear interaction，还是其他量？这些量对 synergy/redundancy 的符号和解释是否稳定？
5. 跨尺度变量存在确定性嵌套时，如何定义不是 taxonomy identity 的 interaction？是否应使用 conditional MI，例如 gene–habitat association conditional on host species/strain，或 residualized latent-state graph？
6. sample weights 用全 feature matrix 构造时是否会把待检验 feature 信息泄漏进权重？是否需要 leave-pair-out distance、cross-fitting 或独立 marker set？
7. SHC within-cluster permutation 的交换性条件何时成立？公共 metagenome 中 study、space、time、host、site 重复测量如何保留？
8. 组成数据下 MI 是否应在 log-ratio latent abundance 上定义？离散化会不会丢失或制造 interaction？
9. 稀疏多分类 contingency table 的 bias correction、effective sample size 和 uncertainty 如何处理？
10. 物理 linkage、taxonomic linkage、mobile-element linkage、shared plasmid 和 horizontal transfer 应是过滤、协变量、分层还是竞争 null？
11. ARACNE 的 data-processing assumption 在有 synergy、共同原因和跨尺度确定性节点时是否成立？
12. 高阶 interaction 应继续采用 background-dependent pair score，还是需要真正的 multivariate synergy 定义？
13. empirical elbow/IQR/MAD threshold 哪些只能做 discovery，哪些能转换为 calibrated FDR/error-control？
14. 如何避免 20k public datasets 主要学习到 protocol、geography 和 database annotation artifacts？
15. 什么结果能证伪 Cinch，而不是所有结果都能被解释为“结构复杂”？

## 十八、你应交付的输出

请按以下顺序输出，不要跳过负面结论：

### A. Executive verdict

- 用一段话判断：Cinch 当前是原创理论、原创组合方法、工程整合，还是探索工作流。
- 分别给出 mathematical novelty、statistical validity、biological interpretability、computational feasibility、experimental testability 的 0–10 分和置信度。
- 明确给出最可能的 desk-reject 原因。

### B. Claim audit

建立表格：当前 claim、证据状态、是否成立、最强反例、允许的改写。必须覆盖 EpiDis novelty、phylogeny correction、high-order epistasis、turnover/HGT interpretation 和 metagenome universality。

### C. Formal theory

提出符号系统，至少包含 sample/study/habitat、四尺度 latent variables、callable masks、taxonomic parent map、physical linkage map、sample weights、composition/absolute-abundance branch、nuisance structure、pair score、high-order/background score和 phenotype Y。给出清楚的 estimand 与 null hypothesis。

### D. Cinch v2/v3 algorithm

给出可执行的 staged DAG、每阶段输入输出、复杂度、并行化/筛选策略、错误控制、停止条件和 provenance。必须说明怎样把 O(p²)/O(p³) 降到可运行规模，同时避免先验筛选导致 double dipping。

### E. Simulation and benchmark

设计包含下列机制的模拟：真实 interaction、无 interaction、composition closure、depth variation、zero inflation、study/batch、habitat confounding、taxonomic nesting、physical linkage、HGT/plasmid、population structure、strain mixture、misassignment 和 missing-not-at-random。提供应比较的方法、评价指标、预期失败模式和 calibration 标准。

### F. >20k public-metagenome study

设计 discovery/replication/external validation；给出 metadata 最低要求、数据 QC、study weighting、habitat ontology、sample independence、held-out policy、negative controls、effect heterogeneity 和 reproducibility package。

### G. Experimental validation

给出 species/strain/gene/mutation 及跨尺度各至少一个实验设计；区分 association confirmation、fitness interaction、mechanistic validation。要求 factorial contrasts、effect-size definition、replication、power、negative controls、rescue 和预注册式候选冻结。

### H. Paper elevation

给出最强但不过度的中心命题、论文标题备选、Figure 1–6 故事线、必须补齐的证据和合理投稿梯度。判断什么条件下适合 The ISME Journal/Nature Communications，什么条件下才配 Nature Microbiology/Nature Ecology & Evolution/Nature Methods，什么条件下仍不足以冲击这些期刊。

### I. Red-team verdict

以最苛刻审稿人身份列出 10 个 major concerns，并标记 `fatal before submission`、`must fix` 或 `can acknowledge`。最后明确回答：这个方向是否值得投入，若值得，最小可发表版本和最高价值版本分别是什么。

## 十九、不可违反的结论纪律

- 不得把 statistical association、outlier 或 community 直接称为 molecular/causal epistasis。
- 不得把 `EpiDis=sqrt(MI_bits)` 包装成新的 pairwise information quantity。
- 不得把 inverse-neighbour weighting 称为已消除 phylogeny。
- 不得把加入/去除 habitat covariate 称为已识别 habitat-specific interaction。
- 不得把 relative abundance 的 interaction 解释成 absolute abundance interaction。
- 不得把同一 species/strain/gene 的确定性归属关系当成跨尺度 epistasis。
- 不得因样本量 >20,000 就忽略 study-level dependence、batch 和多重检验。
- 不得因 toy benchmark 成功就宣称真实生态信号已验证。
- 不得因实验中某个 pair 竞争结果显著就宣称整个多尺度理论成立。
- 若一个核心主张不可识别或没有可证伪 null，请直接说不可成立，并给出替代定义。

请先完成独立评价，再提出拔高方案；不要先替作者补故事后再倒推合理性。
