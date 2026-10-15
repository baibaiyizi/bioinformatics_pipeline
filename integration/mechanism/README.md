# 机制程序与转录因子活性最小模板

本目录只收录两个输入语义清楚、已经验证的生产器。模板必须复制到项目本地并冻结源码 hash 后运行；中心模板不直接读取项目结果，也不代表任意项目可即插即用。

## 1. 跨情境通路程序周转

`../scripts/mechanism_programs.py`消费已经完成的完整 DE 表、完整 GSEA 表和冻结的基因集 RDS，不重新计算 DE 或 GSEA。每个通路在每对情境中的重叠检验背景，是该通路在两个情境中都具有有限 DE `stat` 的成员，而不是全基因组。所有预设“通路 × 情境对”同时进入一次 BH 校正；输出保留全部通路，不按 P 值删行。

输入契约：

- 每张 DE CSV 必须含 `gene_id/SYMBOL/log2FoldChange/stat/pvalue/padj`，并使用同一物种、assembly 和同一效应方向；
- GSEA TSV 必须含 `pathway/pval/padj_context_collection/NES/size/leadingEdge/context/collection`，其中 `leadingEdge` 用分号分隔；
- 基因集 RDS 必须含命名结构 `bundle$libraries[[collection]]`，每个元素为一个通路的基因符号向量；
- `--context NAME=DE_TABLE` 至少提供两次，名称必须与 GSEA 的 `context` 完全一致；
- `--assembly` 和 `--comparison-label` 必须由项目配置读取后显式传入，脚本不会从目录名猜测。

```bash
python3 scripts/mechanism_programs.py \
  --fgsea result/full_fgsea.tsv \
  --gene-sets-rds reference/frozen_gene_sets.rds \
  --context stage1=result/stage1_de.csv \
  --context stage2=result/stage2_de.csv \
  --collection Hallmark \
  --comparison-label "treated - control" \
  --assembly GRCm39 \
  --expected-pathways 50 \
  --output result/mechanism_programs
```

主要产物为四张完整源表和 `receipt.json`：情境层通路效应、通路全部成员效应、两两 leading-edge 周转，以及周转成员明细。Jaccard、NES 方向一致和成员效应方向一致只描述独立情境间的相似性，不能写成纵向持续、遗传传递或因果机制。

运行依赖为 Python 3、`pandas`、`numpy`、`scipy`、`statsmodels` 和可调用的 `Rscript`。`validation.json`记录了本模板用 PFOS 完整真实输入在隔离临时目录重放所得的命令、输入/代码/产物 hash、行数、300 个预设检验的 BH 分母和通路内零背景；临时结果表已删除。

## 2. 小鼠全队列 regulon 活性

`../scripts/run_full_cohort_regulon.R`用本地安装的 `dorothea_mm` A/B/C 有符号先验和 `decoupleR::run_ulm`，从完成 QC 的归一化 bulk RNA 表达矩阵重算每个文库的 TF 活性，再用 limma 比较指定两组。先验从本地 R 包读取并随输入一起记录 hash，不自动更新网络。

输入契约：

- DE CSV 必须含 `gene_id/SYMBOL/baseMean/stat/log2FoldChange/pvalue/padj`；
- 表达 TSV 第一列为与 DE 可匹配的 gene ID，其余列为非负、有限的归一化计数；
- 样本 TSV 必须含唯一的 `sample` 和 `group`，可另含 `biological_unit/sire_id/litter_id`；
- 样本表必须且只能包含显式指定的两组，每组至少两个文库；
- 当前资格验证只覆盖 `--species mouse`，assembly 必须显式传入。

```bash
Rscript scripts/run_full_cohort_regulon.R \
  --de result/full_cohort/DE_results.csv \
  --expression result/full_cohort/normalized_expression.tsv \
  --samples result/full_cohort/samples.tsv \
  --reference-group WT \
  --comparison-group PFOS \
  --species mouse \
  --assembly GRCm39 \
  --output result/full_cohort/regulon
```

输出包括逐文库 TF 活性、两组 TF 活性效应及 BH、全背景 target 效应、冻结的 DoRothEA A/B/C 网络、输入清单和 `receipt.json`。TF 活性来自 RNA 与先验网络的推断，不是 TF 结合或因果调控证据。若 `biological_unit` 或家系信息为 `unknown`，组间比较只能按文库层探索解释。`validation.json`保留 10 文库小鼠项目真实输入只读重放的最小凭据；项目真实结果表未复制到模板库。

运行依赖为 R、`decoupleR`、`dorothea`、`limma`、`jsonlite` 和 `digest`。

## 调控桥接部署资格

PFOS 项目的 `mechanism_regulation.py` 已完成源码审计，但状态为 `source_audited_not_deployed`。现有实现固定了 PFOS 的 ATAC/WGBS 表结构、4 个 ATAC 文库、107 个 chromVAR motif、10 个优先 motif hit 和 DMR 区域汇总覆盖语义。直接复制会把项目特定注释误装成通用调控验证。

只有新项目同时明确以下输入时，才应另建项目本地桥接脚本：样本/个体跨模态映射、共同 assembly、完整 chromVAR motif 结果、实测 motif hit、peak 到 target 的证据类型、WGBS 覆盖是同一 CpG 还是区域汇总，以及 motif 到 TF 的歧义处理规则。家族 motif 不强制映射为单 TF，nearest-gene 不作为已验证 target，未核验的 TOBIAS 不纳入，也不把依赖的多层 P 值合并。

`template_snapshot.tsv`冻结本目录说明、验证凭据和两个已部署脚本。复制到项目后，应由项目自己的快照记录实际采用版本和 hash。
