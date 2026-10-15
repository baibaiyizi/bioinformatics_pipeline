# 多时期RNA：真实扰动参照与同坐标剪接事件


本模板包含两个独立入口，均复制到项目本地运行。一个把公开真实干预的基因响应投影到本项目各时期的完整差异排序；另一个把同坐标剪接事件跨期对齐。它们回答不同问题，不合成“机制总分”。

## 真实扰动签名

```bash
Rscript scripts/perturbation_signatures.R \
  --ranked stage_MO_minus_MY.tsv --vst vst_matrix.tsv --samples samples.tsv \
  --obox-signatures obox_signatures.tsv --pabpn-effects pabpn1_effects.tsv \
  --output result/perturbation_signatures
```

`ranked`必须保留各时期全部可排序基因及`gene_id/gene_symbol/stat/log2FC/padj`；`samples`给出`sample_id/stage/age_group`。签名包含`signature_id/source_stage/direction/gene_symbol`。默认以OBOX作者定义的4个有方向集合为参照，在六时期构成24项计划检验；映射小于10个基因的集合标为NT，不能写成阴性。BH校正保留全部计划检验数。大集合允许至5000成员；实际映射及排除原因留表。

PABPN1输入为公开FPKM组均值构成的效应量，主分析伪计数1，敏感性伪计数0.1，并要求至少一组均值≥1。它只与阶段匹配的项目效应描述性比较，不从外部处理矩阵编造差异P值或FDR。公开数据准备脚本位于`scripts/prepare_perturbation_references.py`；下载与解析留存原始文件、URL及hash。

输出包括完整富集检验、全部leading edge、明确选取的展示成员、逐文库表达、PABPN1效应一致性、英文PDF/PNG及输入签名。图上的一致性表示表达响应相似，不能据此推断本项目发生了相同的基因敲除或同一分子机制。

## 同坐标剪接事件

当前模板把不同阶段的 rMATS 结果按完整事件坐标对齐，并输出逐文库 PSI、原始阶段内 P/FDR、读段支持和可选的真实 Sashimi 图。它不重新拟合 rMATS，也不制造年龄×时期交互显著性。

复制 `scripts/splicing_trajectories.R` 与两张配置表到项目本地后运行：

```bash
Rscript scripts/splicing_trajectories.R \
  --config-dir config/splicing \
  --rnaseq-root /path/to/stage/rnaseq \
  --output result/splicing_trajectories
```

`analysis.tsv` 定义时期顺序、统一比较方向、JC/JCEC范围和 Sashimi 容器；`events.tsv` 为每个事件提供类型、完整坐标及代表性 Sashimi 时期。项目必须从各时期 `run_metadata.tsv` 恢复 b1/b2，不能按文件名推断方向。

主要产物为：

- `event_stage_summary.tsv`：所有阶段、事件和计数模式的标准长表；
- `cross_stage_event_recurrence.tsv`：按完整坐标形成的跨期复现表；
- `target_event_psi_long.tsv` 与 `target_event_read_support.tsv`：保留 NA 位置的逐文库数据；
- `target_event_stage_statistics.tsv`：保留原始精确 P/FDR 和统一方向的 ΔPSI；
- 英文 PSI PDF/PNG，以及在 BAM、BAI、contig 与执行器均满足时生成的真实 Sashimi PDF/PNG。

图中只显示 `* / ** / *** / ns`，精确数值留在源表。Sashimi 是位点读段支持图，不是额外统计检验。`examples/splicing/run_example.sh` 验证方向反转、完整坐标防碰撞、PSI NA 保留以及缺少 BAM 时的诚实状态。


## 复用与验收

统一合成示例入口为`bash run_examples.sh`，路径按脚本位置解析，启动目录不影响输入。

- [扰动示例](examples/perturbation/)与[剪接示例](examples/splicing/)使用明确标注的合成输入，验证计算接口和边界行为；不代表上游论文被独立复现。
- 实际项目运行必须另核对输入、物种、比较方向、库顺序和来源，记录自己的运行凭据与源表。
- `template_snapshot.tsv`冻结本模板的源码；项目的同名快照记录实际复制的文件及hash。项目运行不引用中心模板实时源码。
- 图内使用英文，星号对应明确的原检验/校正范围，精确数值保留在表中；未检验和未显著分别显示。
