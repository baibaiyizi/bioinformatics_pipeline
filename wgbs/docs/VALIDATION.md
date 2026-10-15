> 下文为职责重构前的历史验证，不能代表当前源码已完成真实数据重算。职责重构时另通过源码/内嵌R解析及可搬运小型端到端检查；对应集中维护记录已按研究者要求清理，不能作为当前源码验收。

# WGBS 模板验证记录

日期：2026-09-29。仅修改独立 `model` 仓库；来源项目只读。未升级环境，未启动真实项目全量分析。`provenance.tsv` 记录提取来源及原文件 SHA256；`validation/executed_source_sha256.tsv` 冻结最后一轮实际执行的模板源码、参数和合成输入哈希。

## 已实际运行

| 层级 | 验证与结果 | 证据 |
|---|---|---|
| 语法与合同 | Shell、R/Rmd chunk、Python 语法通过；两模板共 42 项合同/调度检查通过。覆盖空字段、重复 sample、生物学单位重复、配对字段、错误 assembly、缺组、样本重排、08/09、UNVERIFIED/STALE。 | 本文及 `validation/syntax_and_contract.json` |
| 独立复制 | 新目录含空格与单引号；正式入口读取相对配置和模板内辅助脚本成功。help/list/check/dry-run 不创建结果。 | 最终目录 `/tmp/final-wgbsseq numeric's m9433f9u/project` |
| 合成 CpG 核心 | 4 个独立合成样本，各 300 CpG；正式 `run_upstream.sh --steps 12,13` 构建 methylKit/CpG.RData 和覆盖汇总，`run_wgbs.sh --steps 1,3` 运行 DSS、区域聚合、注释、图表与 HTML。 | `examples/synthetic_core/` 可重放；`validation/synthetic_core_execution.log` |
| 差异与方向 | 300 DML 均保留原始 pval、DSS fdr 及额外位点 BH 列 padj_BH；得到 1 个合成 DMR，treated − control = +0.60；含显式组标签。 | `validation/dml_summary.csv`、`dmr_raw.csv`、`dmr_annotated.csv` |
| 报告状态 | 有图核心报告成功；分别构造仅表 COMPLETED、SKIPPED_DISABLED、FAILED、STARTED。后 3 种显示状态并排除旧表。 | `validation/report_states.json` |
| 方向/缓存负例 | RNA 正向元数据通过；缺失、NA、反向拒绝。DMR 原始方向标签不符/缺失拒绝。配置反转使旧缓存 STALE，无签名旧结果 UNVERIFIED。 | 合同与 helper 的隔离检查 |

上述是合成小数据执行，不是原始实验复现。示例 GTF 的坐标与 transcript 人为构造，只检查解析与注释接口。DSS DMR 使用其区域聚合规则；位点 BH 不代表区域 FDR。

## 未运行与适用边界

真实 FASTQ QC、Bismark、MethylDackel、BAM 复用的大规模运行、真实 GRCm39 参考注释、富集/调控/轨迹、精子专题和跨组学：`not_run`。dmrseq 在本次执行为 `SKIPPED_DISABLED`。原始/clean 输入分支只完成源码与语法审计，不能由合成 cached 模式推断已验证。小鼠 GRCm39 独立两组以外设计未实现；技术库/混池内个体不扩增独立 n。

版本见 `validation/dependency_versions.tsv`。实际 R 4.5.1、methylKit 1.36.0、DSS 2.58.0、bsseq 1.46.0、txdbmaker 1.6.2；现代 Bioconductor 的 GTF 建库已改用 txdbmaker。既有 conda 环境定义可参考工作区 `server_migration/envs/wgbs.yml`；模板运行不依赖该路径、不自动更新环境。

DSS 原生方向字段保留：DML 的 `diff`/`stat` 和 DMR 的 `areaStat` 沿用 DSS 的 group1（reference）− group2（comparison）；模板规范化的 `diffMethy`/`diff.Methy` 为 comparison − reference，并附 `group1_label`/`group2_label`。因此同一合成升高区域可以有原生 diff/stat/areaStat 为负、规范化 diffMethy 为正。筛选、Hyper/Hypo 和报告效应统一读取规范化字段，不能把原生统计量的符号混作规范化效应方向。原始 p 和 BH/FDR 不因换向而改变。
