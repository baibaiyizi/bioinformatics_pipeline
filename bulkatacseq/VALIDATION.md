# ATAC 模板验证记录

日期：2026-09-29。真实 liver 项目只读，未生产重算或升级环境。`provenance.tsv` 记录提取来源与源文件 SHA256；新 counts 模式记录为原区间 DESeq2 的新增实现。`validation/executed_source_sha256.tsv` 冻结最终实际运行的源码、配置与合成输入。

## 已实际运行

| 层级 | 验证与结果 | 证据 |
|---|---|---|
| 语法与合同 | Shell、R/Rmd chunk、Python 语法通过；两模板共 42 项样本、TSV 空字段、重复身份、独立单位、组方向、08/09、只读 CLI 与 STALE/UNVERIFIED 检查通过。 | `validation/syntax_and_contract.json` |
| BAM 入口轻量审计 | 复制到含空格/单引号目录，正式 `run_rmd.sh --steps 02` 在显式 ATAC_LIGHT_SETUP=TRUE 下生成输入清单及报告。 | `/tmp/epigenomics current's smoke 56ty2y90/bulkatacseq`；此项不证明 BAM 分析已运行 |
| counts 正式闭环 | 500 个合成整数 consensus peak，2+2 独立样本；生产入口运行实际 DESeq2 Wald、BH、全样本 PCA、火山图、CSV 与 HTML。500 峰检验，72 峰满足默认原始 p 主筛选。 | `examples/consensus_counts/` 可重放；`validation/counts_execution.log` 与 `counts_result_summary.json` |
| 身份与方向 | 最终样本 ID 01–04 的前导零保留；前 30 峰 Fold>1，31–60 峰 Fold<−1。另一次实际反向配置并重排样本后，效应全部翻转，72 个主集合峰不变。分数 counts 拒绝。 | 最终目录 `/tmp/final-bulkatacseq numeric's zhyui1c3/project`；反向目录 `/tmp/atac-counts-direction-11fwbiju/project`；`validation/counts_result_preview.csv` |
| 筛选口径 | p=0.01/FDR=0.2 进入原始 p 主集合；p 缺失时不借用 FDR；负 FDR 不视为显著。BINDetect 已知方向通过，未知条件/未声明样本失败。 | helper 隔离执行检查 |
| 报告与执行状态 | 图文报告及仅表 COMPLETED 成功；SKIPPED_DISABLED/FAILED/STARTED 报告显示状态并过滤旧表。独立审查分别确认 counts --report-only 对最新 FAILED/STARTED 拒绝且旧 index SHA256 不变。 | `validation/report_states.json`、`independent_report_review.json` |

版本见 `validation/dependency_versions.tsv`：R 4.5.1、DESeq2 1.50.2、htmltools 0.5.9、digest 0.6.39；服务器已安装 DiffBind 3.20.0，但本次未运行原生 BAM DiffBind。预检查合格与统计执行完成分别记录。报告生成失败只改变 report 状态，不回写已完成的 05 为 FAILED。

## 未运行与适用边界

FASTQ → BAM → peak/consensus counts、真实 BAM DiffBind、真实参考注释、富集、轨迹、IDR、footprinting/BINDetect、RNA/WGBS 整合、保守性：`not_run`。源码迁移与语法通过不代表这些依赖/参考已经在新项目验证。14 机器学习仍停用。IDR 保留每组恰好 2 个重复的实际限制；不选取子集凑数。

counts 模式仅在原 consensus 区间拟合 DESeq2，不执行 DiffBind summit 重定位，不声称与原生 BAM 检验对象相同。500 峰示例是人工统计夹具；不能据其显著峰数判断真实样本功效、质量或机制。沿用已有工具环境；工作区 `server_migration/envs/macs2.yml`、`idr.yml` 可供环境重建参考，模板运行不读取这些外部路径。
