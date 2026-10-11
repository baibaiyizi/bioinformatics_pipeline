# smallRNA 整理与验证记录

验证日期：2026-10-11。使用现有环境，未安装或升级依赖；最终生产源码在中心模板及三个课题模块中逐文件一致，共 33 份入口、预检、上游和 R/Rmd 源文件，见 [源码快照](SOURCE_SNAPSHOT.tsv)。源码冻结后不自动跟随模板更新。

## 真实项目的代表性验证

| 数据集与阶段 | 实际执行 | 边界 | 图形技术 QA |
|---|---|---|---:|
| plasma describe | qc、sncrna completed | 仅两模块；单库描述 | 42/42 PASS |
| MII differential | mirna completed | MOMII−MYMII，7 个库；保留原设计 | 8/8 PASS |
| MII enrichment | 入口、当前差异凭据衔接及状态报告 completed | **目标来源 skipped；GO/KEGG ORA 全部 not_run** | 6/6 PASS |
| serum differential | capture、capture_bias completed | 供体 2 配对技术分装，方法 2−方法 1 | 16/16 PASS |

运行退出码、逐模块身份和当前输入/源码签名已核对，见 [real_execution.tsv](validation/real_execution.tsv)。报告与状态中的项目名称读取显式配置，不取公共模块目录名 smallrnaseq。上述 72 组 PDF/PNG 的技术 QA 通过；技术通过不替代所有图的科学语义审阅。

MII 的原始 miRNA 差异表及 serum 的 6 份主要捕获统计表与搬迁前结果逐字节一致，见 [科学表对照](validation/scientific_table_comparison.tsv)。MII 的差异凭据门实际验证：当前凭据通过，failed 状态、MD5 变化和输入缺失均被拦截，见 [凭据门记录](validation/stage_receipt_guard.txt)。

**富集限制：multiMiR 1.32.0 已安装，但 `.onLoad` 访问 multimir.org 时 DNS 失败。** 原目标加载函数将命名空间不可用概括成 package unavailable；本次不能声称复现历史 169 个目标基因或完成 ORA。目标来源与各组 GO BP/KEGG 状态均保留 skipped/not_run；没有切换替代来源。详见 [命名空间说明](validation/R_namespace_notes.md)。

## 模板与可搬动性

- Mus/Homo × RUN_ENDOSIRNA=0/1 四个独立合成案例均实际完成 sncrna；启用 endo 时同时完成 endo。32 项检查通过，覆盖只读 help/list/check、缺 FASTQ 的上游拒绝、无效模块拒绝、显式配置覆盖继承环境、当前签名和共同分类。见 [检查表](validation/synthetic_portability_checks.tsv)。
- 合成数据均回收 65 reads；启用 endo 后只有符合条件的 30 reads 转类，20 reads 显式 miRNA 冲突保留原类；长度窗口和宿主匹配边界保留。四案例 108/108 组图技术 QA 通过，见 [图形汇总](validation/synthetic_figure_summary.tsv)。这是合成表输入验证，不是 FASTQ 上游复现。
- 同一鼠 endo 案例复制至第二个独立含空格绝对根，从不同 cwd 重跑后，类别计数、组成审计、转类审计逐字节一致，两份实际路径签名均有效。4 项通过，见 [跨根记录](validation/relocation_checks.tsv)。
- 六个真实数据集从 /tmp 和自身目录分别执行只读预检，12/12 通过；PN5 原始 FASTQ 缺失明确提示，已有 SPORTS 描述输入不受影响。见 [真实预检](validation/real_project_preflight.tsv)。
- Shell/Python 语法检查通过；30 份原 R/Rmd/upstream 文件已核对，变更限于资源定位、运行签名与显式项目标题，科学计算代码保留。见 [源码对照](validation/algorithm_source_comparison.tsv)。

复建合成验证使用 `python3 scripts/verify_portability.py --output <新的临时目录>`。测试在含空格和不同绝对根下执行，结果摘要保存在本目录 validation，完整临时案例位于 /tmp/smallrna-portability-20261011-accepted。

## 未执行与历史边界

本次未重跑任何原始 FASTQ 上游；未重跑全部真实默认 describe；PN5/E2C/8C 新源码仅预检，实际阶段运行保持 not_run。没有新增生物学重复、改变样本选择、比较方向或统计方法。

搬迁前 result 完整保存在各数据集 history/pre_layout_20261011；整理中不同冻结版本的代表性运行也按实际版本归入 history，原状态与签名均未改写。serum/report 原中文离线汇报包原样保留。当前 MII 富集的网络限制不能用历史完成状态补齐。
