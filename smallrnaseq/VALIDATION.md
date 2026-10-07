# 验证记录

验证日期：2026-09-29。仅在模板和隔离的临时目录执行，真实 MII/serum 来源只读；未重算真实项目、未安装/升级依赖。

## 已执行验证

1. **独立复制后最小闭环通过**：由 scripts/make_fixture.py 创建 3+3 个独立合成样本、18 份人工 SPORTS summary/length/output 表；运行模块 1，输出 29 个表达矩阵文件、互斥主类别分母审计、组成 PDF/PNG、中文 HTML 报告、输入 SHA-256 与 sessionInfo。6 个文库互斥类别均回收 100%，ReadDifference=0。未使用真实项目输出制造 fixture。
2. **原模块 4 实际通过**：同一合成输入完成 setup→全局 sncRNA 描述→报告，产生 35 个 PDF；默认 RUN_INFERENCE=0，状态和跳过原因保留。此项验证不等于全部 9 个模块均已运行。
3. **真实 edgeR 函数执行通过**：3+3 个合成独立样本、50 个整数计数特征，植入 Gene1 上调；logFC=4.5579988189，P=1.0123989827e-33，FDR=5.0619949136e-32，SignificanceMethod=PValue。反向对比和逆序样本列仍保持正确方向。默认关闭推断、重复 biological_sample_id、缺 batch、n=1 的拦截及非整数计数拒绝均通过。
4. **只读入口通过**：run_rmd --help/--list/--dry-run/--check 与 sports --help/--dry-run 未新增目录；sports --check 在缺真实 FASTQ 时明确 exit=1，同样未写结果。08/09 解析通过。空 age_group 的 TSV 字段在上游解析中保持为空，species/stage 未错列。
5. **可选状态通过**：显式选择模块 9 而 RUN_ENDOSIRNA=0 时记录 disabled，并成功输出报告；没有默认猜测参考路径。
6. **仅表报告通过**：REPORT_FIGURES=0 的独立新 run 完成 setup→仅表 HTML，未生成组成图；模块状态为 completed。
7. **独立 serum 示例通过**：一对同供体合成工艺文库实际生成 spike RPM、配对观察差与中文边界说明，无分子回收率/人群推断。
8. **语法通过**：Shell 入口与 load_config 的 bash -n、Python 编译、所有 9 个 R/Rmd 的解析（Rmd 在临时文件 purl 后 parse）、报告 R 入口解析通过。未新增 R 测试框架。

环境为现有 PATH R 4.5.1；R 包精确版本见 validation/R_packages.tsv。主要验证目录 `/tmp/smallrna-template-validation`；小型证据摘要保存在 validation，核心输入/输出仍可按 README 的生成器重建。

## 状态边界

| 路径 | 状态 |
|---|---|
| 配置→SPORTS 表输入→模块 1→描述报告 | 合成端到端已运行 |
| 模块 4 全局 sncRNA 描述 | 合成输入已运行 |
| 可选 edgeR 入口 | 独立小型植入效应与边界检查已运行 |
| 模块 2/3/5/6/7/8 | 源码迁移/解析通过，执行 `not_run` |
| 原始 cutadapt/SPORTS/Bowtie 上游 | `not_run`，真实工具与参考由配置提供 |
| endo-siRNA 启用分支 | `not_run`，默认 disabled |
| 在线 miRNA 靶标/富集 | `not_run`，默认 disabled |
| MII 专题文献片段 | 参考档案，`not_run` |
| 真实数据复现、独立生物学重复验证 | `not_run` |

Rmd 执行入口设置 error=FALSE；错误会非零退出并记录当前模块 failed。已有结果目录拒绝重用。模块执行完成与其中可选检验 completed/skipped 的状态分别保留，不把 GO 当作检验已完成。
