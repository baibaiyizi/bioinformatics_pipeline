# 验证记录

验证日期：2026-09-29。仅修改模板目录；来源项目只读。使用现有 Python 3.13.2、Scanpy 1.12.4、AnnData 0.12.19、PyDESeq2 0.5.4、pytest 9.1.1；未安装或升级服务器软件。完整依赖快照见 requirements-validated.txt。

## 已执行

- `python -m pytest tests`：**18 passed，107.26 秒**。覆盖独立单位/混池计数、同样本技术库聚合原始整数 counts、未知身份不合并、单位/条件/批次冲突、缺 batch、批次混杂、原始计数整数门、植入效应与反向对比及行序、重复 library ID、TSV 空字段、配对设计拒绝、可选 I1 和完整 FASTQ 配对、MD5/稳定性/不改原文件、复制模板后 n=1 端到端、帮助/列表/validate 不写结果、推断拟合失败保留 GO + failed。
- 将整目录复制到临时位置，使用独立生成器创建 **3+3 个生物学样本、每样本 2 个文库、共 720 个细胞**；从原始矩阵运行 QC→图谱→marker 注释→描述→**实际 PyDESeq2**→报告→审计全部通过。type_A、type_B 各由 6 个生物学样本构成，每样本 60 个细胞和 2 个技术库。
- 合成 Gene000 在 treated 中植入上调：type_A log2FC=4.177505、P=8.680379e-241、padj=1.302057e-238；type_B log2FC=4.304589、P=4.172700e-273、padj=6.259050e-271。两者实际状态均为 eligibility=GO / execution=completed。这些仅证明计数、方向和执行合同，不是生物学发现。
- 复制后的 n=1 端到端自动降为 skipped_ineligible，保留描述及仅表报告，不生成 DE 表。另一次独立 n=1 合成运行启用 figures=true，实际生成 4 张 PNG 和 4 张 PDF（3 张 UMAP + 捕获组成）及有图报告，审计通过。
- `bash -n run_scrna.sh` 与 Python 源码/测试/示例源码编译检查通过。

测试中 PyDESeq2 对合成数据发出两条 parametric dispersion trend 不收敛提示，按库内实现使用 mean-based trend；实际拟合完成、P 值和方向断言通过。该记录不代表真实数据拟合诊断。

## 能力与边界

| 能力 | 本轮证据 |
|---|---|
| 原始 h5ad → 核心分析 → 报告 | 合成端到端实际运行 |
| 3+3 实际 pseudobulk DE、技术库聚合 | 已运行，独立样本数未扩增 |
| n=1 / 元数据不足 / 混杂 | 测试通过，只描述 |
| FASTQ MD5、稳定性、完整成对与可选 I1 | 小型 gzip FASTQ 实际运行；无写句柄分支源码审查 |
| Cell Ranger 原始数据上游 | `not_run`；可执行入口已接线，真实参考/FASTQ 未运行 |
| Scrublet 启用分支 | `not_run`；默认关闭，失败会终止 |
| 参考映射 / SoupX / 组成推断 | 通用核心未实现或未执行，报告明确说明 |
| 胰岛 velocity / 通讯 / GRN / beta 扩展 | `not_run`；仅隔离参考源码，不宣称可直接运行 |
| 真实数据复现 / 独立重复 | `not_run` |

执行日志临时目录：`/tmp/scrna-template-validation-oZSTuL`。小型可保存验证摘要、被测源码哈希和 pytest 日志在 `validation/`。重新验证使用 README 中 fixture 命令，禁止读取真实项目配置制造合成结果。
