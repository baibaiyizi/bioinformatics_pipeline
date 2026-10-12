# 验证记录

日期：2026-09-29。仅在模板或隔离副本操作，未重算真实空间项目，未升级环境。

## 实际运行

- Python `unittest`：**8 项通过**。覆盖原始 counts 来源、按样本描述且无条件 P/logFC、marker-score 身份、组别冲突、缺 counts/空路径/缺配置、改变输入产生 STALE、无签名旧结果拒绝、失败上游不能供给旧对象、自定义 manifest/result 路径连续分步运行。
- 复制到含空格临时目录，使用独立生成的 **2 样本×40 spots×150 genes** H5AD；X 故意放置非原始值，counts 层保留整数。
- 运行 01–05、07–09、11、13–14 共 11 步：只读输入检查→counts 构建→QC→归一化/PCA→聚类→空间域→marker 分数→样本描述→报告。输出 300 行样本基因描述、5 PNG、5 PDF；不产生条件组间检验。
- 修复仅表报告、失败状态、配置与输入签名绑定；未运行步骤标为 not_run，不通过目录是否存在推断完成。
- NMF 合成测试可能达到 500 次迭代上限；输出 `nmf_fit_status.json` 明确记录实际迭代与重构误差。程序分解仅作计算演示，不作为稳定生物学程序的证据。

脚本可用 `template_audit/replay_spatial.py` 在新的临时副本重放（从仓库根以具备下列依赖的 Python 运行）。8 项模板测试位于 `tests/`，日志和小型执行证据收录在 `../template_audit/validation/`。

## 环境与未验证范围

合成核心运行使用现有 Python 3.13、Scanpy 1.12.4、AnnData 0.12.19 环境；具体版本见 `dependencies.tsv`。通用 base 环境的旧 Scanpy 与当前 AnnData 不匹配，不能当作等价运行环境；定向输入/描述测试只需要 AnnData 及数值依赖。

Visium 原始目录、矩阵加坐标输入、Harmony、Squidpy 步骤 06、segment 步骤 10、邻接/配体受体步骤 12、Cell2location 均 `not_run`。RCTD、SPOTlight 等未接入入口；正式条件组间推断为 `not_implemented`。合成 H5AD 闭环不等于真实空间数据复现。

执行快照的源码/输入哈希保存在合成副本的 run_signature；其后补齐默认 method 字段声明，已由最终 8 项测试复核。最终源码快照另存于 `../template_audit/validation/final_root_source_sha256.tsv`，不将二者混为同一次执行哈希。
