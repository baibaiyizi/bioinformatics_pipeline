# 本次验证

2026-09-29，使用当前 R 环境，在三个独立的含空格临时项目中运行；不读取真实项目结果。

- 峰面积输入：60 个特征 × 6 个独立样本（control/treated 各 3），5 个植入方向效应。
- 无填峰前矩阵：只读检查无输出；完整峰表→预处理→limma→PDF/PNG→中文报告通过；检出率全为 NA，状态为 `not_available_no_prefill_matrix`。
- 有填峰前矩阵：指定全缺失特征被检出率过滤；填峰前后矩阵独立保存。
- 交换 case/control：log2FC 精确翻转（绝对容差 1e-9），P 值一致（相对容差 1e-8）。
- mzML 每个极性缺少部分样本（control 仅 POS、treated 仅 NEG）：预检明确拒绝，不创建结果，避免采集极性与组别混杂。
- 重复 biological_sample_id：预检查失败，不创建结果目录。
- 模块状态按 inventory/xcms/preprocess/plots/annotate/stats/pathway/report 分别记录；执行时校验输入生产步骤，已实际验证拒绝失败上游遗留。
- 现有源码 R parse 与 Shell `bash -n` 通过；依赖实测版本见 `dependencies.tsv`。

原始 mzML 的 XCMS 峰提取、WIFF 外部转换、真实注释标准品、外部通路定量均 `not_run`。上述结果为合成峰表闭环，不是独立真实数据复现。验证目录 `/tmp/model-lcms-validation-ww308d40` 属于临时执行证据，不能作为模板运行依赖。

可重放脚本与本轮输出小表/哈希见 `../template_audit/replay_metabolomics.py` 和 `../template_audit/validation/metabolomics/`。
