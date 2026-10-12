# 模板检查与历史验收

当前模板目录为 `rnaseq`、`wgbs`、`atac`、`scrnaseq`、`smallrnaseq`、`metabolomics`、`spatial_transcriptomics` 等，使用方法见[模板总入口](../README.md)。

```bash
# 从 templates 目录运行，核验当前 Shell / R / Rmd / Python 语法
python3 template_audit/verify_templates.py
```

结果写入 [current_syntax.json](current_syntax.json)，仅证明所列源码语法通过。真实数据重跑、合成输入闭环、输入资格和生产完成状态分别记录，不能互相替代。

2026-09-29 的来源基线、预存修改核验、旧源码检查与验收记录保留于 [history/20260929](history/20260929/)。其中绝对路径表示当时位置，不是当前入口或运行依赖。各模块 `provenance.tsv` 同样记录提取时的来源，不能当作当前复制版本签名。

smallRNA 现采用 2026-10-11 新版 PANDORA 的本地冻结实现；当前验证见其 [VALIDATION.md](../smallrnaseq/docs/VALIDATION.md)，不能沿用旧编号模块的验收结论。

其他模板的具体重放方法和生产上游 `not_run` 边界，仍以各模块 README、VALIDATION 和 [能力矩阵](../analysis_capability_matrix.md) 为准。原有小数据证据保留于 [validation](validation/)。
