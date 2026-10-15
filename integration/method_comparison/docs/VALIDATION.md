# 方法比较模板核验

下列为2026-10-11/12的历史执行记录；集中维护凭据已按研究者要求清理，不代表现役源码重新验收。

2026-10-11 从实际血清分析抽离，项目使用独立本地源码。目录精简后唯一方法比较模块直接位于项目 integration 根，describe、differential 按现有门禁真实重放。

历史 17 张捕获 CSV 的逐值身份、180 个类别×文库计数、单供体三技术配对和观察深度边界仍是本轮核验对象。

report 直接引用原始图表，生成模块 README 与 result/synthesis 的 HTML、SHA 凭据。新项目必须重新通过输入门并实际运行。未重跑 FASTQ 或 SPORTS；Poisson 曲线超出观察深度和分数注释稀释仍为 not_run。

2026-10-12 完成连续编号：01_setup、02_method_similarity、03_spikein_qc、04_differential_capture、05_capture_bias_diagnostics。外参模块的五份本地 Rmd 按 profiles/01–05 调度；没有运行时复制单组学代码。

最终源码下，血清单组学、血浆单组学 describe 与血清联合 --all 均实际退出 0。联合分析 513 张 CSV/TSV（含图源数据）与移动前逐文件 SHA256 相同；当前 153 对描述图和 16 对差异图 PDF/PNG 均存在、MD5 与清单一致、源数据齐全，旧图类型遗漏为 0。describe、differential 各 60 项当前输入签名有效，报告 SHA 凭据有效。工作区与 /tmp 入口检查、从 /tmp 实际生成报告，以及项目/模板各 7 项输入门测试全部通过。

图形抽样涵盖 8 种主要图型；火山图沿用原始布局，极端值压缩基线区域，不能据此声称全部图已达到投稿排版标准。DonorN=1、技术配对、外参独立 QC 与 not_run 边界均保留。
