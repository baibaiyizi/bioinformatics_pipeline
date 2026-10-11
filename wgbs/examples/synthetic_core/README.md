# 合成 methylKit → DSS 示例

四份完全合成的 300 CpG methylKit 文本，每个位点深度 100；control 的甲基化约 20%，treated 约 80%。GTF 的坐标和 transcript 是人为构造的解析夹具，借用一个可映射的鼠基因 ID 仅用于验证注释分支；不能作为真实参考注释或生物学结论。

把完整模板复制至新的临时目录，执行 `bash examples/synthetic_core/run_example.sh`。示例通过显式 cached 模式运行正式上游 12/13（methylKit 合并）与下游 1/3（DSS、DMR、注释），生成图、表和 HTML。仅使用既有 R 环境和 Pandoc，不下载安装。

示例准备过程记录合成 staged 输入来源，不执行、不证明 FASTQ QC、比对或甲基化提取。样本表中的 FASTQ 名称只是未运行上游的合同占位，cached 模式不读取它们。`metadata/synthetic_input_origin.tsv` 与模块状态共同界定实际执行层级。所有四个独立合成样本保留，CpG 和 DMR 不替代生物学 n。
