# 依赖

沿用现有环境，不自动安装或升级依赖。

描述流程需要 Rscript、ggplot2、dplyr、tidyr、stringr、purrr、tibble、knitr、data.table；endo 当前 FASTA 复核需要 Biostrings；差异与富集沿用模块原有 edgeR、multiMiR、OrgDb 等条件依赖。具体版本以实际运行的 sessionInfo 和验证记录为准。

SPORTS 上游需要 Perl、Bowtie1、bowtie-build、gzip、xargs、Python3、sha256sum，剪切需要 cutadapt；可选 FastQC/MultiQC 依配置检查。物种数据库、基因组 FASTA 与 Bowtie1 索引由配置定位。R 描述运行复用已有 SPORTS 表，不代表 FASTQ、剪切或 Bowtie1 已重新执行。

入口和只读预检仅使用 Bash 与 Python3 标准库。
