# 依赖

使用服务器既有 R/工具环境，未升级或安装。R 的实际包版本记录在 validation/R_packages.tsv，每次运行保存 audit/sessionInfo.txt。

最小下游：R、ggplot2、dplyr、tidyr、stringr、purrr、tibble、knitr；Python 3 标准库用于样本表预检查和合成输入，无额外 Python 包。data.table 存在时用于大表读写。

可选专题模块根据源脚本使用 scales、patchwork、ggrepel、pheatmap、ggalluvial、readr 等；缺依赖有显式说明，不等于该图已生成。实际 edgeR 还需 edgeR、limma、statmod；miRNA 靶标/富集需 multiMiR、clusterProfiler、配置的 OrgDb，显式 MSigDB 来源还需 msigdbr。缺包不会自动安装。

原始上游：SPORTS1.1 的 sports.pl、Perl、Bowtie1/bowtie-build、gzip、xargs；预剪切需要 cutadapt；FastQC/MultiQC 为显式可选。复用 workspace/server_migration/envs/smallrna.yml 的环境层思路，但模板执行不读取该工作区路径；在目标服务器提供 PATH 和配置指定的工具/参考即可。原始上游本轮未执行。
