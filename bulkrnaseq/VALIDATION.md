# 验证记录

验证日期：2026-09-29。实际修改只在本模板；真实项目源文件只读，未全量重算、未安装或升级环境。记录区分源码提取、合成数据运行和真实数据复现。

## 已执行

- `bash -n`：4 个入口和配置 helper 通过。R 解析：18 模块、内置 helper、上游嵌入 R 共 33 个输入通过；Python helper 编译检查通过。
- 含空格临时复制目录 `/tmp/bulk template smoke 5tnnzid_/project`：500 个真实 mouse Ensembl ID 的合成整数 counts、4 个独立动物、2+2 组，样本表故意交错排序，运行 `rnaseq_second.sh --all` 返回 0；得到 500 行 DE 表，全部 `sampleA=treatment, sampleB=control`，保留 pvalue 与 padj。
- 同一隔离目录运行 `run_rmd.sh --modules 2,3,4,5`（自动含 01）返回 0：01–05 均 completed；06–18 skipped；输出 17 PNG、17 PDF、9 CSV 和 module_status.tsv；报告 inventory 为 17 图、10 表、6 个目录/状态分区。没有把该合成数据结果当成生物学发现。
- HTML 报告 SHA256：`782fdf092af47e8faf76517400fbc8bae2990edf784eda638d11b333d63b6806`。仅有 module_status.tsv、没有 PNG 的独立目录也成功生成报告。
- 22 项输入/入口检查通过：四入口 help/list/preflight 全通过且未增加文件；重复 sample、重复 biological_id、未知 biological_id、第三组明确拒绝；样本表行序重排不改变 groups 或 samples.txt；更改 counts 后拒绝 STALE；FASTQ 原文件名与 sample 不同仍输出 sample_1/2.fastq.gz；第一阶段 dry-run 未增加文件。
- rMATS 专项：b1=control 和 b1=treatment 均得到统一 treatment-control 差值；样本表重排不改变 metadata 指定的 PSI 次序；原始 PSI/IncLevelDifference/PValue 保留；`0.1,NA,0.3` 不左移；缺 metadata 和错误 sample_key 明确失败。该项是小表方向验证，未运行 rMATS 二进制。
- 模块 05 的实际方向代码以反向 DE 输入验证：交换 baseMeanA/B、翻转 log2FC/stat、输出方向审计正确。
- 首轮暴露并已修复：旧模块调用新 setup 已移除的主题函数；无 PNG 报告禁止生成；根级状态表缺 module title 时异常。修复后对应场景通过。

## 环境

实测 R 4.5.1、Ubuntu 22.04.5，调用现有 R 环境；具体包版本见 `dependencies.tsv`。DESeq2 在这组合成数据上报告 parametric dispersion trend 不适用并采用其自带 local fit；这是实测运行信息，不是生产数据评价。

脚本不自动安装依赖。联网富集依赖及上游工具是否可用于真实项目需在实际复制目录检查。

## 未执行与限制

- `not_run`：真实 FASTQ fastp/HISAT2/featureCounts 全量上游；rMATS/APA/DER 生物数据；06–18 富集、GSVA/TF、PPI、WGCNA、免疫与整合模块。这里仅做来源提取、接口检查和语法检查，不能宣称真实数据复现。
- `MISSING`：模板示例没有携带参考基因组、FASTQ 或生产 counts。直接对未配置的示例运行 preflight 应明确失败。
- `STALE/UNVERIFIED`：签名变化、既有输出无签名均拒绝；更改配置不能以旧文件存在证明执行。FASTQ 签名使用大小/mtime；本轮未执行 FASTQ 全量 checksum 审计。
- 首版仅实现小鼠独立两组设计。配对、重复测量、复杂协变量和技术文库聚合上游不在该模板内实现。

小数据日志、22 项边界清单、模块状态、方向审计、分析参数、R 环境与汇总已保存在 `validation/evidence/`。可重放入口见 `validation/README.md`，从脚本位置定位源码；不依赖原 `/tmp` 目录。未复制大图或真实项目数据。
