# 模板更新交付记录

本轮仅实施 `model` 独立仓库的生信模板计划。真实项目作为只读来源；未全量重算、未安装或升级服务器环境。已有未提交的 `tests/test_sci_agent_config.py` 修改按基线哈希保留。ML、网络药理、分子建模和 AI harness 实现不在范围内。

## 使用与完成边界

七类主要组学模板均以独立目录交付。入口、配置和支持范围见[总 README](../README.md)，推断单位及未实现能力见[能力矩阵](../analysis_capability_matrix.md)。外部工具和参考需按项目明确配置；合成输入是验收材料，不是生产默认数据。

| 模板 | 已执行的最小闭环 | 详细记录 |
|---|---|---|
| bulk RNA | 500 基因×4 独立样本 counts→DESeq2→01–05→HTML；22项身份/迁移检查及DE/rMATS方向小表 | [VALIDATION](../rnaseq/VALIDATION.md)、[重放与证据](../rnaseq/validation/README.md) |
| WGBS | 合成 CpG→methylKit/DSS 下游→报告；精确模块和未运行分支见模板记录 | [VALIDATION](../wgbs/VALIDATION.md) |
| ATAC | 合成 consensus 整数 counts→DESeq2→图表/报告；BAM原生DiffBind为独立输入路径 | [VALIDATION](../atac/VALIDATION.md) |
| scRNA | 3+3 独立样本、每样本2技术库、720细胞→QC/图谱/描述→实际pseudobulk→报告；18项测试 | [VALIDATION](../scrnaseq/VALIDATION.md)、[证据](../scrnaseq/validation/) |
| PANDORA | 3+3 合成SPORTS表→表达/组成描述→图表/报告；整数counts独立edgeR方向验证 | [VALIDATION](../smallrnaseq/VALIDATION.md) |
| LC-MS | 60特征×6样本峰表→QC/limma/图表/报告；有/无填峰前矩阵、反向比较、无效单位与极性 | [VALIDATION](../metabolomics/VALIDATION.md)、[摘要与表](validation/metabolomics/) |
| 空间转录组 | 2样本×40spots×150genes H5AD→11步描述流程→5对图和报告；8项合同检查 | [VALIDATION](../spatial_transcriptomics/VALIDATION.md)、[测试日志](validation/model_spatial_tests.log) |
| multiome / 空间代谢组 | 已有各2项测试；仅counts/长表最小入口 | [multiome日志](validation/model_multiome_tests.log)、[空间代谢组日志](validation/model_spatial_metabo_tests.log) |

所有 FASTQ/真实 mzML/真实 Visium/imzML 的完整生产上游仍为 `not_run`。没有把过去真实项目的结果冒充本次模板复现。表观组的注释和富集是否实际执行，以各模板记录为准；空间正式条件组间推断为 `not_implemented`。

## 可重放检查

从仓库根运行：

```bash
python template_audit/verify_templates.py
# 以下需要对应既有 R/Python 依赖；各自在 /tmp 新建副本，不读真实项目结果
python template_audit/replay_metabolomics.py
python template_audit/replay_spatial.py
python -m unittest discover -s spatial_transcriptomics/tests -v
python -m unittest discover -s singlecell_multiome/tests -v
python -m unittest discover -s spatial_metabolomics/tests -v
```

`replay_spatial.py` 应使用具备已记录 Scanpy/AnnData 版本的解释器。脚本不自动安装依赖。其他模板的 fixture 与重放命令在各自 README/validation 中。

`verify_templates.py` 只做 Shell/R/Python 语法与来源完整性检查，不启动分析。结果写入 [syntax_and_integrity.json](syntax_and_integrity.json)：区分解析通过、来源文件保持不变、用户预存修改保留以及 Git 空白检查。Rmd 仅提取R代码块进行parse，不构造新的R测试框架。

## 来源和证据保留

- [baseline.json](baseline.json)：编辑前仓库状态、来源SHA256和用户预存修改SHA256。真实项目源码基线在最终检查中复核。
- [provenance_index.tsv](provenance_index.tsv)：七类模板来源清单的位置、清单SHA256、提取边界及实际验证记录。各模板 `provenance.tsv` 逐文件记录来源与哈希；该来源哈希不是提取后目标源码哈希。
- [validation/](validation/) 及各模板 validation 目录保留小表、状态、日志或输出哈希；不复制真实数据和大型图。临时执行目录仅作本轮追溯，不是模板运行依赖。
- 检查未执行、必要输入缺失、输入/配置改变或输出无签名，分别保留 `not_run`、`MISSING`、`STALE`、`UNVERIFIED` 等状态。GO资格与实际completed执行状态分开。

本轮修复包括空字段错列、08/09编号解析、样本身份和方向推断、技术库伪重复、PCA挑样本、真实项目路径回退、缺填峰前矩阵却声称检出率、spot条件检验，以及失败上游残留被误用的问题。具体保证范围以实际被测入口为准。

最终统一检查：144 个模板 Shell/R/Python 源文件语法通过；265 个来源基线文件未改动；来源清单另核对 92 个外部源文件哈希；用户预存修改保留，Git 空白检查通过。scRNA 18 项、空间转录组 8 项、multiome 与空间代谢组各 2 项 Python 测试分别通过（共 30 项）；R 小数据与合同检查单列，不混入此计数。
