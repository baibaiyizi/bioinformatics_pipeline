# 联合分析模板

一个物种内的跨组织、跨阶段和跨组学分析使用一个联合入口。复制后编辑 `config/project.yml` 的真实来源、方向与质量状态，再运行：

```bash
bash run.sh --check
bash run.sh --all
bash run.sh --stage liver
```

默认配置含 `INPUT_REQUIRED`，不会搜索真实课题作为替代。相对路径从联合分析根解析；物种、基因组、方向、样本身份和生产状态均需显式核验。

## 标准输入与职责

- RNA：完整合格cohort的DE、表达矩阵、样本表与生产签名；`sampleA=comparison`、`sampleB=reference`。单组学producer为RNA模板的 `scripts/run_full_cohort_de.R`。

  该通用producer必须显式传入`--reference CONTROL --comparison TREATMENT`（替换为样本表实际标签），不会采用真实课题的组名。
- WGBS：含明确组别和基因注释的DMR，坐标为1-based闭区间。DMR区域阈值不等于区域FDR。
- ATAC：差异区域、基因注释、效应和明确方向。最近基因注释不等于远端靶向验证。
- 代谢组：保留原始P、FDR、特征身份、鉴定等级和歧义；缺少反应/基因映射时报告连接缺失。
- 单细胞：描述性细胞类型背景保持自身实验单位，不能以细胞数代替动物数。

`--check` 只读；`--all` 依次执行 `rna`、`wgbs`、`liver`、`programs` 和报告。`--stage` 可独立更新；沿用其他阶段前须匹配输入、源码和产物签名。失配的阶段不出现在当前结论中。重算前的联合产物保存到本地history。

报告直接更新本模块README，HTML位于result/index.html；所有图表保存在所属阶段的result，不复制交付版。报告和图表也核验签名。执行前记录STARTED并失效旧成功报告；渲染失败记录FAILED。科学阶段未齐时总状态PARTIAL，不把局部完成写成全流程完成。区域feature注释与最近基因链接分别保留。

GO程序表是本地注释交汇，不提供富集P值或通路活性。各层效应与冲突分别保留，不输出未经校准的总分。历史上游来源不会因本次连接计算自动变为已复现。

## 环境与复制

沿用现有 Python pandas/numpy/matplotlib/PyYAML/Markdown，R GenomicRanges、AnnotationDbi、GO.db及所选OrgDb；不自动安装或升级。需要smallRNA方法比较时，在中心模板中单独选择 `integration/method_comparison`；复制通用integration不附带这个子模块。项目运行自己的代码副本，来源以 `template_snapshot.tsv` 冻结。
