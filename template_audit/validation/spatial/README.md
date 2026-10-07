# 空间转录组结果索引

依据当前输入签名及执行状态生成；文件存在本身不表示步骤完成。

|模块|执行状态|表格|图文件|
|---|---|---:|---:|
|01_validate_inputs|completed|2|0|
|02_build_objects|completed|1|0|
|03_qc_filter|completed|2|2|
|04_normalize_features|completed|1|0|
|05_integrate_cluster|completed|1|6|
|06_spatial_neighbors_stats|not_run|0|0|
|07_spatial_domains|completed|1|2|
|08_cell_type_annotation|completed|2|0|
|09_deconvolution|completed|2|0|
|10_segment_objects|not_run|0|0|
|11_programs_descriptive|completed|6|0|
|12_neighborhood_communication|not_run|0|0|
|13_multi_sample_summary|completed|4|0|
|14_report_index|completed|0|0|

## 解释边界

- 本流程没有实现正式的条件间推断；第 11 步只输出样本描述，不产生条件 P 值或 logFC。
- spot、segment、切片均不能替代独立生物学样本。
- marker score 是表达评分，不能解释为细胞丰度。
- 空间域和 segment 为计算标签，需要组织学及 marker 审阅。
- Cell2location 是未在本轮验证的可选分支；RCTD、SPOTlight 等方法未接入本入口。
