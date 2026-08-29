# 可审计的公开证据网络药理闭环

入口以疾病、化合物和 target seed 为唯一种子，定向查询 Open Targets、ChEMBL、PubChem、STRING、Reactome、Ensembl、UniProt、RCSB 与 AlphaFold DB。Ensembl 使用官方 BioMart 的精确 ID 查询并记录 release；ChEMBL 活性采用“候选 Ensembl → reviewed UniProt → 同物种单蛋白 target → 化合物 parent ID”的定向路径，不拉取化合物全部无关活动。每次请求保存参数、分页、访问时间、API/数据库版本（来源未报告时明确写 `not_reported`）、响应 SHA-256 和规范化结果；分页不全、ID 一对多及 API 失败都会终止，不能用看似合理的残缺结果继续。

```bash
python scripts/run_network.py --config config/example.yml
python -m unittest discover -s tests -v
```

独立环境由 `/home/h1028/workspace/server_migration/install_sci_network_pharmacology.sh` 建立；脚本只使用显式 conda-forge URL，不修改 base 或全局频道。PrimeKG 只保留为被取代的历史背景，OptimusKG 不安装、不下载，也不进入评分依赖。

核心输出为 `standard_entities.tsv`、`evidence_edges.tsv`、`evidence_multidigraph.graphml`、`candidate_priority.tsv`、`leave_one_lane_out.tsv`、`string_threshold_sensitivity.tsv` 与 `modeling_handoff.tsv`。疾病关联、physical PPI、显式背景 Reactome ORA、结构可用性和直接 compound-target 证据分别计算百分位，再等权汇总；少于两个独立 lane 的候选不会进入正式排序。

人源与小鼠不混合：小鼠项目中的人源疾病证据必须经 Ensembl ortholog 映射并标为 `human_ortholog`。STRING physical/functional 分开存储；没有显式物种背景就不运行 Reactome ORA。网络排名只表示候选优先级，不等于疗效、因果、结合亲和力、复合物稳定性或机制。
