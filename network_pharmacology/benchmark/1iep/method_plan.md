# 1IEP 公开证据网络方法计划

## 知识复用记录

| knowledge_id | candidate_idea | decision | reason | adaptation | validation |
| --- | --- | --- | --- | --- | --- |
| github-mims-harvard-primekg | 异质证据图的实体—关系分层 | background_only | 历史架构可帮助理解 MultiDiGraph，但仓库已声明被取代 | 不下载数据；每条边回到当前官方 API 并保留来源版本 | PrimeKG 不出现在运行依赖、请求或评分中 |
| github-mims-harvard-optimuskg | 大型统一知识图谱与 MultiDiGraph | reject | 当前任务是最小、种子驱动、逐来源可审计检索，完整图谱增加下载与证据折叠 | 仅采用通用 MultiDiGraph 数据结构，不采用其数据或 client | OptimusKG 不安装、不下载；核心 API 独立留痕 |

## 实际反馈

| knowledge_id | outcome | 实际验证 | 可复用边界 |
| --- | --- | --- | --- |
| github-mims-harvard-primekg | not_run | 运行环境与 137 次公开 API 请求均不依赖 PrimeKG | 继续只作历史背景；不能把旧整合边当作当前直接证据 |
| github-mims-harvard-optimuskg | not_run | 未安装、未下载，仍完成 30 target、441 evidence edge 的闭环 | 官方 API 已满足当前问题时继续拒绝完整图谱 |
| synthesis-network-pharmacology-core-api-loop | supported | 记录 Open Targets、ChEMBL 37、PubChem、STRING 12.0、Ensembl 116、UniProt 2026_02、RCSB 和 AlphaFold 的参数、版本、分页及响应校验和；26 个候选满足至少两条 lane | 网络排名仅作候选优先级；本次无显式 Reactome 背景，故不运行 ORA |
