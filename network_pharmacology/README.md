# 网络药理学证据图最小闭环

本流程只接受已经版本化的本地节点表和证据边表，完成 ID/物种/来源审计、证据加权图、网络指标、删边稳健性和候选证据卡。Open Targets、STRING、PrimeKG 等远程数据应先由 `sci_evidence` 取得并保存查询日期、版本和原始响应；本入口不会在未经授权时上传本地数据或隐式联网。

```bash
python scripts/run_network.py --config config/example.yml
python -m unittest discover -s tests -v
```

节点中心性、网络邻近和候选总分只用于**后续验证的优先排序**，不等于治疗有效、靶点因果或分子直接结合。任何机制主张必须回到来源证据与独立实验。
