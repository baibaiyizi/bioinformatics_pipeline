# 空间转录组模板

复制整个目录到新项目，修改 `config/spatial_config.yaml`、样本 manifest 和 marker 表。输入支持 Visium 目录、AnnData 或表达矩阵加坐标。各平台均须确认输入是原始整数 counts；现有通用读取器不代表所有厂商格式均已实测。

```bash
bash run_spatial.sh --help
bash run_spatial.sh --list
bash run_spatial.sh --check    # 只检查，不创建输出
bash run_spatial.sh --steps 1,2,3,4,5,7,8,9,11,13,14
```

默认不做批次整合。`integration.method: harmony` 是显式选择；缺少依赖或批次字段时失败，不自动替换方法。模块依赖的对象必须先由上游生成。配置、源码和输入文件内容共同绑定运行签名；发生改变或旧结果无签名时拒绝复用。下游完成后重跑上游必须使用新的输出树。帮助、列表和预检不会创建分析结果。

## 实际能力

| 步骤 | 当前行为 |
|---|---|
| 1–5 | 输入检查、对象构建、QC、归一化、聚类 |
| 6–8 | 样本内空间统计、空间域、marker 注释 |
| 9 | `marker_score` 导出 marker 分数，**不是细胞丰度或去卷积**；`cell2location` 是单独的可选实现，需要参考和专门环境，尚未在本次运行 |
| 10 | 基于坐标与标签的 segment 描述 |
| 11 | 从原始 counts 归一化产生 NMF 程序，并输出每样本的表达和程序描述汇总 |
| 12–14 | 空间邻接/可选配体受体、样本汇总、结果索引 |

步骤 11 已删除 spot 级组间 Welch 检验，不产生条件 P 值，也不把 scale 后的均值差称为 logFC。输出 `inference_status.json` 明确记录 `cross_sample_inference=not_implemented`。样本汇总中的 section/sample 不自动等于独立受试者，不能据 spot、segment 或切片数扩增生物学 n。

`docs/`、`config/method_options.yaml` 和平台 profile 是方法参考资料；它们不是可执行开关注册表。RCTD、SPOTlight、BANKSY 等未接入此 runner，不能通过修改这些参考文件启动。实际参数入口仅为 `spatial_config.yaml`。

## 配置与输出

- 样本清单：`01rawdata/spatial_info.tsv`；`sample` 唯一，组别显式给出。扩展表可记录 subject、section 和平台身份。
- 标记：`config/marker_genes.tsv`，须按物种和组织审阅。
- 输出：`result/<步骤>/tables/`、`plots/`，AnnData 保留在 processed/integrated 目录。
- 去卷积默认分数文件名为 `marker_scores_by_spot.tsv`，另有方法身份记录。
- 步骤 11 目录为 `11_programs_descriptive`；旧名称与旧 spot-DE 文件不再被消费。

## 依赖和验证

沿用 `environment.yml` 的环境定义；不要把服务器上任意 Python 视为等价环境。默认 base Python 上 Scanpy/AnnData 曾发现不匹配；本次针对性测试直接使用 AnnData 读取已有对象，详见 `VALIDATION.md`。合成 H5AD 的 11 步核心闭环已实际运行；Visium 原始输入、Squidpy/Harmony/Cell2location 仍为 `not_run`。

```bash
python -m unittest discover -s tests -v
```

来源与本次提取/修订记录见 `provenance.tsv`。这次验证覆盖合成 H5AD 核心闭环、样本描述、counts 来源和禁止伪重复检验，不构成真实空间数据全流程复现。
