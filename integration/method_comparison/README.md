# 血清 smallRNA 方法比较模板

这是供体 2 方法 1/2 三个匹配技术分装、方法 3–6 单库、供体 3 独立外参 QC 的已运行示例。复制到课题的 human/integration 根目录，修改 [config/project.sh](config/project.sh) 指定本地 smallRNA 模块。真实供体、配对及外参设计改变时，需要同步调整输入门，不得把技术分装当作人群生物学重复。

```bash
bash run.sh --check
bash run.sh --all
bash run.sh --stage describe
bash run.sh --stage differential
bash run.sh --stage report
```

入口与资源定位不依赖启动目录；代码、缓存和结果均在项目本地。先核验单组学完成状态及当前输入/源码签名，再重建序列特征，并逐值核对标准 RNA 类别计数。方法 2−方法 1，配对模型为 ~PairID + Method。

步骤按 `01_setup、02_method_similarity、03_spikein_qc、04_differential_capture、05_capture_bias_diagnostics` 排序；外参类别图使用本地 `spikein_profiles/01–05` 源码，结果写入 `03_spikein_qc/profiles/01–05`，不在运行时复制单组学源码。

结果保存在 result/describe 和 result/differential；report 更新本模块 README，并在 result/synthesis 保存完整 HTML 和报告凭据，直接链接正式图表，不复制精选图或关键表。报告本身不参与科学拟合签名。

深度检出沿用 Poisson 期望模型，仅在实际观察深度以内计算；越界保持 not_run。该曲线不是精确无放回抽样。供体 3 序列层 200 次无放回稀释独立用于外参 QC，注释层分数计数不取整稀释。

[核验说明](docs/VALIDATION.md)记录的是实际血清示例，不能视为未来复制项目已运行。源码身份见 [SOURCE_SNAPSHOT.tsv](SOURCE_SNAPSHOT.tsv)。
