# 空间代谢组最小闭环

流程边界：`imzML/ibd 或长表 → 输入审计 → 像素/光谱 QC → m/z 分箱与 feature matrix → TIC 归一化 → ROI 汇总 → 样本内空间统计 → 注释证据表`。

```bash
python scripts/run_spatial_metabolomics.py --config config/example.yml
python -m unittest discover -s tests -v
```

`mode: imzml` 需要 `pyimzML`；缺少依赖时流程明确返回 dependency error。合成 smoke test 使用 `mode: long_table`，并不宣称替代原始 imzML 验证。Cardinal/SpatialData 属于后续增强层；当前入口保留可被二者读取的坐标、feature 和 ROI 表。

统计边界：Moran's I 只描述单个样本内的空间结构；像素和 ROI 不是生物学重复。跨组差异必须先汇总到样本层并交 `sci_analyst` 建模。未经标准品确认的注释不能标为 level 1。
