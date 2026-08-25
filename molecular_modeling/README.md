# 分子对接与分子动力学 preflight

当前工作区没有 Vina/Open Babel/PoseBusters/OpenMM/MDAnalysis 的完整环境，因此本目录首先提供一个可运行、可测试的**结构与执行契约审计闭环**。它不会生成伪造的 docking 或 MD 结果。

```bash
python scripts/preflight.py --config config/example.yml
python -m unittest discover -s tests -v
```

preflight 检查：受体 PDB 基本完整性、配体文件/质子化/电荷方法、对接盒、redocking 与阳性/诱饵对照、多随机种子、MD 独立重复、力场/水模型/配体参数化声明，以及本地执行器和分析依赖。输出 `preflight.json`、逐阶段状态表和执行计划。

只有状态为 `ready` 的阶段才能进入真实运行。`docking score` 不等于结合亲和力，单条或单次短 MD 不等于复合物稳定或机制证明；这些边界写入每次审计产物。
