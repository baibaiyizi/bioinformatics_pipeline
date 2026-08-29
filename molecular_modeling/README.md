# AutoDock Vina—GROMACS 分阶段闭环

唯一入口按固定阶段执行：`preflight → dock → md-prepare → md-smoke → md-production → analyze`。对接使用 AutoDock Vina，Meeko 负责 PDBQT 准备和结果还原，PoseBusters 与 RDKit 重原子对称校正 RMSD 负责 redocking 硬门；禁止用 Open Babel 猜测 PDBQT 键级。

```bash
python scripts/run_modeling.py preflight --config config/example.yml
python scripts/run_modeling.py dock --config config/example.yml
python scripts/run_modeling.py md-prepare --config config/example.yml
python scripts/run_modeling.py md-smoke --config config/example.yml
python scripts/run_modeling.py analyze --config config/example.yml
python -m unittest discover -s tests -v
```

公开短基准及输入溯源位于 `benchmark/1iep/`。独立环境、GROMACS 2026.3 源码构建和分级验证分别由工作区 `server_migration/install_sci_molecular_modeling.sh`、`build_gromacs_2026_3.sh`、`verify_sci_modeling_stack.sh` 执行；不安装或调用 OpenMM、DiffDock。

MD 参数化固定为 AmberTools `ff14SB + GAFF2/AM1-BCC + TIP3P`，完整体系经 ParmEd 转为 GROMACS。受体由 Reduce 建立氢键网络，HIS 再依据实际 HD1/HE2 原子显式记录为 HID/HIE/HIP；无法判定就阻断。金属、LINK/共价体系、非标准残基、质子化不明确、缺 native/阳性/decoy、少于三个 docking seed 或少于三个 MD replica 都会失败可见。生产运行没有默认时长，必须显式给出时长和批准 pose，并通过真实 CUDA 设备、`nvidia-smi` 与 CUDA 版 GROMACS 三重门控。

短 smoke 会运行三条独立 seed、checkpoint 续跑并检查 NaN/爆炸，再用 MDAnalysis 输出 RMSD、RMSF、Rg、口袋接触、氢键、关键距离和 replica 一致性。MDAnalysis 2.10 尚不支持 GROMACS 2026.3 的 TPX v138，因此分析固定读取生成 GROMACS topology 的同源 Amber `system.prmtop` 和 XTC；二者保持相同原子顺序，且 prmtop 保留氢键分析所需的键连接。它只证明环境与解析闭环可执行；`docking score` 不等于亲和力，单条或短轨迹不等于稳定性或机制。
