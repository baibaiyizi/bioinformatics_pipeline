# 1IEP–imatinib–ABL1 短闭环基准

本目录只用于验证公开结构检索、Vina 三 seed redocking、PoseBusters、AmberTools/ParmEd 参数化、GROMACS 三重复短 smoke、checkpoint 续跑和轨迹解析是否能在当前环境执行。它不是药效排序、生产 MD、自由能或机制研究。

## 公开输入

- `1iep_holo.pdb`：RCSB PDB 的公开 `1IEP` 原始 PDB；仅使用链 A。
- `native_imatinib.sdf`：AutoDock Vina `v1.2.7` 官方 basic docking 教程的共晶 STI/imatinib 参考结构。
- `positive_dasatinib.sdf`：PubChem CID `3062316` 的公开 3D SDF，只作已知 ABL 抑制剂的流程阳性对照。
- `decoy_caffeine.sdf`：PubChem CID `2519` 的公开 3D SDF，只作无关化学结构的执行性 decoy；它不是性质匹配的 DUD-E decoy。

文件校验和见 `checksums.sha256`。共晶口袋盒参数与 Vina 官方教程一致：中心 `(15.190, 53.903, 16.917)` Å，尺寸 `20 × 20 × 20` Å。

## 结构边界

原始 1IEP 在链 A 缺少 N 端残基 223–224 和 C 端残基 499–515；这些残基不构成教程共晶口袋，本基准只接受其作为能力测试。该决定不能自动迁移到真实项目。金属、共价配体、非标准蛋白残基或无法明确质子化的体系仍会阻断自动参数化。

## 验收边界

- native ligand 至少两个 seed 的最佳 PoseBusters 有效 pose 达到重原子对称校正 RMSD `≤2 Å`，才生成 MD handoff。
- 阳性对照和 decoy 只检查流程与几何，不以这三个小分子建立分类性能主张。
- 三条短轨迹只验证运行、checkpoint 和解析，不证明稳定性、结合强度或机制。
