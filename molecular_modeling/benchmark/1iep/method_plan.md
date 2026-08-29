# 1IEP 短闭环方法计划

## 知识复用记录

| knowledge_id | candidate_idea | decision | reason | adaptation | validation |
| --- | --- | --- | --- | --- | --- |
| github-ccsb-scripps-autodock-vina | Vina 固定盒三 seed redocking | adopt | 官方 1IEP 教程直接匹配本基准 | 增加三 seed、重原子对称 RMSD 和硬门 | 至少两个 seed `≤2 Å` 且 PoseBusters 有效 |
| github-forlilab-meeko | 受体/配体 PDBQT 准备和 pose 还原 | adopt | 与 Vina 同步维护且保留分子结构信息 | 用 Meeko 输出 SDF，禁止 Open Babel 猜 PDBQT 键级 | CLI 实跑、RDKit 可读和原子身份一致 |
| github-maabuu-posebusters | pose 物理合理性检查 | adopt | 为 redocking 提供独立于 docking score 的几何门 | native 使用 `redock`，阳性/decoy 使用 `dock`；逐 pose 判定 | 最佳有效单 pose 通过全部布尔检查并保存逐项结果 |
| gitlab-gromacs-gromacs | 三重复短 MD、checkpoint 与轨迹 | adapt | GROMACS 是目标生产引擎，但当前 GPU 设备门未通过 | 仅执行 CPU smoke；生产保持 `resource-limited` | regression tests、三 replica、续跑、无 NaN/爆炸和轨迹解析 |
| github-mims-harvard-primekg | 历史知识图谱背景 | background_only | 上游已声明被 OptimusKG 取代 | 不进入运行依赖或评分 | 只保留迁移事实和历史来源边界 |
| github-mims-harvard-optimuskg | 大型统一知识图谱 | reject | 当前闭环要求逐来源版本化和最小下载，直接图谱不优于核心 API | 不安装、不下载 | not_run |

## 实际反馈

| knowledge_id | outcome | 实际验证 | 可复用边界 |
| --- | --- | --- | --- |
| github-ccsb-scripps-autodock-vina | supported | 三个 native seed 的对称校正 RMSD 为 0.274、0.277、0.280 Å，全部通过 2 Å 门 | 只支持当前 1IEP 设置的 pose 恢复，不代表亲和力或筛选性能 |
| github-forlilab-meeko | supported | 受体、三类配体及九组 pose 的准备/导出成功，SDF 可由 RDKit 复读 | 新化学类型仍须重新审计模板、质子化和键级 |
| github-maabuu-posebusters | supported | native、阳性对照和 decoy 的九组输出均有通过全部布尔门的 pose | 必须与 native RMSD、对照和输入化学审计联合使用 |
| gitlab-gromacs-gromacs | supported | 2026.3 源构建 106/106 测试通过；39,969 原子体系完成三 replica、checkpoint 续跑、无 NaN/爆炸和轨迹解析 | 只支持 CPU 短闭环；GPU production 为 `resource-limited`，1 ps 不支持稳定性或机制结论 |
| github-mims-harvard-primekg | not_run | 未安装、未下载、未进入评分 | 继续只作历史背景 |
| github-mims-harvard-optimuskg | not_run | 未安装、未下载、未进入评分 | 当前种子驱动任务继续拒绝引入完整图谱 |
