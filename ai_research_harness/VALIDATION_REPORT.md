# AI Research Harness 验收报告

检查日期：2026-09-03。

## 实现状态

- 12 个 `ai_research` TOML 可解析且名称唯一；共享契约存在。
- Python 环境锁定为 3.12–3.13，当前解析为 Python 3.13；`openai-codex` 与 bundled runtime 均为 `0.147.0`。
- `doctor` 的版本、路径、agent、契约、`ai_learning`、32 题 schema/hash、memory anchor、认证与 `gpt-5.6-sol` 检查全部通过。
- wheel 与 sdist 构建通过，wheel 包含全部 Python 模块、冻结资格集和三套 memory bundle。

## 离线与协议测试

- Harness：39/39 pytest 通过；另通过 `compileall` 和两个 Bash 入口语法检查。
- 覆盖状态迁移、权威 receipt、真实 artifact 重哈希、成功命令但 checker 失败、completed/failed/declined/interrupted、部分文件、候选隔离、近重复、gold leakage、hash 篡改、原子晋升、显式 rollback、context compaction、中断后从追加 trace 恢复 thread、sandbox 越权和声明式演化限制。
- 最终代码 live SDK smoke：`run-acc9622ff5b2` 完成；thread=`01a06561-9a11-7012-9fd5-111768eeaed7`，turn=`01a06561-9d15-7321-8ff1-5476efc4b014`。收到 2 个权威 receipt（`item/completed`、`turn/completed`），checker 全过，W 原子提交为 `done`，phase boundary 的 rho 选择记录存在。

## 32 项冻结成对 gate

- 决策：`NOT_ADMITTED`；全局 champion 仍为 `base-v1`。
- 结果：candidate 1 胜、1 负、30 平，净胜 0；task-level bootstrap 95% CI `[-0.09375, 0.09375]`；McNemar exact `p=1.0`。
- domain delta：evidence `+1`、empirical_ml `-1`，其余领域 `0`；无关键项回归，但“各领域不下降”与“净胜至少 1”未通过。
- 两个声明的来源失败没有同时满足“incumbent 失败且 candidate 通过”：`evidence-03` 两者均通过，`empirical-03` 两者均被词法 checker 判失败。
- 成本：incumbent 753,706 tokens / 645,236 ms；candidate 751,159 tokens / 579,155 ms；candidate memory 增加 4 卡、4,542 bytes。
- 固定资格集 SHA-256：`8423c2e5946fefcfdd8deb1d52df671b069ecd3ae0165077c12eacbff294a988`。
- 决策、outcomes 与 84 条成对 trial trace：`/home/h1028/.codex/ai_research_memory/reports/20260903T033757Z-recuris-v1-vs-base-v1/`。
- 两个 arm 在同一 `x86_64` Linux 主机运行：2 路 Intel Xeon Platinum 8173M、112 logical CPUs、Python 3.13.2；本 gate 不使用本地 GPU 计算。

冻结后审计发现 deterministic lexical rubric 存在否定语境和同义词假阴性，例如“不能称为显著提升”仍命中 forbidden phrase，“预算不等”未命中“不匹配”。不追改 v1 资格题，不重跑择优；该问题记录为下一版评测设计的 `evaluation` failure。因此本结果既不能写成候选改进，也不是对 Recuris 思想的普适否定。

配对 bootstrap 使用固定 seed `260824876`。`openai-codex==0.147.0` 的稳定 turn API 不暴露模型采样 seed，因此模型 seed 无法冻结；协议以相同模型、指令、输入、工具、sandbox、预算和顺序交换约束两 arm，并对首次不一致题重复两次取多数。没有为追求 seed 控制启用实验接口。

## 现有 SCI 回归

既有 SCI smoke 的 48 个测试均在其已有依赖环境中通过：通用/ML/single-cell multiome/空间代谢/网络药理使用当前工作区 Python，分子建模 11 项使用既有 `sci_molecular_modeling` conda 环境。直接以 base Python 跑总脚本时，唯一环境错误是 base 缺少 RDKit；未修改 SCI Agent 或其依赖配置。

## 最终边界

双入口、E/W/rho/C、恢复、演化隔离和冻结 gate 已真实运行；`recuris-v1` 仍只在 candidate 目录，未进入正式任务。只有未来新的预注册候选实际通过同一门控，才能写 `PROMOTED/IMPROVED`。
