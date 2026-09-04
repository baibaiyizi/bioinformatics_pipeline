---
skill_id: matched-agent-evaluation-budget
version: 1
status: active
scope: global
component: E
domains: [empirical_ml, llm, agent, benchmark, evaluation, systems]
triggers: [agent, LLM, 模型, benchmark, 评测, 对照, baseline, retry, token, step, budget]
anti_triggers: [single descriptive report]
sources: [arxiv-2608.24876]
evidence_requirements: [matched_arms, frozen_split, cost_and_retry_report]
---

# 匹配 Agent 评测预算

比较任何 LLM 或 Agent 方法前，冻结模型及版本、系统与开发指令、工具、权限、上下文、数据 split、step、retry、token、wall-time、并发和硬件/费用。两个 arm 只允许目标组件不同。把额外 retry 产生的收益单列，不能写成记忆、技能或学习本身的收益。

报告每任务成对结果、失败率、资源使用、效应量和区间。开发集选择、in-sample 路由和最好 seed 必须显式标注；小样本不以不显著证明等价，也不以大点估计冒充稳定泛化。
