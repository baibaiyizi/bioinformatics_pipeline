---
skill_id: component-local-memory-admission
version: 1
status: active
scope: global
component: E
domains: [agent, harness, memory, evaluation, governance]
triggers: [memory, 记忆, skill, 技能, 演化, evolve, gate, failure, 失败, regression]
anti_triggers: [ordinary content editing]
sources: [arxiv-2608.24876, github-gen-verse-recuris, url-seongyeon1-recuris-review]
evidence_requirements: [failure_localization, isolated_candidate, frozen_regression_gate]
---

# 组件定位后才允许记忆晋升

先把失败定位到经验技能 E、工作记忆 W、调用策略 rho、checker C，或不可由记忆修补的 model/tool/data/environment/evaluation。一次候选只改责任组件，并始终与 incumbent 隔离。候选需通过 schema、hash、容量、近重复、gold leakage、来源失败修复和冻结回归门后才可晋升。

未通过 gate 的候选不得成为 provisional base。自动晋升只操作声明式内容，使用不可变版本、原子 champion 指针和追加 ledger；可执行代码、模型、权限、预算、冻结题和 gate 阈值不属于自动演化范围。
