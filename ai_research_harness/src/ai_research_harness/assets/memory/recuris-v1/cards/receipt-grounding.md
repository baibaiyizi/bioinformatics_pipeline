---
skill_id: receipt-grounded-completion
version: 1
status: active
scope: global
component: E
domains: [engineering, reproduction, agent, workflow, evaluation]
triggers: [完成, done, artifact, 结果, output, command, 测试, test, run, reproduce]
anti_triggers: [pure brainstorming]
sources: [arxiv-2608.24876, github-gen-verse-recuris]
evidence_requirements: [authoritative_completed_event, checker_result, artifact_hash]
---

# 用 receipt 提交完成状态

模型只能提出完成。Harness 必须收到权威 completed event，并由对应 checker 验证真实文件、hash、测试语义或目标输出后，才能把原子目标从 `pending` 提交为 `done`。命令退出码 0 只证明进程返回成功；模型文字、合成 receipt、文件存在或测试骨架均不能单独证明研究目标完成。

失败、拒绝、权限不足和输出为空必须保留为 `blocked` 或未完成，并记录唯一阻塞原因和下一步，不生成替代结果。
