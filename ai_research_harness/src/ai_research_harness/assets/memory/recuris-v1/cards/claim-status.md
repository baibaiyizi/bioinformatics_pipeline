---
skill_id: claim-status-boundary
version: 1
status: active
scope: global
component: E
domains: [evidence, research, reproduction, paper, code]
triggers: [论文, paper, 代码, code, 复现, reproduction, claim, evidence, 评论]
anti_triggers: [creative fiction]
sources: [arxiv-2608.24876, github-gen-verse-recuris]
evidence_requirements: [source_locator, version_or_commit, claim_status]
---

# 区分来源主张、实现与复现

先给每个关键陈述分配状态：想法、来源作者主张、文档声称、源码中存在、当前真实执行结果或独立复现。作者论文和 README 只能支持各自声称；固定 commit 的源码审计只支持实现观察；当前环境中的命令、日志和 artifact 才支持实际执行。独立复现还要求独立环境、固定协议和匹配输入，不能由作者公告、stars、issue 或重复运行安装检查替代。

输出时将三条证据 lane 分开：论文原文与 locator；代码路径、commit 和测试；第三方评论及其独立性。没有证据就保留 `unresolved`，不要升级措辞。
