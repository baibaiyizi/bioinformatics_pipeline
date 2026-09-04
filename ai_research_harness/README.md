# Codex AI Research Harness

这是 `ai_research` 的可审计 E/W/rho/C 运行入口。普通任务可直接调用 Codex 的 `ai_research` agent；需要 working memory、receipt、恢复、失败定位、冻结 gate 和自动 champion 晋升时使用本目录。

## 安装

```bash
cd /home/h1028/workspace/model/ai_research_harness
uv sync --extra test
./run_ai_research.sh doctor
```

环境固定为 Python 3.12–3.13 与 `openai-codex==0.147.0`。SDK 使用其自带的稳定 Codex runtime 和现有 Codex 登录；指定的 `gpt-5.6-sol` 不可用时 `doctor` 返回 4，不切换模型。

Harness 仅在自己的 App Server 子进程中禁用继承的可选 MCP、插件、web search、通知和 Codex 内建 memory，避免未授权网络、外部工具改变冻结预算或被本机的跨平台路径阻塞；不会修改 `~/.codex/config.toml`。Python SDK 的实验接口被显式关闭。

## 入口

```bash
./run_ai_research.sh init --project /absolute/project
./run_ai_research.sh run --project /absolute/project --prompt "核验这篇论文和代码"
./run_ai_research.sh resume --project /absolute/project --run-id run-xxxx
./run_ai_research.sh status --project /absolute/project --run-id run-xxxx
./run_ai_research.sh evolve --project /absolute/project --run-id run-xxxx --scope project
./run_ai_research.sh gate --candidate recuris-v1 --scope global
./run_ai_research.sh rollback --version base-v1 --scope global
```

退出码：0 完成或 gate 通过；2 阻塞/需要输入；3 gate 未通过；4 环境或协议错误。

## 状态与权限

- 项目运行写入 `<project>/.ai_research/`；全局 memory 默认在 `~/.codex/ai_research_memory/`，可用 `AI_RESEARCH_MEMORY_ROOT` 在测试或隔离环境中重定向。
- SDK 固定使用 `ApprovalMode.deny_all`。研究任务只使用 read-only 或项目 workspace-write；权限升级会阻塞。
- 模型只能提出状态。Harness 仅从 `item/completed` 形成 receipt，并在 checker 通过后提交 `done`。
- 自动演化仅允许声明式 Markdown/YAML；候选不能包含 Python/Shell，不能更改模型、权限、预算、资格题或 gate 阈值。
- 全局候选运行冻结 32 项资格集；`init` 会把同一资格集及 SHA-256 建为项目初始回归锚点。若要加入项目特异题，必须在生成候选前冻结 `.ai_research/regression_anchors/qualification.yaml` 并同步更新 `qualification.sha256`；缺失或 hash 不匹配时返回 `NOT_ADMITTED`。

## 结果解释

`PROMOTED` 表示 candidate 满足工程资格门，不表示统计学普适提升、论文结果复现或外部独立验证。`NOT_ADMITTED` 保留 incumbent；未通过 gate 的 candidate 不参与后续正式任务。

当前内置 `recuris-v1` 的首次冻结 gate 已完成，结果为 `NOT_ADMITTED`；详情见 [VALIDATION_REPORT.md](VALIDATION_REPORT.md)。资格运行固定为单次回答一个 turn、零工具动作、零任务内重试、100,000 token 上限、120 秒上限与最多三个并发 SDK turn。报告同时保存匿名 reviewer、每次 trial、token/时间成本、memory 增量和 task-level 配对统计。
