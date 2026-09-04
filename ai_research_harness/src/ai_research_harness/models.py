from __future__ import annotations

from datetime import datetime, timezone
from enum import StrEnum
from pathlib import Path
from typing import Any, Literal

from pydantic import BaseModel, ConfigDict, Field, field_validator, model_validator


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


class GoalStatus(StrEnum):
    pending = "pending"
    done = "done"
    blocked = "blocked"


class CheckResult(BaseModel):
    model_config = ConfigDict(extra="forbid")

    checker_id: str
    passed: bool
    detail: str
    evidence: dict[str, Any] = Field(default_factory=dict)


class Receipt(BaseModel):
    model_config = ConfigDict(extra="forbid")

    receipt_id: str
    kind: Literal[
        "agent_message",
        "command_execution",
        "file_change",
        "mcp_tool_call",
        "web_search",
        "artifact",
        "turn",
    ]
    authoritative: bool
    status: Literal["completed", "failed", "declined", "interrupted"]
    item_id: str | None = None
    artifact_path: str | None = None
    sha256: str | None = None
    details: dict[str, Any] = Field(default_factory=dict)
    created_at: str = Field(default_factory=utc_now)

    @model_validator(mode="after")
    def artifact_has_digest(self) -> "Receipt":
        if self.kind == "artifact" and self.status == "completed":
            if not self.artifact_path or not self.sha256:
                raise ValueError("completed artifact receipts require path and sha256")
        return self


class Goal(BaseModel):
    model_config = ConfigDict(extra="forbid")

    goal_id: str
    objective: str
    acceptance_checks: list[str]
    status: GoalStatus = GoalStatus.pending
    receipts: list[Receipt] = Field(default_factory=list)
    check_results: list[CheckResult] = Field(default_factory=list)
    blocker: str | None = None

    @field_validator("acceptance_checks")
    @classmethod
    def checks_not_empty(cls, value: list[str]) -> list[str]:
        if not value:
            raise ValueError("a goal requires at least one acceptance checker")
        if len(value) != len(set(value)):
            raise ValueError("acceptance checker IDs must be unique")
        return value

    def propose_done(self, receipts: list[Receipt], checks: list[CheckResult]) -> None:
        if self.status != GoalStatus.pending:
            raise ValueError("only pending goals can become done")
        completed = [r for r in receipts if r.authoritative and r.status == "completed"]
        by_id = {check.checker_id: check for check in checks}
        missing = [checker for checker in self.acceptance_checks if checker not in by_id]
        failed = [checker for checker in self.acceptance_checks if checker in by_id and not by_id[checker].passed]
        if not completed:
            raise ValueError("done requires an authoritative completed receipt")
        if missing or failed:
            raise ValueError(f"acceptance checks incomplete: missing={missing}, failed={failed}")
        self.receipts.extend(receipts)
        self.check_results.extend(checks)
        self.status = GoalStatus.done
        self.blocker = None

    def block(self, reason: str, receipts: list[Receipt] | None = None) -> None:
        if self.status != GoalStatus.pending:
            raise ValueError("only pending goals can become blocked")
        if not reason.strip():
            raise ValueError("blocked goals require a reason")
        self.status = GoalStatus.blocked
        self.blocker = reason.strip()
        if receipts:
            self.receipts.extend(receipts)


class WorkingMemory(BaseModel):
    model_config = ConfigDict(extra="forbid")

    schema_version: Literal[1] = 1
    run_id: str
    task_id: str
    project_root: str
    phase: str = "claim_contract"
    goals: list[Goal]
    active_constraints: list[str] = Field(default_factory=list)
    selected_skill_ids: list[str] = Field(default_factory=list)
    thread_id: str | None = None
    turn_ids: list[str] = Field(default_factory=list)
    created_at: str = Field(default_factory=utc_now)
    updated_at: str = Field(default_factory=utc_now)

    @field_validator("project_root")
    @classmethod
    def project_is_absolute(cls, value: str) -> str:
        if not Path(value).is_absolute():
            raise ValueError("project_root must be absolute")
        return value

    @model_validator(mode="after")
    def unique_goal_ids(self) -> "WorkingMemory":
        ids = [goal.goal_id for goal in self.goals]
        if len(ids) != len(set(ids)):
            raise ValueError("goal IDs must be unique")
        return self

    def pending_goals(self) -> list[Goal]:
        return [goal for goal in self.goals if goal.status == GoalStatus.pending]

    def touch(self) -> None:
        self.updated_at = utc_now()


class SkillCard(BaseModel):
    model_config = ConfigDict(extra="forbid")

    skill_id: str
    version: int = Field(ge=1)
    status: Literal["active", "deprecated"] = "active"
    scope: Literal["global", "project"]
    component: Literal["E"] = "E"
    domains: list[str]
    triggers: list[str]
    anti_triggers: list[str] = Field(default_factory=list)
    sources: list[str]
    evidence_requirements: list[str]
    path: Path
    body: str
    token_estimate: int = Field(ge=1)


class MemoryManifest(BaseModel):
    model_config = ConfigDict(extra="forbid")

    schema_version: Literal[1] = 1
    version_id: str
    scope: Literal["global", "project"]
    parent_version: str | None = None
    cards: list[str]
    working_memory_spec: str
    invocation_policy: str
    checker_config: str
    source_failure_ids: list[str] = Field(default_factory=list)
    created_at: str


class TraceRecord(BaseModel):
    model_config = ConfigDict(extra="forbid")

    timestamp: str = Field(default_factory=utc_now)
    run_id: str
    event: str
    thread_id: str | None = None
    turn_id: str | None = None
    item_id: str | None = None
    state_before_sha256: str | None = None
    state_after_sha256: str | None = None
    selected_skill_ids: list[str] = Field(default_factory=list)
    receipt: dict[str, Any] | None = None
    checker_results: list[dict[str, Any]] = Field(default_factory=list)
    commit_decision: str | None = None
    details: dict[str, Any] = Field(default_factory=dict)


class QualificationTask(BaseModel):
    model_config = ConfigDict(extra="forbid")

    task_id: str
    domain: Literal["evidence", "empirical_ml", "reproducibility", "systems", "theory", "security", "hci"]
    critical: bool = False
    prompt: str
    required_groups: list[list[str]]
    forbidden: list[str] = Field(default_factory=list)
    gold_fragments: list[str] = Field(default_factory=list)
    semantic_review: bool = False
    rubric: str = ""


class ArmOutcome(BaseModel):
    model_config = ConfigDict(extra="forbid")

    task_id: str
    arm: Literal["incumbent", "candidate"]
    passed: bool
    deterministic_passed: bool
    semantic_passed: bool | None = None
    response: str
    reasons: list[str] = Field(default_factory=list)
    usage: dict[str, Any] = Field(default_factory=dict)
    duration_ms: int | None = None
    trace: dict[str, Any] = Field(default_factory=dict)


class GateDecision(BaseModel):
    model_config = ConfigDict(extra="forbid")

    candidate_id: str
    incumbent_id: str
    status: Literal["PROMOTED", "NOT_ADMITTED"]
    source_failures_repaired: bool
    critical_regressions: list[str]
    domain_deltas: dict[str, int]
    wins: int
    losses: int
    ties: int
    net_wins: int
    bootstrap_ci95: tuple[float, float]
    mcnemar_exact_p: float
    metrics: dict[str, Any] = Field(default_factory=dict)
    checks: list[CheckResult]
    created_at: str = Field(default_factory=utc_now)
