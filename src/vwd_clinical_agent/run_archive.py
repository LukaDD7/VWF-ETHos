from __future__ import annotations

from datetime import datetime, timezone
import json
from pathlib import Path
from typing import Any

from pydantic import BaseModel


def _jsonable(value: Any) -> Any:
    if isinstance(value, BaseModel):
        return value.model_dump(mode="json")
    if isinstance(value, Path):
        return str(value)
    if isinstance(value, datetime):
        return value.isoformat()
    if isinstance(value, list):
        return [_jsonable(item) for item in value]
    if isinstance(value, dict):
        return {key: _jsonable(item) for key, item in value.items()}
    return value


class RunArchive:
    """Complete on-disk archive for one LangGraph run."""

    def __init__(self, run_id: str, archive_root: str | Path):
        self.run_id = run_id
        self.root = Path(archive_root) / run_id
        self.root.mkdir(parents=True, exist_ok=True)
        self.cases_dir = self.root / "cases"
        self.cases_dir.mkdir(exist_ok=True)
        self.checkpoint_db = self.root / "checkpoints.sqlite"
        self.manifest_path = self.root / "run_manifest.json"

    def case_dir(self, case_id: str) -> Path:
        path = self.cases_dir / case_id
        path.mkdir(parents=True, exist_ok=True)
        return path

    def write_debug_event(self, case_id: str, event: dict[str, Any]) -> None:
        path = self.case_dir(case_id) / "debug_trace.jsonl"
        with path.open("a", encoding="utf-8") as handle:
            handle.write(json.dumps(_jsonable(event), ensure_ascii=False, default=str) + "\n")

    def write_state_history(self, case_id: str, history: list[dict[str, Any]]) -> None:
        path = self.case_dir(case_id) / "state_history.jsonl"
        with path.open("w", encoding="utf-8") as handle:
            for snapshot in history:
                handle.write(json.dumps(_jsonable(snapshot), ensure_ascii=False, default=str) + "\n")

    def write_final_state(self, case_id: str, state: dict[str, Any]) -> None:
        path = self.case_dir(case_id) / "final_state.json"
        path.write_text(json.dumps(_jsonable(state), ensure_ascii=False, indent=2, default=str), encoding="utf-8")

    def write_report(self, case_id: str, final_opinion: Any) -> None:
        data = _jsonable(final_opinion)
        json_path = self.case_dir(case_id) / "final_report.json"
        json_path.write_text(json.dumps(data, ensure_ascii=False, indent=2), encoding="utf-8")
        markdown = self._report_markdown(case_id, data)
        (self.case_dir(case_id) / "final_report.md").write_text(markdown, encoding="utf-8")

    def write_manifest(self, metadata: dict[str, Any]) -> None:
        payload = {
            "run_id": self.run_id,
            "created_at": datetime.now(timezone.utc).isoformat(),
            **metadata,
        }
        self.manifest_path.write_text(json.dumps(_jsonable(payload), ensure_ascii=False, indent=2), encoding="utf-8")

    @staticmethod
    def _report_markdown(case_id: str, data: dict[str, Any]) -> str:
        lines = [
            f"# VWD Research Report — {case_id}",
            "",
            f"- Status: {data.get('confidence', 'low')}",
            f"- Abstention: {data.get('abstention')}",
            f"- Expert review required: {data.get('expert_review_required')}",
            "",
            "## Opinion",
            "",
            str(data.get("opinion", "")),
            "",
            "## Candidate subtypes",
            "",
            ", ".join(data.get("candidate_subtypes", [])) or "N/A",
            "",
            "## Supporting evidence",
            "",
        ]
        lines.extend(f"- {ref}" for ref in data.get("supporting_evidence", []))
        lines.extend(["", "## Missing information", ""])
        lines.extend(f"- {item}" for item in data.get("missing_information", []))
        lines.extend(["", "## Limitations", ""])
        lines.extend(f"- {item}" for item in data.get("limitations", []))
        return "\n".join(lines) + "\n"
