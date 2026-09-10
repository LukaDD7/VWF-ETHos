#!/usr/bin/env python3
"""Verify generated FHIR files and optionally require the official HL7 validator."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import re
import subprocess
import sys
from typing import Any


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from src.vwd_clinical_agent.mechanistic_tasks import ProtocolRegistry  # noqa: E402


def _emit(payload: Any, *, stream: Any = sys.stdout) -> None:
    print(json.dumps(payload, ensure_ascii=False, indent=2), file=stream)


def _expected(registry: ProtocolRegistry) -> dict[str, dict[str, Any]]:
    expected = {
        "PlanDefinition-vwd-mechanism-guided-acquisition-v1.json": (
            registry.fhir_plan_definition()
        )
    }
    expected.update(
        {
            f"ActivityDefinition-{protocol_id}.json": registry.fhir_activity_definition(protocol_id)
            for protocol_id in registry.protocols
        }
    )
    return expected


def validate_export(registry: ProtocolRegistry, fhir_dir: Path) -> list[Path]:
    files: list[Path] = []
    for name, expected in _expected(registry).items():
        path = fhir_dir / name
        actual = json.loads(path.read_text(encoding="utf-8"))
        if actual != expected:
            raise ValueError(f"FHIR export is stale or differs from registry: {path}")
        if actual.get("resourceType") not in {"PlanDefinition", "ActivityDefinition"}:
            raise ValueError(f"Unexpected FHIR resource type in {path}")
        for required in ("id", "url", "version", "status"):
            if not actual.get(required):
                raise ValueError(f"FHIR resource {path} lacks required field {required}")
        files.append(path)
    return files


def validate_with_hl7_jar(files: list[Path], validator_jar: Path) -> list[dict[str, Any]]:
    if not validator_jar.is_file():
        raise FileNotFoundError(f"HL7 validator jar does not exist: {validator_jar}")
    reports: list[dict[str, Any]] = []
    for path in files:
        completed = subprocess.run(
            [
                "java",
                "-jar",
                str(validator_jar),
                str(path),
                "-version",
                "5.0.0",
                "-tx",
                "n/a",
            ],
            cwd=ROOT,
            check=False,
            capture_output=True,
            text=True,
        )
        output = "\n".join(item for item in (completed.stdout, completed.stderr) if item)
        error_counts = [
            int(match)
            for match in re.findall(r"\bErrors?\s*[:=]\s*(\d+)", output, flags=re.IGNORECASE)
        ]
        if completed.returncode != 0 or any(count > 0 for count in error_counts):
            raise ValueError(
                f"Official HL7 FHIR R5 validation failed for {path}: {output[-3000:]}"
            )
        reports.append(
            {
                "file": str(path),
                "returncode": completed.returncode,
                "reported_error_counts": error_counts,
            }
        )
    return reports


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--registry",
        type=Path,
        default=ROOT / "protocols/vwd_mechanistic_v1/registry.json",
    )
    parser.add_argument(
        "--fhir-dir",
        type=Path,
        default=ROOT / "protocols/vwd_mechanistic_v1/fhir",
    )
    parser.add_argument("--validator-jar", type=Path)
    parser.add_argument(
        "--require-official",
        action="store_true",
        help="Fail unless --validator-jar is supplied and official validation passes.",
    )
    args = parser.parse_args()
    try:
        registry = ProtocolRegistry.load(args.registry)
        files = validate_export(registry, args.fhir_dir)
        if args.require_official and args.validator_jar is None:
            raise ValueError("--require-official requires --validator-jar")
        official_reports = (
            validate_with_hl7_jar(files, args.validator_jar)
            if args.validator_jar is not None
            else []
        )
        _emit(
            {
                "ok": True,
                "registry_digest": registry.digest,
                "files": [str(path) for path in files],
                "registry_export_parity": True,
                "official_hl7_fhir_r5_validated": bool(official_reports),
                "official_reports": official_reports,
            }
        )
        return 0
    except Exception as exc:
        _emit(
            {"ok": False, "error": {"type": type(exc).__name__, "message": str(exc)}},
            stream=sys.stderr,
        )
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
