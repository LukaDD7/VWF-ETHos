from __future__ import annotations

from abc import ABC, abstractmethod
from hashlib import sha256
from pathlib import Path
import json
import time
from typing import Any, Literal, Protocol

from pydantic import BaseModel, Field, ValidationError

from .fhir import FHIRBundle, FHIRResource, operation_outcome, provenance, utc_now


TOOL_POLICY_PATH = Path(__file__).with_name("tool_policy.json")


class ToolRequest(BaseModel):
    operation: str
    patient_id: str
    variant_id: str
    parameters: dict[str, Any] = Field(default_factory=dict)
    input: FHIRBundle = Field(default_factory=FHIRBundle)


class ToolResponse(BaseModel):
    tool: str
    operation: str
    status: Literal["success", "not_found", "error", "disabled"]
    output: FHIRBundle
    diagnostics: list[str] = Field(default_factory=list)
    retry_count: int = 0
    latency_ms: int = 0
    cache_hit: bool = False


class ToolContract(Protocol):
    name: str
    version: str

    def invoke(self, request: ToolRequest) -> ToolResponse:
        ...


class BaseBiomedicalTool(ABC):
    name: str = "base_tool"
    version: str = "0"
    endpoint: str = "local://"
    timeout: float = 30.0
    max_retries: int = 2

    @abstractmethod
    def run(self, request: ToolRequest) -> tuple[list[FHIRResource], str]:
        """Return FHIR resources and a status."""

    def invoke(self, request: ToolRequest) -> ToolResponse:
        started = time.monotonic()
        try:
            resources, status = self.run(request)
            diagnostics: list[str] = []
        except NotImplementedError as exc:
            resources, status, diagnostics = [operation_outcome("warning", "not-supported", str(exc))], "disabled", [str(exc)]
        except Exception as exc:
            resources, status, diagnostics = [operation_outcome("error", "exception", str(exc))], "error", [str(exc)]
        bundle = FHIRBundle.of(resources)
        elapsed = int((time.monotonic() - started) * 1000)
        return ToolResponse(
            tool=self.name,
            operation=request.operation,
            status=status,  # type: ignore[arg-type]
            output=bundle,
            diagnostics=diagnostics,
            latency_ms=elapsed,
        )

    def provenance_for(self, resource_id: str, request: ToolRequest, response: Any) -> FHIRResource:
        return provenance(
            resource_id=resource_id,
            tool_name=self.name,
            endpoint=self.endpoint,
            request_payload=request.model_dump(mode="json"),
            response_payload=response,
            occurred=utc_now(),
        )


class ToolRegistry:
    def __init__(self, tools: list[BaseBiomedicalTool] | None = None, cache_dir: str | Path | None = None):
        self.tools: dict[str, BaseBiomedicalTool] = {tool.name: tool for tool in tools or []}
        self.cache_dir = Path(cache_dir) if cache_dir else None
        if self.cache_dir:
            self.cache_dir.mkdir(parents=True, exist_ok=True)

    def register(self, tool: BaseBiomedicalTool) -> None:
        self.tools[tool.name] = tool

    @staticmethod
    def _policy_violation(tool_name: str, operation: str) -> str | None:
        policy = json.loads(TOOL_POLICY_PATH.read_text(encoding="utf-8"))
        allowed_operations = policy["allowed_tools"].get(tool_name)
        if allowed_operations is None:
            return f"Tool is not registered in code-as-search policy: {tool_name}"
        if operation not in allowed_operations:
            return f"Operation {operation!r} is not allowed for tool {tool_name!r}"
        return None

    def invoke(self, tool_name: str, request: ToolRequest) -> ToolResponse:
        tool = self.tools.get(tool_name)
        violation = self._policy_violation(tool_name, request.operation)
        if violation:
            return ToolResponse(
                tool=tool_name,
                operation=request.operation,
                status="error",
                output=FHIRBundle.of([operation_outcome("error", "business-rule", violation)]),
                diagnostics=[violation],
            )
        if tool is None:
            return ToolResponse(
                tool=tool_name,
                operation=request.operation,
                status="error",
                output=FHIRBundle.of([operation_outcome("error", "not-found", f"Unknown tool: {tool_name}")]),
                diagnostics=[f"Unknown tool: {tool_name}"],
            )
        key = sha256(f"{tool.name}|{request.model_dump_json()}".encode("utf-8")).hexdigest()
        cache_path = self.cache_dir / f"{key}.json" if self.cache_dir else None
        if cache_path and cache_path.exists():
            try:
                cached = ToolResponse.model_validate_json(cache_path.read_text(encoding="utf-8"))
                cached.cache_hit = True
                return cached
            except ValidationError:
                cache_path.unlink(missing_ok=True)
        response = tool.invoke(request)
        if response.status in {"success", "not_found", "disabled"} and cache_path:
            cache_path.write_text(response.model_dump_json(indent=2), encoding="utf-8")
        return response
