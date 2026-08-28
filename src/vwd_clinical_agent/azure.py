from __future__ import annotations

from dotenv import load_dotenv
import httpx
import json
import os
from typing import Any, Protocol


class LLMProvider(Protocol):
    name: str
    version: str

    def complete_json(self, system_prompt: str, user_prompt: str) -> dict[str, Any]:
        """Return a validated JSON object."""


class DeterministicLLMProvider:
    name = "deterministic_policy"
    version = "v0"

    def complete_json(self, system_prompt: str, user_prompt: str) -> dict[str, Any]:
        try:
            state = json.loads(user_prompt)
        except json.JSONDecodeError:
            state = {}
        evidence = state.get("evidence", [])
        conflicts = state.get("evidence_conflicts", [])
        missing = state.get("evidence_missing", [])
        actions = state.get("recommended_actions", [])
        variants = state.get("variants", [])
        summary = (
            f"Retrospective review integrated {len(evidence)} FHIR evidence resources across "
            f"{len(variants)} reported variant(s), identified {len(conflicts)} structured evidence issue(s), "
            f"and found {len(missing)} missing evidence item(s). "
        )
        if actions:
            summary += "The evidence-constrained second-level plan prioritizes " + ", ".join(actions) + ". "
        if state.get("second_level_status") == "not_available":
            summary += (
                "All requested second-level tests were explicitly unavailable in this retrospective environment; "
                "no values were imputed. "
            )
        summary += "The system abstains from autonomous diagnosis and requires expert review."
        return {
            "summary": summary,
            "abstention": True,
            "expert_review_required": True,
            "candidate_subtypes": [],
        }


class AzureOpenAIProvider:
    name = "azure_openai"
    version = "chat-completions-rest-v1"

    def __init__(
        self,
        endpoint: str,
        api_key: str,
        deployment: str,
        api_version: str = "2024-10-21",
        timeout: float = 30.0,
        use_json_mode: bool = True,
    ):
        self.endpoint = endpoint.rstrip("/")
        self.api_key = api_key
        self.deployment = deployment
        self.api_version = api_version
        self.timeout = timeout
        self.use_json_mode = use_json_mode

    @staticmethod
    def _normalize_deployment(value: str) -> str:
        return value.removeprefix("azure/").strip()

    @classmethod
    def from_env(cls) -> "AzureOpenAIProvider":
        load_dotenv()
        endpoint = os.getenv("AZURE_OPENAI_ENDPOINT")
        api_key = os.getenv("AZURE_OPENAI_API_KEY")
        deployment = os.getenv("AZURE_OPENAI_DEPLOYMENT") or os.getenv("SUBAGENT_MODEL_NAME")
        if not endpoint or not api_key or not deployment:
            raise RuntimeError(
                "Azure OpenAI requires AZURE_OPENAI_ENDPOINT, AZURE_OPENAI_API_KEY, "
                "and AZURE_OPENAI_DEPLOYMENT or SUBAGENT_MODEL_NAME."
            )
        return cls(
            endpoint=endpoint,
            api_key=api_key,
            deployment=cls._normalize_deployment(deployment),
            api_version=os.getenv("AZURE_OPENAI_API_VERSION", "2024-10-21"),
            timeout=float(os.getenv("AZURE_OPENAI_TIMEOUT", "30")),
            use_json_mode=os.getenv("AZURE_OPENAI_JSON_MODE", "1") not in {"0", "false", "False"},
        )

    def complete_json(self, system_prompt: str, user_prompt: str) -> dict[str, Any]:
        payload: dict[str, Any] = {
            "messages": [
                {"role": "system", "content": system_prompt},
                {"role": "user", "content": user_prompt},
            ],
            "temperature": 0,
            "max_completion_tokens": 1200,
        }
        if self.use_json_mode:
            payload["response_format"] = {"type": "json_object"}
        url = (
            f"{self.endpoint}/openai/deployments/{self.deployment}/chat/completions"
            f"?api-version={self.api_version}"
        )
        response = httpx.post(
            url,
            headers={"api-key": self.api_key, "Content-Type": "application/json"},
            json=payload,
            timeout=self.timeout,
        )
        if response.status_code >= 400:
            raise RuntimeError(
                f"Azure OpenAI request failed with HTTP {response.status_code}: {response.text}"
            )
        content = response.json()["choices"][0]["message"]["content"]
        try:
            parsed = json.loads(content)
        except json.JSONDecodeError as exc:
            raise ValueError("Azure OpenAI returned non-JSON content") from exc
        if not isinstance(parsed, dict):
            raise ValueError("Azure OpenAI JSON response must be an object")
        return parsed
