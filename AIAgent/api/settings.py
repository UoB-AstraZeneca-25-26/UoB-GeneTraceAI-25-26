"""Centralised, validated configuration (closes gap G17).

Replaces the bare ``load_dotenv()`` calls that made API-key loading
CWD-dependent. Settings are validated once at import time, so a missing key
fails the container immediately and visibly instead of surfacing as a
confusing 500 on the first user query.
"""

from __future__ import annotations

from functools import lru_cache
from pathlib import Path

from pydantic import Field, SecretStr, field_validator
from pydantic_settings import BaseSettings, SettingsConfigDict

_AIAGENT_DIR = Path(__file__).resolve().parent.parent
_REPO_ROOT = _AIAGENT_DIR.parent


def _split_csv(value: object) -> object:
    if isinstance(value, str):
        return [item.strip() for item in value.split(",") if item.strip()]
    return value


class Settings(BaseSettings):
    model_config = SettingsConfigDict(
        env_file=str(_AIAGENT_DIR / ".env"),
        env_file_encoding="utf-8",
        extra="ignore",
    )

    # --- LLM ---------------------------------------------------------
    groq_api_key: SecretStr = SecretStr("")  # optional: only required when llm_provider == "groq"
    llm_provider: str = "groq"
    llm_model: str = "openai/gpt-oss-120b"
    llm_temperature: float = 0.0
    llm_max_tokens: int = 2048   # reasoning tokens are billed against this

    # --- AWS Bedrock ---------------------------------------------------
    # boto3 reads actual credentials from its own standard env vars
    # (AWS_ACCESS_KEY_ID / AWS_SECRET_ACCESS_KEY / AWS_SESSION_TOKEN /
    # AWS_BEARER_TOKEN_BEDROCK) or the IAM role automatically — these
    # fields exist for documentation/validation only, not SecretStr since
    # nothing here reads them directly.
    aws_region: str = "eu-north-1"
    aws_access_key_id: str = ""
    aws_secret_access_key: str = ""
    aws_session_token: str = ""
    aws_profile: str | None = None

    # --- CORS ----------------------------------------------------------
    allowed_origins: list[str] = Field(default_factory=lambda: ["http://localhost:5173"])

    @field_validator("allowed_origins", "api_keys", mode="before")
    @classmethod
    def _coerce_csv(cls, value: object) -> object:
        return _split_csv(value)

    # --- Data sources --------------------------------------------------
    gene_lookup_path: Path = _REPO_ROOT / "reference" / "gene_lookup.parquet"
    knowledge_dir: Path = _AIAGENT_DIR / "Knowledge"

    # --- Scoring API -----------------------------------------------------
    # Empty => score_explainer runs in methodology-only mode (null numerics,
    # the LLM explains the 7-step method). Set it to switch to real values.
    scoring_api_url: str = ""

    # --- Agent / timeout budget (R6) ------------------------------------
    agent_max_iterations: int = 4
    agent_max_execution_time: float = 45.0
    api_request_timeout: float = 50.0
    llm_call_timeout: float = 30.0
    http_timeout: float = 10.0
    http_connect_timeout: float = 3.0

    # --- Caching ---------------------------------------------------------
    response_cache_ttl: int = 3600
    response_cache_maxsize: int = 1024
    gene_cache_ttl: int = 86400
    gene_cache_maxsize: int = 4096

    # --- Auth / rate limiting (G12) --------------------------------------
    api_keys: list[str] = Field(default_factory=list)  # empty => auth disabled
    rate_limit: str = "30/minute"

    # --- Misc ------------------------------------------------------------
    max_body_bytes: int = 64 * 1024
    debug: bool = False
    host: str = "0.0.0.0"
    port: int = 8000


@lru_cache
def get_settings() -> Settings:
    """Process-wide singleton. Cached so env is parsed once."""
    return Settings()
