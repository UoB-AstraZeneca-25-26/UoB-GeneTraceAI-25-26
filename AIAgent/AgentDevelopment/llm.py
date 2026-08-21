"""LLM backend configuration for the CellLineFinder agent.

Swap provider by changing env vars — no code change required.
"""

import os

from dotenv import load_dotenv
from langchain_core.language_models.chat_models import BaseChatModel

load_dotenv()


def get_llm() -> BaseChatModel:
    """Return a configured chat model based on env vars.

    Reads env fresh on every call — no global state.
    """
    provider = os.getenv("LLM_PROVIDER", "groq").lower()
    model = os.getenv("LLM_MODEL", "qwen-qwq-32b")
    temperature = float(os.getenv("LLM_TEMPERATURE", "0.0"))
    # R6 timeout budget: the LLM call is the innermost layer and must be
    # strictly shorter than the agent's max_execution_time (45s).
    timeout = float(os.getenv("LLM_CALL_TIMEOUT", "30.0"))
    # 2048, not 512: the default model (openai/gpt-oss-120b) is a reasoning
    # model whose thinking tokens are billed against max_tokens. A 512 budget
    # was spent entirely on reasoning (finish_reason="length", 0 content), so
    # every answer came back empty. 2048 leaves room for ~1100 reasoning
    # tokens plus the 7-step score_explainer JSON, and still fits the 8K TPM
    # limit alongside a ~3.3K prompt.
    max_tokens = int(os.getenv("LLM_MAX_TOKENS", "2048"))
    # Optional, and provider-specific: on a reasoning model "low" cuts ~1000
    # thinking tokens (and ~1s) per call. Left unset by default because a
    # non-reasoning model rejects the parameter outright.
    reasoning_effort = os.getenv("LLM_REASONING_EFFORT", "").strip()

    if provider == "groq":
        from langchain_groq import ChatGroq

        extra = {"reasoning_effort": reasoning_effort} if reasoning_effort else {}
        return ChatGroq(
            model=model,
            temperature=temperature,
            api_key=os.getenv("GROQ_API_KEY"),
            max_tokens=max_tokens,
            timeout=timeout,
            **extra,
        )

    if provider == "bedrock":
        from langchain_aws import ChatBedrockConverse

        kwargs: dict = {
            "model_id": model,
            "region_name": os.getenv("AWS_REGION", "eu-north-1"),
            "temperature": temperature,
            "max_tokens": max_tokens,
        }

        # Only set profile if explicitly configured — otherwise boto3's
        # standard credential chain (env vars, bearer token, IAM role)
        # resolves credentials on its own.
        profile = os.getenv("AWS_PROFILE")
        if profile:
            kwargs["credentials_profile_name"] = profile

        return ChatBedrockConverse(**kwargs)

    if provider == "huggingface":
        raise NotImplementedError("huggingface provider is a future stub")

    if provider == "openai":
        raise NotImplementedError("openai provider is a future stub")

    raise ValueError(f"Unknown LLM_PROVIDER: {provider!r}")
