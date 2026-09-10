"""LLM backend — Amazon Bedrock (Claude Haiku 4.5).

Provider is fixed to Bedrock. Swap the model via LLM_MODEL env var.
Auth is handled by the IAM role attached to the host (no API key needed).
"""

import os

from langchain_aws import ChatBedrock
from langchain_core.language_models.chat_models import BaseChatModel


def get_llm() -> BaseChatModel:
    """Return a configured Bedrock chat model."""
    return ChatBedrock(
        model_id=os.getenv("LLM_MODEL", "anthropic.claude-haiku-4-5"),
        model_kwargs={
            "temperature": float(os.getenv("LLM_TEMPERATURE", "0.0")),
            "max_tokens": int(os.getenv("LLM_MAX_TOKENS", "1024")),
        },
        region_name=os.getenv("AWS_REGION", "us-east-1"),
    )
