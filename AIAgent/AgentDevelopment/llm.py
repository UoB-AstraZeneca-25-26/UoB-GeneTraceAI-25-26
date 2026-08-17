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

    if provider == "groq":
        from langchain_groq import ChatGroq

        return ChatGroq(
            model=model,
            temperature=temperature,
            api_key=os.getenv("GROQ_API_KEY"),
            max_tokens=512,
        )

    if provider == "huggingface":
        raise NotImplementedError("huggingface provider is a future stub")

    if provider == "openai":
        raise NotImplementedError("openai provider is a future stub")

    raise ValueError(f"Unknown LLM_PROVIDER: {provider!r}")
