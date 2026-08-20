"""
final_pipeline/api/bedrock.py
------------------------------
Native Amazon Bedrock Agent action-group adapter.

Bedrock invokes the Lambda directly (by ARN — stable across every redeploy),
passing its own event shape rather than an API Gateway HTTP event. This module
parses that event, calls the same ranking functions the HTTP routes use, and
returns the response envelope Bedrock expects.

Wiring in Bedrock: create an action group whose executor is this Lambda's ARN,
with an OpenAPI schema exposing /gene, /genes, /exclude, /gene/detail, /health.
No API Gateway, no URL to drift.
"""

from __future__ import annotations

import json
import logging

from .ranking import detail, ensure_loaded, exclude, is_ready, rank

logger = logging.getLogger(__name__)


def is_bedrock_event(event: dict) -> bool:
    """Bedrock agent events carry messageVersion + an agent block."""
    return isinstance(event, dict) and "messageVersion" in event and "agent" in event


def _params(event: dict) -> dict:
    """Flatten Bedrock's parameter list [{name,type,value}, ...] into a dict.

    Also merges any application/json requestBody properties, so the same
    handler works whether the agent sends query params or a body.
    """
    out: dict = {}
    for p in event.get("parameters", []) or []:
        out[p["name"]] = p["value"]

    body = (event.get("requestBody", {}) or {}).get("content", {})
    props = (body.get("application/json", {}) or {}).get("properties", []) or []
    for p in props:
        out[p["name"]] = p["value"]
    return out


def _lineage_list(raw) -> list[str]:
    if not raw:
        return []
    if isinstance(raw, list):
        return [str(x).strip() for x in raw if str(x).strip()]
    return [t.strip() for t in str(raw).split(",") if t.strip()]


def _dispatch(api_path: str, params: dict) -> dict:
    if api_path == "/health":
        return {"status": "ok"}

    if api_path == "/gene":
        gene = params["gene"]
        return rank([gene], _lineage_list(params.get("lineage")),
                    float(params.get("floor", 0.0)),
                    int(params.get("top_n", 30)))

    if api_path == "/genes":
        genes = [g.strip() for g in str(params["genes"]).split(",") if g.strip()]
        if len(genes) < 2:
            raise ValueError("Provide at least 2 comma-separated genes.")
        return rank(genes, _lineage_list(params.get("lineage")),
                    float(params.get("floor", 0.5)),
                    int(params.get("top_n", 30)))

    if api_path == "/exclude":
        return exclude(params["gene_a"], params["gene_b"],
                       _lineage_list(params.get("lineage")),
                       int(params.get("top_n", 30)))

    if api_path == "/gene/detail":
        return detail(params["gene"], params["model_id"])

    raise ValueError(f"Unknown apiPath: {api_path!r}")


# Maps the console "function details" function names to the OpenAPI-style paths,
# so both action-group definition styles hit the same dispatch logic.
_FUNCTION_TO_PATH = {
    "health": "/health",
    "gene": "/gene",
    "genes": "/genes",
    "exclude": "/exclude",
    "detail": "/gene/detail",
    "gene_detail": "/gene/detail",
}


def handle_bedrock(event: dict) -> dict:
    """Run a Bedrock action-group request and return its response envelope.

    Supports both action-group definition styles:
      * OpenAPI schema   -> event has apiPath / httpMethod
      * function details -> event has function (mapped to a path here)
    """
    action_group = event.get("actionGroup", "")
    function = event.get("function")
    http_method = event.get("httpMethod", "GET")

    if function is not None:  # console "function details" style
        api_path = _FUNCTION_TO_PATH.get(function, "")
    else:                     # OpenAPI schema style
        api_path = event.get("apiPath", "")

    ensure_loaded()
    if not is_ready():
        status, payload = 503, {"error": "Ranking data not loaded."}
    else:
        try:
            payload = _dispatch(api_path, _params(event))
            status = 200
        except (ValueError, KeyError) as exc:
            status, payload = 422, {"error": str(exc)}
        except Exception as exc:  # noqa: BLE001 - surface any failure to the agent
            logger.exception("bedrock dispatch failed")
            status, payload = 500, {"error": str(exc)}

    body = json.dumps(payload)

    if function is not None:  # function-details response envelope
        return {
            "messageVersion": "1.0",
            "response": {
                "actionGroup": action_group,
                "function": function,
                "functionResponse": {
                    "responseBody": {"TEXT": {"body": body}}
                },
            },
        }

    return {  # OpenAPI response envelope
        "messageVersion": "1.0",
        "response": {
            "actionGroup": action_group,
            "apiPath": api_path,
            "httpMethod": http_method,
            "httpStatusCode": status,
            "responseBody": {
                "application/json": {"body": body}
            },
        },
    }
