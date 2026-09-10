"""
Lambda entry point.

One function image serves two callers, auto-detected per invocation:
  * Amazon Bedrock Agent  — invokes this Lambda directly by ARN (stable across
    every redeploy). Handled natively, no API Gateway involved.
  * API Gateway (HTTP)    — the website's REST calls, wrapped by Mangum/FastAPI.

Lambda handler:  final_pipeline.api.handler.handler
"""

from mangum import Mangum

from .app import app
from .bedrock import handle_bedrock, is_bedrock_event

# API Gateway serves this under the /prod stage, so the event path arrives as
# "/prod/genes". Strip that stage prefix before Starlette routing, otherwise
# every route 404s as "/prod/genes" instead of "/genes".
_http_handler = Mangum(app, lifespan="on", api_gateway_base_path="/prod")


def handler(event, context):
    if is_bedrock_event(event):
        return handle_bedrock(event)
    return _http_handler(event, context)
