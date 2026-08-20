"""
Lambda entry point — wraps the FastAPI app with Mangum.

AWS Lambda invokes this as:  handler.handler
"""

from mangum import Mangum
from .app import app

# API Gateway serves this under the /prod stage, so the event path arrives as
# "/prod/genes". Strip that stage prefix before Starlette routing, otherwise
# every route 404s as "/prod/genes" instead of "/genes".
handler = Mangum(app, lifespan="on", api_gateway_base_path="/prod")
