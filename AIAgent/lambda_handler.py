import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from mangum import Mangum
from api.app import app

# lifespan="auto" runs the app's startup (httpx client, gene index, LLM
# prewarm -- see api/app.py) once per cold start and keeps it for every
# warm invocation that reuses this container.
handler = Mangum(app, lifespan="auto")