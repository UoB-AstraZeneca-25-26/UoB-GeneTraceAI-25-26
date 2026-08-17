"""Canonical entrypoint (closes gap G7).

    python main.py                 # dev server
    uvicorn api.app:app            # production (see Dockerfile)
"""

from api.app import app
from api.settings import get_settings

if __name__ == "__main__":
    import uvicorn

    settings = get_settings()
    uvicorn.run(app, host=settings.host, port=settings.port)
