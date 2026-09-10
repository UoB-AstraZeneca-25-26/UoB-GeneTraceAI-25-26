FROM python:3.12-slim
WORKDIR /app
ENV PYTHONUNBUFFERED=1 PIP_NO_CACHE_DIR=1

COPY AIAgent/requirements.lock.txt .
RUN pip install --no-cache-dir -r requirements.lock.txt

COPY AIAgent/ ./AIAgent/
COPY reference/gene_lookup.parquet ./reference/gene_lookup.parquet

RUN useradd -m app && chown -R app /app
USER app

WORKDIR /app/AIAgent
EXPOSE 8000

HEALTHCHECK --interval=30s --timeout=3s CMD python -c "import urllib.request;urllib.request.urlopen('http://localhost:8000/health')"

# --workers 1, not 2: Lambda runs exactly one invocation at a time per
# container instance (concurrency comes from separate instances, not
# in-process workers), so a second uvicorn worker buys nothing here except
# importing the entire langchain/langgraph/pandas/pyarrow stack a second
# time. That doubled import cost is what pushed cold start past Lambda's
# 10s init-phase budget after this lockfile picked up much newer (heavier)
# langchain/langgraph/langchain-aws releases than whatever unpinned
# versions happened to resolve back when this last worked -- the container
# then fell through to finishing startup during the first real invocation,
# consuming its 60s timeout instead and killing every request until fixed.
CMD ["uvicorn", "api.app:app", "--host", "0.0.0.0", "--port", "8000", "--workers", "1"]
