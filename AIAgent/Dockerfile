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

CMD ["uvicorn", "api.app:app", "--host", "0.0.0.0", "--port", "8000", "--workers", "2"]
