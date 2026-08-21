FROM python:3.12-slim

COPY --from=public.ecr.aws/awsguru/aws-lambda-adapter:0.9.1 /lambda-adapter /opt/extensions/lambda-adapter
ENV PORT=8000

WORKDIR /app

COPY AIAgent/requirements.lock.txt .
RUN pip install --no-cache-dir -r requirements.lock.txt
COPY AIAgent/ ./AIAgent/
COPY reference/gene_lookup.parquet ./reference/gene_lookup.parquet
RUN useradd -m app && chown -R app /app
USER app
WORKDIR /app/AIAgent
CMD ["uvicorn", "api.app:app", "--host", "0.0.0.0", "--port", "8000"]