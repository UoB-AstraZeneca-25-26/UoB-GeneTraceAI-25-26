FROM public.ecr.aws/lambda/python:3.12

WORKDIR ${LAMBDA_TASK_ROOT}

COPY final_pipeline/requirements.txt .
RUN pip install --no-cache-dir -r requirements.txt

COPY final_pipeline/ ./final_pipeline/

# No data files are baked in -- everything (predictions, gene/cell/sample_info
# lookups, per-source omics parquets) is fetched from S3 at cold start (see
# ranking.py's _resolve_data_path, DATA_BUCKET env var, s3://<bucket>/processed/).
# The image never depends on what happens to exist on the machine running
# `docker build` -- a fresh git clone with no data files builds an identical
# image. .dockerignore also excludes any parquet/db that might be present
# locally, so this holds even if someone forgets to remove local data first.

CMD ["final_pipeline.api.handler.handler"]
