FROM public.ecr.aws/lambda/python:3.12

WORKDIR ${LAMBDA_TASK_ROOT}

COPY final_pipeline/requirements.txt .
RUN pip install --no-cache-dir -r requirements.txt

COPY final_pipeline/ ./final_pipeline/

# Data files — large; override via env vars if mounting EFS instead
COPY final_pipeline/outputs/predictions_with_confidence.parquet ./final_pipeline/outputs/predictions_with_confidence.parquet
COPY final_pipeline/reference/gene_lookup.parquet              ./final_pipeline/reference/gene_lookup.parquet
COPY final_pipeline/reference/cell_line_lookup.parquet         ./final_pipeline/reference/cell_line_lookup.parquet
COPY final_pipeline/outputs/celllineselector.db                ./final_pipeline/outputs/celllineselector.db

CMD ["final_pipeline.api.handler.handler"]
