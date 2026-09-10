import sys, time, warnings, importlib.util
warnings.filterwarnings('ignore')
spec = importlib.util.spec_from_file_location('m2', r'C:\Disertation\UoB-GeneTraceAI-25-26\architecture\02_Proteinomics\06_protein_score_M2.py')
m2 = importlib.util.module_from_spec(spec)
spec.loader.exec_module(m2)
import duckdb

t0 = time.time()
con = duckdb.connect(r'C:\Disertation\UoB-GeneTraceAI-25-26\architecture\outputs\celllineselector.db', read_only=True)
print('connected', time.time()-t0, flush=True)

result, isoform_qc, thr = m2.compute_protein_z(con, table_names=[m2.CCLE_SOURCE, m2.PROCAN_SOURCE],
                                                alpha=m2.DEFAULT_ALPHA, null_pctile=m2.DEFAULT_NULL_PCTILE,
                                                matched_pctile=m2.DEFAULT_MATCHED_PCTILE, min_gap=m2.DEFAULT_MIN_GAP,
                                                seed=0, verbose=True)
print('compute_protein_z done', time.time()-t0, 's', flush=True)
print('result shape:', result.shape, flush=True)

out_dir = r'C:\Disertation\UoB-GeneTraceAI-25-26\architecture\02_Proteinomics\results'
result.to_parquet(out_dir + r'\protein_z_m2.parquet', index=False)
isoform_qc.to_parquet(out_dir + r'\protein_isoform_qc_m2.parquet', index=False)
import json
with open(out_dir + r'\m2_thresholds.json', 'w') as f:
    json.dump(thr, f, indent=2)
print('TOTAL', time.time()-t0, 's', flush=True)
con.close()
