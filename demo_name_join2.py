import pandas as pd

sample_info = pd.read_parquet(r'C:\Disertation\UoB-GeneTraceAI-25-26\data\parquet\9_DepMap_sample_info.parquet')
cellosaurus = pd.read_parquet(r'C:\Disertation\UoB-GeneTraceAI-25-26\data\parquet\7_cellosaurus.parquet')

cello_named = cellosaurus[['Identifier (cell line name)','Accession (CVCL_xxxx)']].rename(
    columns={'Identifier (cell line name)': 'cell_line_name',
             'Accession (CVCL_xxxx)': 'cvcl_accession'})

si = sample_info[['DepMap_ID','cell_line_name','RRID']].copy()

# ── Name join (exact on cell_line_name) ──────────────────────────────
name_join = si.dropna(subset=['cell_line_name']).merge(
    cello_named, on='cell_line_name', how='left')

# ── RRID join ────────────────────────────────────────────────────────
rrid_join = si.merge(
    cellosaurus[['Accession (CVCL_xxxx)','Synonyms']].rename(
        columns={'Accession (CVCL_xxxx)':'RRID'}),
    on='RRID', how='left')

print('=== Match rate comparison ===')
print(f'Name join  — rows fed in: {len(name_join)}, matched: {name_join["cvcl_accession"].notna().sum()}, missed: {name_join["cvcl_accession"].isna().sum()} ({name_join["cvcl_accession"].isna().sum()/len(name_join)*100:.1f}%)')
print(f'RRID join  — rows fed in: {len(rrid_join)}, matched: {rrid_join["RRID"].notna().sum()},   missed: {rrid_join["RRID"].isna().sum()} ({rrid_join["RRID"].isna().sum()/len(rrid_join)*100:.1f}%)')

# ── Show what the name join misses but RRID gets ─────────────────────
name_missed = name_join[name_join['cvcl_accession'].isna()]['DepMap_ID']
rrid_hit    = rrid_join[rrid_join['RRID'].notna()]['DepMap_ID']
missed_by_name_but_rrid_works = set(name_missed) & set(rrid_hit)

print(f'\nCell lines name-join MISSES but RRID-join GETS: {len(missed_by_name_but_rrid_works)}')
examples = sample_info[sample_info['DepMap_ID'].isin(list(missed_by_name_but_rrid_works)[:8])][['DepMap_ID','cell_line_name','RRID']]
print(examples.to_string(index=False))

# ── Show the synonym problem: same cell line, many names in real data ─
print('\n=== The synonym problem: one RRID, many names in the wild ===')
# Pull synonyms for a well-known cell line
example_rrid = 'CVCL_0031'  # MCF-7
row = cellosaurus[cellosaurus['Accession (CVCL_xxxx)'] == example_rrid]
if len(row):
    print(f'RRID: {example_rrid}')
    print(f'Cellosaurus name : {row["Identifier (cell line name)"].values[0]}')
    print(f'Synonyms         : {row["Synonyms"].values[0]}')
    print()

example_rrid2 = 'CVCL_0089'  # MHH-CALL-3 from ACH-000032 used in earlier example
row2 = cellosaurus[cellosaurus['Accession (CVCL_xxxx)'] == example_rrid2]
if len(row2):
    print(f'RRID: {example_rrid2}')
    print(f'Cellosaurus name : {row2["Identifier (cell line name)"].values[0]}')
    print(f'Synonyms         : {row2["Synonyms"].values[0]}')

print()
print('=== The duplicate name problem: U-251 MG in sample_info ===')
u251 = sample_info[sample_info['cell_line_name'] == 'U-251 MG'][['DepMap_ID','cell_line_name','RRID','primary_disease','lineage']]
print(u251.to_string(index=False))
print()
print('Same name, different DepMap_IDs, one has no RRID.')
print('A name-based join cannot tell them apart.')
print('An RRID-based join routes each to the correct Cellosaurus entry (or null).')