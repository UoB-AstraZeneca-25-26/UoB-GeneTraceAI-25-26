import {
  CellLineDetail,
  MetadataItem,
  MultiRanked,
  RankedCellLine,
  ResultMeta,
  ScoreTier,
  SelectivityRanked,
  SingleRanked,
  TrackScores,
} from '../types';

const BASE = 'https://9368clqa34.execute-api.eu-north-1.amazonaws.com';
const GENE_URL = `${BASE}/gene`;
const EXCLUDE_URL = `${BASE}/exclude`;
const GENES_URL = `${BASE}/genes`;

// ---------- /gene ----------

export interface GeneApiLine {
  rank: number;
  model_id: string;
  name: string;
  score: number;
  tier: string;
  n_layers: number;
  driver_alteration: boolean;
}

export interface GeneApiResponse {
  gene: string;
  ensg: string;
  total: number;
  showing: number;
  lines: GeneApiLine[];
}

// ---------- /prod/exclude ----------
// score_<gene> keys are dynamic — read them via gene_high / gene_low.

export interface ExcludeApiLine {
  model_id: string;
  name: string;
  selectivity: number;
  [scoreKey: `score_${string}`]: number;
}

export interface ExcludeApiResponse {
  gene_high: string;
  gene_low: string;
  formula: string;
  total_ranked: number;
  lines: ExcludeApiLine[];
}

// ---------- /genes (joint multi-gene) ----------
// Per-gene scores may be NaN (no evidence). "None" names fall back to model_id.

export interface GenesApiLine {
  model_id: string;
  name: string;
  joint_score: number;
  scores: Record<string, number | null>;
}

export interface GenesApiResponse {
  genes: string[];
  floor: number;
  total_passing: number;
  lines: GenesApiLine[];
}

// ---------- /gene?...&cell_line=<id> (per-line detail) ----------

export interface CellLineDetailApiResponse {
  gene: string;
  ensg: string;
  cell_line: { model_id: string; name: string };
  rank: number;
  total: number;
  score: number;
  tier: string;
  n_layers: number;
  driver_alteration: boolean;
  p_mutation: number | null;
  p_fusion: number | null;
  has_cna_alteration: boolean;
  metadata: Record<string, string | number | null>;
  rna_alternatives: { model_id: string; name: string; similarity: number }[];
}

// ---------- fetchers ----------

/* The /genes endpoint emits bare NaN/Infinity, which are invalid JSON and make
   JSON.parse throw — losing the whole payload, not just the bad cell. This repairs
   those tokens to null, but ONLY in value positions (preceded by : [ , and followed
   by , } ]) so quoted strings are never touched. The proper fix is server-side. */
function parseLenientJson<T>(text: string): T {
  const repaired = text.replace(/(?<=[:[,]\s*)(?:-?Infinity|NaN)(?=\s*[,}\]])/g, 'null');
  return JSON.parse(repaired) as T;
}

async function getJson<T>(url: string, signal?: AbortSignal): Promise<T> {
  const res = await fetch(url, { signal });
  if (!res.ok) throw new Error(`API error: ${res.status} ${res.statusText}`);
  return (await res.json()) as T;
}

async function getJsonLenient<T>(url: string, signal?: AbortSignal): Promise<T> {
  const res = await fetch(url, { signal });
  if (!res.ok) throw new Error(`API error: ${res.status} ${res.statusText}`);
  return parseLenientJson<T>(await res.text());
}

export async function fetchGeneRanking(
  gene: string,
  signal?: AbortSignal
): Promise<GeneApiResponse> {
  const symbol = gene.trim();
  if (!symbol) throw new Error('No gene symbol provided.');
  const data = await getJson<GeneApiResponse>(
    `${GENE_URL}?query=${encodeURIComponent(symbol)}`,
    signal
  );
  if (!data || !Array.isArray(data.lines)) throw new Error('Unexpected /gene response shape.');
  return data;
}

export async function fetchSelectivityRanking(
  geneA: string,
  geneB: string,
  signal?: AbortSignal
): Promise<ExcludeApiResponse> {
  const a = geneA.trim();
  const b = geneB.trim();
  if (!a || !b) throw new Error('Both a target and an exclusion gene are required.');
  const data = await getJson<ExcludeApiResponse>(
    `${EXCLUDE_URL}?gene_a=${encodeURIComponent(a)}&gene_b=${encodeURIComponent(b)}`,
    signal
  );
  if (!data || !Array.isArray(data.lines)) throw new Error('Unexpected /exclude response shape.');
  return data;
}

export async function fetchMultiRanking(
  genes: string[],
  signal?: AbortSignal
): Promise<GenesApiResponse> {
  const clean = genes.map((g) => g.trim()).filter(Boolean);
  if (clean.length < 2) throw new Error('Provide at least two target genes for a joint ranking.');
  // Encode each symbol but keep commas literal to match the endpoint's format.
  const query = clean.map(encodeURIComponent).join(',');
  const data = await getJsonLenient<GenesApiResponse>(`${GENES_URL}?query=${query}`, signal);
  if (!data || !Array.isArray(data.lines)) throw new Error('Unexpected /genes response shape.');
  return data;
}

// ---------- adapters ----------

function normalizeTier(raw: string): ScoreTier {
  const t = (raw || '').toUpperCase();
  return t === 'HIGH' || t === 'MEDIUM' || t === 'LOW' ? t : 'LOW';
}

/** Spread near-identical scores across 5..100% for a readable bar. */
function relativize(values: number[], value: number): number {
  const min = Math.min(...values);
  const max = Math.max(...values);
  const span = max - min;
  return span > 0 ? Math.round((5 + ((value - min) / span) * 95) * 10) / 10 : 100;
}

export function toSingleRanked(resp: GeneApiResponse): SingleRanked[] {
  const scores = resp.lines.map((l) => l.score);
  return resp.lines.map((l) => ({
    mode: 'single',
    rank: l.rank,
    modelId: l.model_id,
    cellLine: l.name.toUpperCase(),
    score: l.score,
    confidenceScore: Math.round(l.score * 1000) / 10,
    relativeScore: relativize(scores, l.score),
    tier: normalizeTier(l.tier),
    nLayers: l.n_layers,
    driverAlteration: l.driver_alteration,
  }));
}

export function toSelectivityRanked(resp: ExcludeApiResponse): SelectivityRanked[] {
  const highKey: `score_${string}` = `score_${resp.gene_high}`;
  const lowKey: `score_${string}` = `score_${resp.gene_low}`;
  const sels = resp.lines.map((l) => l.selectivity);
  return resp.lines.map((l, i) => ({
    mode: 'selectivity',
    rank: i + 1, // no rank field; lines arrive pre-sorted by selectivity
    modelId: l.model_id,
    cellLine: l.name.toUpperCase(),
    selectivity: l.selectivity,
    scoreHigh: Number(l[highKey] ?? 0),
    scoreLow: Number(l[lowKey] ?? 0),
    geneHigh: resp.gene_high,
    geneLow: resp.gene_low,
    relativeScore: relativize(sels, l.selectivity),
  }));
}

export function metaFromGene(resp: GeneApiResponse): ResultMeta {
  return {
    mode: 'single',
    primaryGene: resp.gene,
    ensg: resp.ensg,
    total: resp.total,
    showing: resp.showing,
  };
}

export function metaFromExclude(resp: ExcludeApiResponse): ResultMeta {
  return {
    mode: 'selectivity',
    primaryGene: resp.gene_high,
    excludedGene: resp.gene_low,
    formula: resp.formula,
    total: resp.total_ranked,
    showing: resp.lines.length,
  };
}

export function toMultiRanked(resp: GenesApiResponse): MultiRanked[] {
  const joints = resp.lines.map((l) => l.joint_score);
  return resp.lines.map((l, i) => {
    const geneScores: Record<string, number | null> = {};
    for (const g of resp.genes) {
      const v = l.scores?.[g];
      // NaN survives as null after lenient parse; guard anyway.
      geneScores[g] = typeof v === 'number' && Number.isFinite(v) ? v : null;
    }
    const name = l.name && l.name !== 'None' ? l.name.toUpperCase() : l.model_id;
    return {
      mode: 'multi',
      rank: i + 1, // no rank field; lines arrive pre-sorted by joint_score
      modelId: l.model_id,
      cellLine: name,
      jointScore: l.joint_score,
      genes: resp.genes,
      geneScores,
      relativeScore: relativize(joints, l.joint_score),
    };
  });
}

export function metaFromGenes(resp: GenesApiResponse): ResultMeta {
  return {
    mode: 'multi',
    primaryGene: resp.genes.join(' + '),
    genes: resp.genes,
    floor: resp.floor,
    total: resp.total_passing,
    showing: resp.lines.length,
  };
}

/**
 * High-level entry point. Routes by what's selected:
 *   2+ targets            -> /genes  (joint multi-gene ranking; exclusions ignored)
 *   1 target + exclusion  -> /exclude (selectivity)
 *   1 target              -> /gene   (single)
 */
export async function fetchRanking(
  targets: string[],
  exclusion: string | undefined,
  signal?: AbortSignal
): Promise<{ results: RankedCellLine[]; meta: ResultMeta }> {
  const clean = targets.map((t) => t.trim()).filter(Boolean);
  if (clean.length === 0) throw new Error('Select at least one target gene.');

  if (clean.length >= 2) {
    const resp = await fetchMultiRanking(clean, signal);
    return { results: toMultiRanked(resp), meta: metaFromGenes(resp) };
  }

  if (exclusion && exclusion.trim()) {
    const resp = await fetchSelectivityRanking(clean[0], exclusion, signal);
    return { results: toSelectivityRanked(resp), meta: metaFromExclude(resp) };
  }

  const resp = await fetchGeneRanking(clean[0], signal);
  return { results: toSingleRanked(resp), meta: metaFromGene(resp) };
}

// ---------- cell-line detail ----------

export async function fetchCellLineDetail(
  gene: string,
  modelId: string,
  signal?: AbortSignal
): Promise<CellLineDetailApiResponse> {
  const g = gene.trim();
  const id = modelId.trim();
  if (!g || !id) throw new Error('Both a gene and a cell-line model id are required.');
  const data = await getJson<CellLineDetailApiResponse>(
    `${GENE_URL}?query=${encodeURIComponent(g)}&cell_line=${encodeURIComponent(id)}`,
    signal
  );
  if (!data || !data.cell_line) throw new Error('Unexpected cell-line detail response shape.');
  return data;
}

// Human-readable labels for the metadata block, in display order. Keys not in
// this map are skipped so the panel stays clean and predictable.
const METADATA_LABELS: [string, string][] = [
  ['lineage', 'Lineage'],
  ['lineage_subtype', 'Lineage subtype'],
  ['primary_disease', 'Primary disease'],
  ['subtype', 'Disease subtype'],
  ['cellosaurus_ncit_disease', 'NCIt disease'],
  ['primary_or_metastasis', 'Origin'],
  ['sample_collection_site', 'Collection site'],
  ['default_growth_pattern', 'Growth pattern'],
  ['sex', 'Sex'],
  ['age', 'Age'],
];

function cleanMetadata(meta: Record<string, string | number | null>): MetadataItem[] {
  const items: MetadataItem[] = [];
  for (const [key, label] of METADATA_LABELS) {
    const raw = meta?.[key];
    if (raw === null || raw === undefined || raw === '' || raw === 'None') continue;
    items.push({ label, value: String(raw) });
  }
  return items;
}

export function toCellLineDetail(resp: CellLineDetailApiResponse): CellLineDetail {
  // Only the layers this endpoint actually reports get a track. Expression,
  // proteomics and signature aren't returned per-line, so they stay absent
  // rather than being shown as "no signal".
  const tracks: TrackScores = {};
  if (typeof resp.p_mutation === 'number' && Number.isFinite(resp.p_mutation)) {
    tracks.mutation = {
      present: true,
      score: resp.p_mutation,
      finding: resp.driver_alteration ? 'driver' : undefined,
    };
  }
  if (typeof resp.p_fusion === 'number' && Number.isFinite(resp.p_fusion)) {
    tracks.fusion = { present: true, score: resp.p_fusion };
  }
  tracks.cna = {
    present: resp.has_cna_alteration,
    score: resp.has_cna_alteration ? 1 : 0,
    finding: resp.has_cna_alteration ? 'altered' : 'none',
  };

  const name = resp.cell_line.name && resp.cell_line.name !== 'None'
    ? resp.cell_line.name.toUpperCase()
    : resp.cell_line.model_id;

  return {
    gene: resp.gene,
    ensg: resp.ensg,
    modelId: resp.cell_line.model_id,
    cellLine: name,
    rank: resp.rank,
    total: resp.total,
    score: resp.score,
    tier: resp.tier,
    nLayers: resp.n_layers,
    driverAlteration: resp.driver_alteration,
    pMutation: resp.p_mutation,
    pFusion: resp.p_fusion,
    hasCnaAlteration: resp.has_cna_alteration,
    tracks,
    metadata: cleanMetadata(resp.metadata),
    alternatives: (resp.rna_alternatives || []).map((a) => ({
      modelId: a.model_id,
      name: a.name && a.name !== 'None' ? a.name.toUpperCase() : a.model_id,
      similarity: a.similarity,
    })),
  };
}