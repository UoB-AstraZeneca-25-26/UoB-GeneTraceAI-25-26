import {
  RankedCellLine,
  ResultMeta,
  ScoreTier,
  SelectivityRanked,
  SingleRanked,
} from '../types';

const BASE = 'https://9368clqa34.execute-api.eu-north-1.amazonaws.com';
const GENE_URL = `${BASE}/gene`;
const EXCLUDE_URL = `${BASE}/exclude`;

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

// ---------- fetchers ----------

async function getJson<T>(url: string, signal?: AbortSignal): Promise<T> {
  const res = await fetch(url, { signal });
  if (!res.ok) throw new Error(`API error: ${res.status} ${res.statusText}`);
  return (await res.json()) as T;
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

/**
 * High-level entry point used by the app. Routes to the selectivity endpoint
 * when an exclusion is supplied, otherwise the single-gene endpoint.
 */
export async function fetchRanking(
  target: string,
  exclusion: string | undefined,
  signal?: AbortSignal
): Promise<{ results: RankedCellLine[]; meta: ResultMeta }> {
  if (exclusion && exclusion.trim()) {
    const resp = await fetchSelectivityRanking(target, exclusion, signal);
    return { results: toSelectivityRanked(resp), meta: metaFromExclude(resp) };
  }
  const resp = await fetchGeneRanking(target, signal);
  return { results: toSingleRanked(resp), meta: metaFromGene(resp) };
}