import {
  mockGeneResponse,
  mockExcludeResponse,
  mockGenesResponse,
  mockCellLineDetailResponse,
  mockDelay,
} from './mockApi';
import {
  CellLineDetail,
  MetadataItem,
  MultiRanked,
  RankedCellLine,
  ResultMeta,
  SelectivityRanked,
  SingleRanked,
  TrackScores,
} from '../types';

// In dev, use the same-origin '/api' path so the Vite proxy handles the request
// (sidestepping the missing CORS headers). In a production build, call the API
// directly — which requires the API to send CORS headers, or same-origin hosting.
const BASE = import.meta.env.DEV
  ? '/api'
  : 'https://9368clqa34.execute-api.eu-north-1.amazonaws.com';

// The whole API lives under the /prod stage.
const API = `${BASE}/prod`;
const GENE_URL = `${API}/gene`;
const GENES_URL = `${API}/genes`;
const EXCLUDE_URL = `${API}/exclude`;
const DETAIL_URL = `${API}/gene/detail`;

// How many lines to request from the server (client slices further for display).
const TOP_N = 30;

// ---- Offline mock mode -------------------------------------------------------
// Flip to false to use the live API. When true, all fetchers return local mock
// data in the same shapes as the live endpoints.
export const USE_MOCK = false;

// ---------- shared line metadata ----------
export interface LineMetadata {
  lineage?: string;
  lineage_subtype?: string;
  primary_disease?: string;
  sex?: string;
}

// ---------- /prod/gene and /prod/genes (joint shape) ----------
// Single-gene and multi-gene ranking now return the same shape. Per-gene scores
// may be NaN (no evidence); "None"/null names fall back to model_id.
export interface JointApiLine {
  model_id: string;
  name: string | null;
  joint_score: number;
  limiting_gene?: string;
  scores: Record<string, number | null>;
  metadata?: LineMetadata;
}

export interface JointApiResponse {
  genes: string[];
  lineage?: string[];
  floor: number;
  total_passing: number;
  lines: JointApiLine[];
}

// ---------- /prod/exclude (selectivity) ----------
export interface ExcludeApiLine {
  model_id: string;
  name: string | null;
  score_a: number;
  score_b: number;
  selectivity: number;
  metadata?: LineMetadata;
}

export interface ExcludeApiResponse {
  gene_a: string;
  gene_b: string;
  lineage?: string[];
  total_ranked: number;
  lines: ExcludeApiLine[];
}

// ---------- /prod/gene/detail (per-line detail) ----------
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

// ---------- fetch helpers ----------

/* Some endpoints can emit bare NaN/Infinity, which are invalid JSON and make
   JSON.parse throw. This repairs those tokens to null, only in value positions,
   so quoted strings are never touched. The proper fix is server-side. */
function parseLenientJson<T>(text: string): T {
  const repaired = text.replace(/(?<=[:[,]\s*)(?:-?Infinity|NaN)(?=\s*[,}\]])/g, 'null');
  return JSON.parse(repaired) as T;
}

async function errorFrom(res: Response): Promise<Error> {
  let detail = '';
  try {
    detail = (await res.text()).trim().slice(0, 300);
  } catch {
    /* body unavailable */
  }
  return new Error(`API error: ${res.status} ${res.statusText}${detail ? ` — ${detail}` : ''}`);
}

async function getJson<T>(url: string, signal?: AbortSignal): Promise<T> {
  const res = await fetch(url, { signal });
  if (!res.ok) throw await errorFrom(res);
  return (await res.json()) as T;
}

async function getJsonLenient<T>(url: string, signal?: AbortSignal): Promise<T> {
  const res = await fetch(url, { signal });
  if (!res.ok) throw await errorFrom(res);
  return parseLenientJson<T>(await res.text());
}

// ---------- fetchers ----------

export async function fetchGeneRanking(gene: string, signal?: AbortSignal): Promise<JointApiResponse> {
  const symbol = gene.trim();
  if (!symbol) throw new Error('No gene symbol provided.');
  if (USE_MOCK) {
    await mockDelay();
    return mockGeneResponse(symbol);
  }
  const data = await getJsonLenient<JointApiResponse>(
    `${GENE_URL}?gene=${encodeURIComponent(symbol)}&top_n=${TOP_N}`,
    signal
  );
  if (!data || !Array.isArray(data.lines)) throw new Error('Unexpected /gene response shape.');
  return data;
}

export async function fetchMultiRanking(genes: string[], signal?: AbortSignal): Promise<JointApiResponse> {
  const clean = genes.map((g) => g.trim()).filter(Boolean);
  if (clean.length < 2) throw new Error('Provide at least two target genes for a joint ranking.');
  if (USE_MOCK) {
    await mockDelay();
    return mockGenesResponse(clean);
  }
  const query = clean.map(encodeURIComponent).join(',');
  const data = await getJsonLenient<JointApiResponse>(`${GENES_URL}?genes=${query}&top_n=${TOP_N}`, signal);
  if (!data || !Array.isArray(data.lines)) throw new Error('Unexpected /genes response shape.');
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
  if (USE_MOCK) {
    await mockDelay();
    return mockExcludeResponse(a, b);
  }
  const data = await getJson<ExcludeApiResponse>(
    `${EXCLUDE_URL}?gene_a=${encodeURIComponent(a)}&gene_b=${encodeURIComponent(b)}&top_n=${TOP_N}`,
    signal
  );
  if (!data || !Array.isArray(data.lines)) throw new Error('Unexpected /exclude response shape.');
  return data;
}

// ---------- helpers ----------

/** Spread near-identical scores across 5..100% for a readable bar. */
function relativize(values: number[], value: number): number {
  const min = Math.min(...values);
  const max = Math.max(...values);
  const span = max - min;
  return span > 0 ? Math.round((5 + ((value - min) / span) * 95) * 10) / 10 : 100;
}

/** Normalize a raw lineage string ("central_nervous_system" -> "Central nervous system"). */
function lineageOf(raw?: string | null): string {
  const v = (raw ?? '').trim().replace(/_/g, ' ');
  return v ? v.charAt(0).toUpperCase() + v.slice(1) : 'Unknown';
}

const displayName = (name: string | null, modelId: string) =>
  name && name !== 'None' ? name.toUpperCase() : modelId;

// ---------- adapters ----------

export function toSingleRanked(resp: JointApiResponse): SingleRanked[] {
  const joints = resp.lines.map((l) => l.joint_score);
  return resp.lines.map((l, i) => ({
    mode: 'single',
    rank: i + 1,
    modelId: l.model_id,
    cellLine: displayName(l.name, l.model_id),
    lineage: lineageOf(l.metadata?.lineage),
    score: l.joint_score,
    relativeScore: relativize(joints, l.joint_score),
  }));
}

export function toMultiRanked(resp: JointApiResponse): MultiRanked[] {
  const joints = resp.lines.map((l) => l.joint_score);
  return resp.lines.map((l, i) => {
    const geneScores: Record<string, number | null> = {};
    for (const g of resp.genes) {
      const v = l.scores?.[g];
      geneScores[g] = typeof v === 'number' && Number.isFinite(v) ? v : null;
    }
    return {
      mode: 'multi',
      rank: i + 1,
      modelId: l.model_id,
      cellLine: displayName(l.name, l.model_id),
      lineage: lineageOf(l.metadata?.lineage),
      jointScore: l.joint_score,
      genes: resp.genes,
      geneScores,
      limitingGene: l.limiting_gene ?? null,
      relativeScore: relativize(joints, l.joint_score),
    };
  });
}

export function toSelectivityRanked(resp: ExcludeApiResponse): SelectivityRanked[] {
  const sels = resp.lines.map((l) => l.selectivity);
  return resp.lines.map((l, i) => ({
    mode: 'selectivity',
    rank: i + 1,
    modelId: l.model_id,
    cellLine: displayName(l.name, l.model_id),
    lineage: lineageOf(l.metadata?.lineage),
    selectivity: l.selectivity,
    scoreHigh: l.score_a,
    scoreLow: l.score_b,
    geneHigh: resp.gene_a,
    geneLow: resp.gene_b,
    relativeScore: relativize(sels, l.selectivity),
  }));
}

export function metaFromGene(resp: JointApiResponse): ResultMeta {
  return {
    mode: 'single',
    primaryGene: resp.genes[0] ?? '',
    genes: resp.genes,
    total: resp.total_passing,
    showing: resp.lines.length,
  };
}

export function metaFromGenes(resp: JointApiResponse): ResultMeta {
  return {
    mode: 'multi',
    primaryGene: resp.genes.join(' + '),
    genes: resp.genes,
    floor: resp.floor,
    total: resp.total_passing,
    showing: resp.lines.length,
  };
}

export function metaFromExclude(resp: ExcludeApiResponse): ResultMeta {
  return {
    mode: 'selectivity',
    primaryGene: resp.gene_a,
    excludedGene: resp.gene_b,
    formula: `score_${resp.gene_a} × (1 − score_${resp.gene_b})`,
    total: resp.total_ranked,
    showing: resp.lines.length,
  };
}

/**
 * High-level entry point. Routes by what's selected:
 *   2+ targets            -> /prod/genes   (joint multi-gene; exclusions ignored)
 *   1 target + exclusion  -> /prod/exclude (selectivity)
 *   1 target              -> /prod/gene    (single, joint shape with one gene)
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
  if (USE_MOCK) {
    await mockDelay();
    return mockCellLineDetailResponse(g, id);
  }
  const data = await getJson<CellLineDetailApiResponse>(
    `${DETAIL_URL}?gene=${encodeURIComponent(g)}&model_id=${encodeURIComponent(id)}`,
    signal
  );
  if (!data || !data.cell_line) throw new Error('Unexpected cell-line detail response shape.');
  return data;
}

// Human-readable labels for the metadata block, in display order.
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
    items.push({ label, value: String(raw).replace(/_/g, ' ') });
  }
  return items;
}

export function toCellLineDetail(resp: CellLineDetailApiResponse): CellLineDetail {
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

  const name = displayName(resp.cell_line.name, resp.cell_line.model_id);

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
      name: displayName(a.name, a.model_id),
      similarity: a.similarity,
    })),
  };
}
