// Offline mock data for UI work while the live endpoints are down.
// Returns the exact same response shapes as the real API, so every adapter,
// verdict, tier and view runs identically. Toggle via USE_MOCK in api.ts.
import type {
  GeneApiResponse,
  ExcludeApiResponse,
  ExcludeApiLine,
  GenesApiResponse,
  CellLineDetailApiResponse,
} from './api';

interface PoolLine {
  model_id: string;
  name: string;
  lineage: string;
  lineage_subtype: string;
  primary_disease: string;
  subtype: string;
  sex: string;
  age: string;
  origin: string;
  site: string;
  ncit: string;
}

const POOL: PoolLine[] = [
  mk('ACH-000219', 'a-375', 'skin', 'melanoma', 'skin cancer', 'melanoma', 'female', '54', 'primary', 'skin', 'amelanotic melanoma'),
  mk('ACH-000322', 'ht-144', 'skin', 'melanoma', 'skin cancer', 'melanoma', 'male', '29', 'metastasis', 'soft tissue', 'melanoma'),
  mk('ACH-000014', 'hs 294t', 'skin', 'melanoma', 'skin cancer', 'melanoma', 'male', '56', 'metastasis', 'lymph node', 'melanoma'),
  mk('ACH-000008', 'a101d', 'skin', 'melanoma', 'skin cancer', 'melanoma', 'female', '43', 'primary', 'skin', 'melanoma'),
  mk('ACH-000620', 'jhh-1', 'liver', 'hepatocellular', 'liver cancer', 'HCC', 'male', '57', 'primary', 'liver', 'hepatocellular carcinoma'),
  mk('ACH-000758', 'mkn74', 'gastric', 'adenocarcinoma', 'gastric cancer', 'adenocarcinoma', 'male', '37', 'metastasis', 'liver', 'gastric adenocarcinoma'),
  mk('ACH-000552', 'a549', 'lung', 'NSCLC', 'lung cancer', 'adenocarcinoma', 'male', '58', 'primary', 'lung', 'lung adenocarcinoma'),
  mk('ACH-000681', 'calu-3', 'lung', 'NSCLC', 'lung cancer', 'adenocarcinoma', 'male', '25', 'metastasis', 'pleural effusion', 'lung adenocarcinoma'),
  mk('ACH-000585', 'hct116', 'colorectal', 'adenocarcinoma', 'colorectal cancer', 'adenocarcinoma', 'male', '48', 'primary', 'colon', 'colorectal carcinoma'),
  mk('ACH-000971', 'sw480', 'colorectal', 'adenocarcinoma', 'colorectal cancer', 'adenocarcinoma', 'male', '50', 'primary', 'colon', 'colon adenocarcinoma'),
  mk('ACH-001307', '8505c', 'thyroid', 'anaplastic', 'thyroid cancer', 'anaplastic', 'female', '78', 'primary', 'thyroid', 'anaplastic thyroid carcinoma'),
  mk('ACH-000883', 'sw 1783', 'CNS/brain', 'glioma', 'brain cancer', 'astrocytoma', 'female', '32', 'primary', 'brain', 'astrocytoma'),
  mk('ACH-000659', 'sclc-21h', 'lung', 'SCLC', 'lung cancer', 'small cell', 'male', '75', 'primary', 'lung', 'small cell lung carcinoma'),
  mk('ACH-000976', 'hucct1', 'biliary', 'cholangiocarcinoma', 'bile duct cancer', 'cholangiocarcinoma', 'male', '56', 'primary', 'bile duct', 'cholangiocarcinoma'),
  mk('ACH-000336', 'oci-aml3', 'blood', 'AML', 'leukemia', 'acute myeloid', 'male', '57', 'primary', 'bone marrow', 'acute myeloid leukemia'),
  mk('ACH-000601', 'panc-1', 'pancreas', 'ductal', 'pancreatic cancer', 'ductal adenocarcinoma', 'male', '56', 'primary', 'pancreas', 'pancreatic ductal adenocarcinoma'),
];

function mk(
  model_id: string, name: string, lineage: string, lineage_subtype: string,
  primary_disease: string, subtype: string, sex: string, age: string,
  origin: string, site: string, ncit: string
): PoolLine {
  return { model_id, name, lineage, lineage_subtype, primary_disease, subtype, sex, age, origin, site, ncit };
}

// Deterministic 0..1 from a string, so a given query is stable across renders.
function rand(seed: string): number {
  let h = 2166136261;
  for (let i = 0; i < seed.length; i++) {
    h ^= seed.charCodeAt(i);
    h = Math.imul(h, 16777619);
  }
  return ((h >>> 0) % 100000) / 100000;
}

const round = (x: number, d = 4) => Math.round(x * 10 ** d) / 10 ** d;
const clampGene = (g: string) => (g || '').trim().toUpperCase();

function tierFor(rankIndex: number): string {
  if (rankIndex < 3) return 'HIGH';
  if (rankIndex < 8) return 'MEDIUM';
  return 'LOW';
}

// ---------- /gene ----------
export function mockGeneResponse(geneRaw: string): GeneApiResponse {
  const gene = clampGene(geneRaw);
  const ranked = [...POOL]
    .map((p) => ({ p, r: rand(gene + p.model_id) }))
    .sort((a, b) => b.r - a.r);

  const lines = ranked.map(({ p }, i) => ({
    rank: i + 1,
    model_id: p.model_id,
    name: p.name,
    // clustered near the ceiling on purpose -> exercises the low-separation verdict
    score: round(Math.max(0.5, 0.9998 - i * 0.0009 - rand('s' + gene + p.model_id) * 0.0003)),
    tier: tierFor(i),
    n_layers: 1 + Math.floor(rand('l' + gene + p.model_id) * 3),
    driver_alteration: rand('d' + gene + p.model_id) > 0.55,
  }));

  return { gene, ensg: 'ENSG' + (10000000 + Math.floor(rand(gene) * 8999999)), total: 1523, showing: lines.length, lines };
}

// ---------- /exclude ----------
export function mockExcludeResponse(aRaw: string, bRaw: string): ExcludeApiResponse {
  const a = clampGene(aRaw);
  const b = clampGene(bRaw);
  const lines = POOL.map((p) => {
    const hi = round(0.85 + rand(a + p.model_id) * 0.149);
    const lo = round(rand(b + p.model_id) * 0.5);
    const sel = round(hi * (1 - lo));
    const line: ExcludeApiLine = { model_id: p.model_id, name: p.name, selectivity: sel };
    line[`score_${a}`] = hi;
    line[`score_${b}`] = lo;
    return line;
  }).sort((x, y) => y.selectivity - x.selectivity);

  return { gene_high: a, gene_low: b, formula: `score_${a} × (1 − score_${b})`, total_ranked: 1493, lines };
}

// ---------- /genes (joint multi-gene) ----------
export function mockGenesResponse(genesRaw: string[]): GenesApiResponse {
  const genes = genesRaw.map(clampGene).filter(Boolean);
  const lines = POOL.map((p) => {
    const scores: Record<string, number> = {};
    let present = 0;
    let minPresent = 1;
    genes.forEach((g) => {
      // ~40% of gene/line combos are unscored -> NaN, so partial coverage dominates
      const has = rand('m' + g + p.model_id) > 0.4;
      if (has) {
        const s = round(Math.max(0.5, 0.999 - rand('v' + g + p.model_id) * 0.06));
        scores[g] = s;
        minPresent = Math.min(minPresent, s);
        present++;
      } else {
        scores[g] = NaN;
      }
    });
    return { model_id: p.model_id, name: p.name, joint_score: present > 0 ? round(minPresent) : 0, scores, present };
  })
    .filter((l) => l.present > 0)
    .sort((x, y) => y.joint_score - x.joint_score)
    .map(({ present: _p, ...line }) => line);

  return { genes, floor: 0.5, total_passing: 501, lines };
}

// ---------- /gene?...&cell_line=<id> ----------
export function mockCellLineDetailResponse(geneRaw: string, modelId: string): CellLineDetailApiResponse {
  const gene = clampGene(geneRaw);
  const p = POOL.find((x) => x.model_id === modelId) ?? POOL[0];
  const score = round(0.5 + rand(gene + p.model_id) * 0.49);
  const tier = score > 0.95 ? 'HIGH' : score > 0.85 ? 'MEDIUM' : score > 0.7 ? 'CONTEXT' : 'LOW';
  const driver = rand('d' + gene + p.model_id) > 0.5;

  const alts = POOL.filter((x) => x.model_id !== p.model_id)
    .map((x) => ({ x, sim: round(0.55 + rand(p.model_id + x.model_id) * 0.35, 3) }))
    .sort((m, n) => n.sim - m.sim)
    .slice(0, 5)
    .map(({ x, sim }) => ({ model_id: x.model_id, name: x.name, similarity: sim }));

  return {
    gene,
    ensg: 'ENSG' + (10000000 + Math.floor(rand(gene) * 8999999)),
    cell_line: { model_id: p.model_id, name: p.name },
    rank: 1 + Math.floor(rand('r' + gene + p.model_id) * 1522),
    total: 1523,
    score,
    tier,
    n_layers: 1 + Math.floor(rand('l' + gene + p.model_id) * 3),
    driver_alteration: driver,
    p_mutation: driver ? round(0.6 + rand('pm' + p.model_id) * 0.39, 3) : (rand('pm' + p.model_id) > 0.5 ? round(rand('pm2' + p.model_id) * 0.4, 3) : null),
    p_fusion: rand('pf' + p.model_id) > 0.75 ? round(rand('pf2' + p.model_id), 3) : null,
    has_cna_alteration: rand('cna' + gene + p.model_id) > 0.6,
    metadata: {
      lineage: p.lineage,
      lineage_subtype: p.lineage_subtype,
      primary_disease: p.primary_disease,
      subtype: p.subtype,
      sex: p.sex,
      age: p.age,
      default_growth_pattern: '2d: adherent',
      primary_or_metastasis: p.origin,
      sample_collection_site: p.site,
      cellosaurus_ncit_disease: p.ncit,
    },
    rna_alternatives: alts,
  };
}

export const mockDelay = (ms = 350) => new Promise((r) => setTimeout(r, ms));