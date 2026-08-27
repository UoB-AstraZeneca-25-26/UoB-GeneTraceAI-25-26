// Offline mock data. Returns the same response shapes as the live API so every
// adapter and view runs identically. Toggle via USE_MOCK in api.ts.
import type {
  JointApiResponse,
  ExcludeManyApiResponse,
  ExcludeManyApiLine,
  CellLineDetailApiResponse,
  LineMetadata,
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

function mk(
  model_id: string, name: string, lineage: string, lineage_subtype: string,
  primary_disease: string, subtype: string, sex: string, age: string,
  origin: string, site: string, ncit: string
): PoolLine {
  return { model_id, name, lineage, lineage_subtype, primary_disease, subtype, sex, age, origin, site, ncit };
}

const POOL: PoolLine[] = [
  mk('ACH-000219', 'a-375', 'skin', 'melanoma', 'skin cancer', 'melanoma', 'female', '54', 'primary', 'skin', 'amelanotic melanoma'),
  mk('ACH-000322', 'ht-144', 'skin', 'melanoma', 'skin cancer', 'melanoma', 'male', '29', 'metastasis', 'soft tissue', 'melanoma'),
  mk('ACH-000014', 'hs 294t', 'skin', 'melanoma', 'skin cancer', 'melanoma', 'male', '56', 'metastasis', 'lymph node', 'melanoma'),
  mk('ACH-000008', 'a101d', 'skin', 'melanoma', 'skin cancer', 'melanoma', 'female', '43', 'primary', 'skin', 'melanoma'),
  mk('ACH-000707', 'sk-mel-28', 'skin', 'melanoma', 'skin cancer', 'melanoma', 'male', '51', 'primary', 'skin', 'malignant melanoma'),
  mk('ACH-000620', 'jhh-1', 'liver', 'hepatocellular_carcinoma', 'liver cancer', 'HCC', 'male', '57', 'primary', 'liver', 'hepatocellular carcinoma'),
  mk('ACH-000758', 'mkn74', 'gastric', 'adenocarcinoma', 'gastric cancer', 'adenocarcinoma', 'male', '37', 'metastasis', 'liver', 'gastric adenocarcinoma'),
  mk('ACH-000552', 'a549', 'lung', 'NSCLC', 'lung cancer', 'adenocarcinoma', 'male', '58', 'primary', 'lung', 'lung adenocarcinoma'),
  mk('ACH-000681', 'calu-3', 'lung', 'NSCLC', 'lung cancer', 'adenocarcinoma', 'male', '25', 'metastasis', 'pleural effusion', 'lung adenocarcinoma'),
  mk('ACH-000788', 'nci-h1975', 'lung', 'NSCLC', 'lung cancer', 'adenocarcinoma', 'female', '68', 'primary', 'lung', 'lung adenocarcinoma'),
  mk('ACH-000585', 'hct116', 'colorectal', 'adenocarcinoma', 'colorectal cancer', 'adenocarcinoma', 'male', '48', 'primary', 'colon', 'colorectal carcinoma'),
  mk('ACH-000971', 'sw480', 'colorectal', 'adenocarcinoma', 'colorectal cancer', 'adenocarcinoma', 'male', '50', 'primary', 'colon', 'colon adenocarcinoma'),
  mk('ACH-000148', 'ht-29', 'colorectal', 'adenocarcinoma', 'colorectal cancer', 'adenocarcinoma', 'female', '44', 'primary', 'colon', 'colorectal adenocarcinoma'),
  mk('ACH-001307', '8505c', 'thyroid', 'anaplastic', 'thyroid cancer', 'anaplastic', 'female', '78', 'primary', 'thyroid', 'anaplastic thyroid carcinoma'),
  mk('ACH-000883', 'sw 1783', 'central_nervous_system', 'glioma', 'brain cancer', 'astrocytoma', 'female', '32', 'primary', 'brain', 'astrocytoma'),
  mk('ACH-000410', 'u-87 mg', 'central_nervous_system', 'glioblastoma', 'brain cancer', 'GBM', 'male', '69', 'primary', 'brain', 'glioblastoma'),
  mk('ACH-000659', 'sclc-21h', 'lung', 'SCLC', 'lung cancer', 'small cell', 'male', '75', 'primary', 'lung', 'small cell lung carcinoma'),
  mk('ACH-000976', 'hucct1', 'biliary_tract', 'cholangiocarcinoma', 'bile duct cancer', 'cholangiocarcinoma', 'male', '56', 'primary', 'bile duct', 'cholangiocarcinoma'),
  mk('ACH-000336', 'oci-aml3', 'blood', 'AML', 'leukemia', 'acute myeloid', 'male', '57', 'primary', 'bone marrow', 'acute myeloid leukemia'),
  mk('ACH-000956', 'k562', 'blood', 'CML', 'leukemia', 'chronic myeloid', 'female', '53', 'primary', 'bone marrow', 'chronic myelogenous leukemia'),
  mk('ACH-000113', 'jurkat', 'blood', 'T-ALL', 'leukemia', 'acute lymphoblastic', 'male', '14', 'primary', 'peripheral blood', 'T-cell acute lymphoblastic leukemia'),
  mk('ACH-000601', 'panc-1', 'pancreas', 'ductal', 'pancreatic cancer', 'ductal adenocarcinoma', 'male', '56', 'primary', 'pancreas', 'pancreatic ductal adenocarcinoma'),
  mk('ACH-000019', 'mcf7', 'breast', 'ductal', 'breast cancer', 'invasive ductal', 'female', '69', 'metastasis', 'pleural effusion', 'breast carcinoma'),
  mk('ACH-000212', 'mda-mb-231', 'breast', 'triple_negative', 'breast cancer', 'TNBC', 'female', '51', 'metastasis', 'pleural effusion', 'breast adenocarcinoma'),
  mk('ACH-000001', 'nih:ovcar-3', 'ovary', 'serous', 'ovarian cancer', 'high-grade serous', 'female', '60', 'metastasis', 'ascites', 'ovarian serous adenocarcinoma'),
  mk('ACH-000031', 'du145', 'prostate', 'adenocarcinoma', 'prostate cancer', 'adenocarcinoma', 'male', '69', 'metastasis', 'brain', 'prostate carcinoma'),
  mk('ACH-000445', 'caki-1', 'kidney', 'clear_cell', 'renal cancer', 'ccRCC', 'male', '49', 'metastasis', 'skin', 'renal cell carcinoma'),
  mk('ACH-000234', 'detroit 562', 'upper_aerodigestive', 'pharynx', 'head and neck cancer', 'squamous cell', 'female', '50', 'metastasis', 'pleural effusion', 'pharyngeal carcinoma'),
  mk('ACH-001421', 'weri-rb-1', 'eye', 'retinoblastoma', 'eye cancer', 'retinoblastoma', 'female', '2', 'primary', 'eye', 'retinoblastoma'),
];

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
const metaOf = (p: PoolLine): LineMetadata => ({
  lineage: p.lineage,
  lineage_subtype: p.lineage_subtype,
  primary_disease: p.primary_disease,
  sex: p.sex,
});

// ---------- joint shape (/gene single, /genes multi) ----------
function mockJointResponse(genesRaw: string[]): JointApiResponse {
  const genes = genesRaw.map(clampGene).filter(Boolean);
  const single = genes.length === 1;

  const lines = POOL.map((p) => {
    const scores: Record<string, number> = {};
    let present = 0;
    let minPresent = 1;
    let limiter: string | null = null;
    genes.forEach((g) => {
      // single-gene queries always score; multi drops ~35% -> partial coverage
      const has = single || rand('m' + g + p.model_id) > 0.35;
      if (has) {
        const s = round(Math.max(0.5, 0.9 + rand('v' + g + p.model_id) * 0.0999));
        scores[g] = s;
        if (s < minPresent) {
          minPresent = s;
          limiter = g;
        }
        present++;
      } else {
        scores[g] = NaN;
      }
    });
    return {
      model_id: p.model_id,
      name: p.name,
      joint_score: present > 0 ? round(minPresent) : 0,
      limiting_gene: limiter ?? genes[0],
      scores,
      metadata: metaOf(p),
      present,
    };
  })
    .filter((l) => l.present > 0)
    .sort((x, y) => y.joint_score - x.joint_score)
    .map(({ present: _p, ...line }) => line);

  return { genes, lineage: [], floor: single ? 0 : 0.5, total_passing: single ? 1039 : 248, lines };
}

export const mockGeneResponse = (gene: string): JointApiResponse => mockJointResponse([gene]);
export const mockGenesResponse = (genes: string[]): JointApiResponse => mockJointResponse(genes);

// ---------- /exclude/many ----------
export function mockExcludeManyResponse(aRaw: string, bsRaw: string[]): ExcludeManyApiResponse {
  const a = clampGene(aRaw);
  const bs = bsRaw.map(clampGene).filter(Boolean);
  const lines: ExcludeManyApiLine[] = POOL.map((p) => {
    const sa = round(0.85 + rand(a + p.model_id) * 0.149);
    const exclusion_scores: Record<string, number> = {};
    let selectivity = sa;
    bs.forEach((b) => {
      const sb = round(rand(b + p.model_id) * 0.5);
      exclusion_scores[b] = sb;
      selectivity *= 1 - sb;
    });
    return {
      model_id: p.model_id,
      name: p.name,
      score_a: sa,
      exclusion_scores,
      selectivity: round(selectivity),
      metadata: metaOf(p),
    };
  }).sort((x, y) => y.selectivity - x.selectivity);

  return { gene_a: a, excluded_genes: bs, lineage: [], total_ranked: 912, lines };
}

// ---------- /gene/detail ----------
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

  const exprLevel = round(0.3 + rand('rna' + gene + p.model_id) * 0.69, 3);
  const protLevel = rand('prot' + gene + p.model_id) > 0.3 ? round(0.2 + rand('prot2' + gene + p.model_id) * 0.79, 3) : null;

  // Which datasets contributed — DepMap always present for expression.
  const exprSources = [
    'DepMap',
    ...(rand('hpa' + gene + p.model_id) > 0.4 ? ['HPA'] : []),
    ...(rand('geo' + gene + p.model_id) > 0.6 ? ['GEO'] : []),
  ];
  const protSources =
    protLevel != null
      ? [
          ...(rand('procan' + gene + p.model_id) > 0.35 ? ['ProCan'] : []),
          ...(rand('ccle' + gene + p.model_id) > 0.45 ? ['CCLE'] : []),
        ]
      : [];

  // Raw per-source measurements. Units differ per source; `max` sets each bar's
  // own scale. Values track the modality level so the numbers stay coherent.
  const EXPR_UNIT: Record<string, { unit: string; max: number }> = {
    DepMap: { unit: 'log2(TPM+1)', max: 16 },
    HPA: { unit: 'nTPM', max: 140 },
    GEO: { unit: 'z-score', max: 4 },
  };
  const PROT_UNIT: Record<string, { unit: string; max: number }> = {
    ProCan: { unit: 'log2 intensity', max: 10 },
    CCLE: { unit: 'norm. quant', max: 6 },
  };
  const exprMeasurements = exprSources.map((s) => {
    const { unit, max } = EXPR_UNIT[s];
    const v = exprLevel * max * (0.7 + rand('em' + s + gene + p.model_id) * 0.3);
    return { source: s, value: round(Math.min(v, max), 2), unit, max };
  });
  const protMeasurements = protSources.map((s) => {
    const { unit, max } = PROT_UNIT[s];
    const v = (protLevel ?? 0) * max * (0.7 + rand('pm' + s + gene + p.model_id) * 0.3);
    return { source: s, value: round(Math.min(v, max), 2), unit, max };
  });

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
    expression_level: exprLevel,
    proteomics_level: protLevel,
    expression_sources: exprSources,
    proteomics_sources: protSources,
    expression_measurements: exprMeasurements,
    proteomics_measurements: protMeasurements,
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
