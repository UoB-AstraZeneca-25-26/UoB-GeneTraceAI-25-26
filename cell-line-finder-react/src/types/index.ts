export type CellLineData = {
  cellLine: string;
  tissue: string;
  gene: string;
  tpmExpression: number;
  proteinExpression: number;
  hasSomaticVariant: boolean;
};

export type QueryParams = {
  targets: string[];
  exclusions: string[];
  lineage: string;
  topK: number;
};

export type ScoredResult = {
  cellLine: string;
  tissue: string;
  meanTpm: number;
  mutations: number;
  confidenceScore: number;
  exclusionPenalty: number;
};

export type ViewType = 'query' | 'results' | 'profile' | 'about';

// ---- Live API result model ----

export type ScoreTier = 'HIGH' | 'MEDIUM' | 'LOW';

export type ResultMode = 'single' | 'selectivity';

type BaseRanked = {
  rank: number;
  modelId: string;   // DepMap ACH id
  cellLine: string;  // display name
  relativeScore: number; // min-max normalized within the set, for bar width
};

/** Result of GET /gene?query=<SYMBOL> */
export type SingleRanked = BaseRanked & {
  mode: 'single';
  score: number;           // raw 0..1
  confidenceScore: number; // score * 100
  tier: ScoreTier;
  nLayers: number;
  driverAlteration: boolean;
};

/** Result of GET /prod/exclude?gene_a=<HIGH>&gene_b=<LOW> */
export type SelectivityRanked = BaseRanked & {
  mode: 'selectivity';
  selectivity: number; // score_high * (1 - score_low)
  geneHigh: string;
  geneLow: string;
  scoreHigh: number;
  scoreLow: number;
};

export type RankedCellLine = SingleRanked | SelectivityRanked;

export type ResultMeta = {
  mode: ResultMode;
  primaryGene: string;   // gene_high (selectivity) or gene (single)
  ensg?: string;
  excludedGene?: string; // gene_low (selectivity only)
  formula?: string;      // selectivity only
  total: number;
  showing: number;
};