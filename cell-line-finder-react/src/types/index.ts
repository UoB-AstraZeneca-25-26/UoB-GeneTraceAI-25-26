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
  tracks?: TrackScores;    // per-layer detail; populated by a future /cellline endpoint
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

// ---- Evidence tracks (adopted from GeneTraceAI's layer model) ----

export type TrackKey = 'expression' | 'proteomics' | 'cna' | 'mutation' | 'fusion' | 'signature';

/** What the track's number actually is — so views never imply a magnitude the
 *  pipeline doesn't have. level/ratio are real gradients; call is categorical. */
export type TrackKind = 'level' | 'ratio' | 'call';

export type EvidenceTrack = {
  key: TrackKey;
  label: string;
  color: string;
  kind: TrackKind;
};

export type TrackDatum = { present: boolean; score: number; finding?: string };
export type TrackScores = Partial<Record<TrackKey, TrackDatum>>;

// ---- Tier presentation (separated from score, never derived from it) ----

export type TierTone = 'high' | 'medium' | 'low' | 'unknown';
export type TierPresentation = { label: string; tone: TierTone };

// ---- Verdict / abstention states ----

export type VerdictState = 'RANKED' | 'LOW_SEPARATION' | 'NO_EVIDENCE';
export type Verdict = {
  state: VerdictState;
  headline: string;
  detail: string[];
  qualifyRanking?: boolean;
};