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

export type ResultMode = 'single' | 'selectivity' | 'multi';

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

/** Result of GET /genes?query=<A,B,...> (joint multi-gene ranking) */
export type MultiRanked = BaseRanked & {
  mode: 'multi';
  jointScore: number;
  genes: string[];                           // queried genes, in order
  geneScores: Record<string, number | null>; // per-gene; null = NaN/no evidence
};

export type RankedCellLine = SingleRanked | SelectivityRanked | MultiRanked;

export type ResultMeta = {
  mode: ResultMode;
  primaryGene: string;   // gene_high (selectivity), gene (single), or joined list (multi)
  ensg?: string;
  excludedGene?: string; // gene_low (selectivity only)
  formula?: string;      // selectivity only
  genes?: string[];      // multi only
  floor?: number;        // multi only
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

export type TierTone = 'high' | 'medium' | 'low' | 'context' | 'unknown';
export type TierPresentation = { label: string; tone: TierTone };

// ---- Verdict / abstention states ----

export type VerdictState = 'RANKED' | 'LOW_SEPARATION' | 'NO_EVIDENCE' | 'PARTIAL_COVERAGE';
export type Verdict = {
  state: VerdictState;
  headline: string;
  detail: string[];
  qualifyRanking?: boolean;
};

// ---- Cell-line detail (GET /gene?query=<GENE>&cell_line=<MODEL_ID>) ----

/** What the Inspect action carries: which gene's detail to fetch, for which model. */
export type InspectTarget = { gene: string; modelId: string; cellLine: string };

export type MetadataItem = { label: string; value: string };
export type SimilarLine = { modelId: string; name: string; similarity: number };

export type CellLineDetail = {
  gene: string;
  ensg: string;
  modelId: string;
  cellLine: string;
  rank: number;
  total: number;
  score: number;
  tier: string;
  nLayers: number;
  driverAlteration: boolean;
  pMutation: number | null;
  pFusion: number | null;
  hasCnaAlteration: boolean;
  tracks: TrackScores;         // per-layer signals reported by this endpoint
  metadata: MetadataItem[];    // ordered, cleaned for display
  alternatives: SimilarLine[]; // RNA-similar models
};