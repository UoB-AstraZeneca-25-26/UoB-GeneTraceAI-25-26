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

export type ViewType = 'query' | 'results' | 'profile' | 'about' | 'assistant' | 'guide' | 'reference';

// ---- Live API result model ----

export type ScoreTier = 'HIGH' | 'MEDIUM' | 'LOW';

export type ResultMode = 'single' | 'selectivity' | 'multi';

export type BaseRanked = {
  rank: number;
  modelId: string;   // DepMap ACH id (primary identifier)
  cellLine: string;  // display name / alias (secondary)
  lineage: string;   // tissue lineage, e.g. "skin"
  relativeScore: number;  // min-max normalized within the set, for bar width
  lineageScore: number;   // 0-1 score relative to same-lineage lines (within-lineage)
};

/** Result of /prod/gene (joint shape with a single gene). tier / n_layers /
 *  driver are no longer in the list response — they live in the detail endpoint. */
export type SingleRanked = BaseRanked & {
  mode: 'single';
  score: number; // joint_score for the single gene
};

/** Result of GET /prod/exclude/many?gene_a=<HIGH>&exclude=<LOW1>&exclude=<LOW2>... */
export type SelectivityRanked = BaseRanked & {
  mode: 'selectivity';
  selectivity: number; // score_high * Π(1 - score_low_i)
  geneHigh: string;
  scoreHigh: number;
  excludedGenes: string[];
  exclusionScores: Record<string, number>;
};

/** Result of GET /genes?query=<A,B,...> (joint multi-gene ranking) */
export type MultiRanked = BaseRanked & {
  mode: 'multi';
  jointScore: number;
  genes: string[];                           // queried genes, in order
  geneScores: Record<string, number | null>; // per-gene; null = NaN/no evidence
  limitingGene: string | null;               // gene that caps the joint score
};

export type RankedCellLine = SingleRanked | SelectivityRanked | MultiRanked;

export type ResultMeta = {
  mode: ResultMode;
  primaryGene: string;   // gene_high (selectivity), gene (single), or joined list (multi)
  ensg?: string;
  excludedGenes?: string[]; // gene_b list (selectivity only)
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
export type InspectTarget = {
  gene: string;
  modelId: string;
  cellLine: string;
  ranked?: RankedCellLine; // the ranked row this line came from — carries how it scored in this query
};

export type MetadataItem = { label: string; value: string };
export type SimilarLine = { modelId: string; name: string; similarity: number; lineage?: string; primaryDisease?: string };

/** A raw per-source measurement behind a modality's level, e.g. DepMap TPM.
 *  Units differ per source, so `max` gives a per-source scale for its bar. */
export type SourceMeasurement = { source: string; value: number; unit: string; max: number };

export type CellLineDetail = {
  gene: string;
  ensg: string;
  modelId: string;
  cellLine: string;
  rank: number;
  total: number;
  score: number;
  lineageScore: number | null;  // 0-1 within-lineage score (null if backend omits)
  lineageRank: number | null;   // rank within the same lineage
  lineageTotal: number | null;  // number of lines in the same lineage
  tier: string;
  nLayers: number;
  driverAlteration: boolean;
  pMutation: number | null;
  pFusion: number | null;
  hasCnaAlteration: boolean;
  expressionLevel: number | null;   // within-gene percentile (0..1), or null if unmeasured
  proteomicsLevel: number | null;   // within-gene percentile (0..1), or null if unmeasured
  expressionSources: string[];      // which expression datasets contributed (DepMap/HPA/GEO)
  proteomicsSources: string[];      // which proteomics datasets contributed (ProCan/CCLE)
  expressionMeasurements: SourceMeasurement[]; // raw per-source expression values, if reported
  proteomicsMeasurements: SourceMeasurement[]; // raw per-source proteomics values, if reported
  tracks: TrackScores;         // per-layer signals reported by this endpoint
  metadata: MetadataItem[];    // ordered, cleaned for display
  alternatives: SimilarLine[]; // RNA-similar models
};

// ---- Assistant / chat ----

export type ChatRole = 'user' | 'assistant';
export type ChatMessage = { role: ChatRole; content: string };

/** Snapshot of the current query + top results, handed to the agent so it can
 *  explain the ranking with real numbers. */
export type AssistantContext = {
  mode: ResultMode | null;
  genes: string[];
  exclusions: string[];
  lineage: string;
  primaryGene: string | null;
  total: number | null;
  formula: string | null;
  topLines: { rank: number; cellLine: string; modelId: string; score: number; tier: string | null }[];
};

// ---- Agent-backed methodology / gene alias enrichment ----

export type MethodologyStepStatus = 'active' | 'skipped' | 'inverted';

export type MethodologyStep = {
  key: string;
  formula?: string;
  value?: string;
  status?: MethodologyStepStatus;
};

export type GeneAliasInfo = {
  found?: boolean;
  full_name?: string;
  synonyms?: string[];
};
