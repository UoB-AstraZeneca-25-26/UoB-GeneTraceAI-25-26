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

export type ScoreTier = 'HIGH' | 'MEDIUM' | 'LOW' | 'CONTEXT' | 'NO_EVIDENCE';

export type ResultMode = 'single' | 'selectivity' | 'multi' | 'jointSelectivity';

export type BaseRanked = {
  rank: number;
  modelId: string;   // DepMap ACH id (primary identifier)
  cellLine: string;  // display name / alias (secondary)
  lineage: string;   // tissue lineage, e.g. "skin"
  relativeScore: number;  // min-max normalized within the set, for bar width
  lineageScore: number;   // 0-1 score relative to same-lineage lines (within-lineage)
  // Cell-line context, for deciding whether a top-scoring line is actually a
  // sensible experimental model -- not every line has every field.
  primaryDisease?: string;
  subtype?: string;
  growthPattern?: string;  // e.g. "2d: adherent" -- practical culture-protocol info
  sex?: string;
  rrid?: string;           // Cellosaurus RRID, for ordering/citing the line
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

/** Result of GET /prod/genes/exclude/many?genes=<A,B,...>&exclude=<LOW1>&exclude=<LOW2>...
 *  Combines rank()'s N-target weakest-link joint score with exclude_many()'s
 *  N-exclusion selectivity penalty: combinedScore = jointScore * Pi(1 - exclusionScore_i). */
export type JointSelectivityRanked = BaseRanked & {
  mode: 'jointSelectivity';
  jointScore: number;                        // min across target genes (weakest-link)
  combinedScore: number;                     // jointScore * Pi(1 - exclusion scores) -- the ranking metric
  genes: string[];                           // queried target genes, in order
  geneScores: Record<string, number | null>; // per-target-gene; null = NaN/no evidence
  limitingGene: string | null;               // target gene that caps the joint score
  excludedGenes: string[];
  exclusionScores: Record<string, number>;
};

export type RankedCellLine = SingleRanked | SelectivityRanked | MultiRanked | JointSelectivityRanked;

export type ResultMeta = {
  mode: ResultMode;
  primaryGene: string;   // gene_high (selectivity), gene (single), or joined list (multi / jointSelectivity)
  ensg?: string;
  excludedGenes?: string[]; // gene_b list (selectivity), exclude list (jointSelectivity)
  formula?: string;      // selectivity, jointSelectivity
  genes?: string[];      // multi, jointSelectivity
  floor?: number;        // multi, jointSelectivity
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
  viaAlternative?: boolean; // true when reached by clicking an RNA-similar alternative, not from the ranked results table
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
  mutationDriver: boolean; // p_mutation >= 0.5 -- this layer's own share of driverAlteration
  fusionDriver: boolean;   // p_fusion >= 0.5 -- this layer's own share of driverAlteration
  hasCnaAlteration: boolean;
  expressionLevel: number | null;   // within-gene percentile (0..1), or null if unmeasured
  proteomicsLevel: number | null;   // within-gene percentile (0..1), or null if unmeasured
  expressionBySource: Record<string, number | null>; // per-source level (0..1), e.g. {DepMap: 0.8, HPA: null}
  proteomicsBySource: Record<string, number | null>; // per-source level (0..1), e.g. {ProCan: 0.4, CCLE: 0.6}
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

// ---- Lineage-grouped ranking (GET /prod/gene/by-lineage?gene=<GENE>) ----
// Top 10 lineages by their best candidate's score, top 5 lines per lineage.
// Not affected by top_n -- sizing is fixed by this rule.

export type LineageGroupedLine = {
  modelId: string;
  cellLine: string;
  lineage: string;
  coreScore: number;
  rankWithinLineage: number; // 1..5 within this lineage (ties broken by modelId)
  rankGlobal: number;        // pooled rank across every scoreable line for this gene (ties share a rank)
  primaryDisease?: string;
  subtype?: string;
  growthPattern?: string;
  sex?: string;
  rrid?: string;
};

export type LineageGroupedResult = {
  gene: string;
  lineagesReturned: number;      // <= 10
  totalScoreable: number;        // all scoreable lines for this gene, not just the ones shown
  lineageUnassignedCount: number; // scoreable lines excluded because their cell line has no lineage on record
  lines: LineageGroupedLine[];   // <= 50
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
