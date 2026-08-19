import {
  EvidenceTrack,
  RankedCellLine,
  ScoreTier,
  TierPresentation,
  TierTone,
  TrackKey,
  Verdict,
} from '../types';

/* The six evidence layers every cell line is scored on. `kind` records what the
   number actually is, so the UI never implies a magnitude the pipeline lacks:
     level -> within-gene percentile of a measured quantity (real gradient)
     ratio -> deviation from diploid copy number (real gradient)
     call  -> categorical evidence, shown full/half strength, never interpolated */
export const TRACKS: EvidenceTrack[] = [
  { key: 'expression', label: 'Expression', color: '#4f46e5', kind: 'level' },
  { key: 'proteomics', label: 'Proteomics', color: '#059669', kind: 'level' },
  { key: 'cna', label: 'Copy number', color: '#0891b2', kind: 'ratio' },
  { key: 'mutation', label: 'Mutation', color: '#e11d48', kind: 'call' },
  { key: 'fusion', label: 'Fusion', color: '#9333ea', kind: 'call' },
  { key: 'signature', label: 'Signature', color: '#d97706', kind: 'call' },
];

export const TRACK_BY_KEY = Object.fromEntries(TRACKS.map((t) => [t.key, t])) as Record<
  TrackKey,
  EvidenceTrack
>;

export const pct = (x: number) => Math.round((x || 0) * 100);

/* The tier is an ordinal band from the pipeline. It is NOT derived from the score
   here — deriving a band from the number would re-introduce exactly the
   calibrated-confidence claim the scoring layer refuses to make. */
export function confTier(tier: ScoreTier | string): TierPresentation {
  switch (String(tier || '').toUpperCase()) {
    case 'HIGH':
      return { label: 'High confidence', tone: 'high' };
    case 'MEDIUM':
      return { label: 'Moderate', tone: 'medium' };
    case 'LOW':
      return { label: 'Low', tone: 'low' };
    default:
      return { label: tier || 'Unvalidated', tone: 'unknown' };
  }
}

export const TIER_PILL_CLASS: Record<TierTone, string> = {
  high: 'bg-emerald-50 text-emerald-700 border-emerald-200',
  medium: 'bg-amber-50 text-amber-700 border-amber-200',
  low: 'bg-slate-100 text-slate-600 border-slate-200',
  unknown: 'bg-slate-50 text-slate-400 border-slate-200',
};

export const TIER_BAR_CLASS: Record<TierTone, string> = {
  high: 'bg-emerald-500',
  medium: 'bg-amber-500',
  low: 'bg-slate-400',
  unknown: 'bg-slate-300',
};

/* Client-side abstention check. The pipeline doesn't (yet) return a verdict block,
   but the ranking metric is a within-gene percentile, so when the shown models
   cluster tightly at the ceiling the rank order carries little signal — the same
   failure mode GeneTraceAI flags for housekeeping genes like GAPDH. Threshold is
   deliberately conservative and easy to tune. Returns null when the ranking looks
   fine, so callers render nothing in the healthy case. */
const LOW_SEPARATION_RANGE = 0.02;

export function deriveVerdict(results: RankedCellLine[]): Verdict | null {
  if (results.length === 0) {
    return {
      state: 'NO_EVIDENCE',
      headline: 'No models returned for this query.',
      detail: ['The ranking service returned an empty set — there is no cell-line evidence to show.'],
    };
  }

  const metric = (r: RankedCellLine) => (r.mode === 'selectivity' ? r.selectivity : r.score);
  const vals = results.map(metric).sort((a, b) => b - a);
  const top = vals[0];
  const median = vals[Math.floor(vals.length / 2)];
  const range = top - vals[vals.length - 1];

  if (results.length > 1 && top > 0.9 && range < LOW_SEPARATION_RANGE) {
    return {
      state: 'LOW_SEPARATION',
      headline: 'Top models are near-indistinguishable.',
      detail: [
        `The shown models span only ${(range * 100).toFixed(2)} points (top ${top.toFixed(
          4
        )}, median ${median.toFixed(4)}).`,
        'When scores cluster this tightly at the ceiling, rank order carries little signal — treat these as a set of comparable candidates, not a strict ordering.',
      ],
      qualifyRanking: true,
    };
  }

  return null;
}