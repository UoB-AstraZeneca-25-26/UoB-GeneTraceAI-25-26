import React from 'react';
import { Verdict } from '../../types';
import { AlertTriangle, Ban, Info } from 'lucide-react';

const STYLES: Record<Verdict['state'], string> = {
  LOW_SEPARATION: 'bg-amber-50 border-amber-200 text-amber-800',
  NO_EVIDENCE: 'bg-slate-50 border-slate-200 text-slate-600',
  RANKED: 'bg-emerald-50 border-emerald-200 text-emerald-800',
};

const ICONS: Record<Verdict['state'], React.ElementType> = {
  LOW_SEPARATION: AlertTriangle,
  NO_EVIDENCE: Ban,
  RANKED: Info,
};

/**
 * Surfaces an abstention / low-confidence verdict ABOVE the results, rather than
 * drawing a confident-looking table with a buried footnote. Renders nothing when
 * the ranking is healthy (verdict === null).
 */
export const VerdictBanner: React.FC<{ verdict: Verdict | null }> = ({ verdict }) => {
  if (!verdict) return null;
  const Icon = ICONS[verdict.state];

  return (
    <div className={`rounded-lg border p-4 ${STYLES[verdict.state]}`}>
      <div className="flex items-start gap-2">
        <Icon className="w-4 h-4 mt-0.5 shrink-0" />
        <div className="space-y-1.5">
          <div className="text-[10px] font-bold uppercase tracking-wider opacity-70">
            {verdict.state.replace('_', ' ')}
          </div>
          <div className="text-sm font-semibold">{verdict.headline}</div>
          <ul className="text-xs space-y-1 list-disc pl-4 opacity-90">
            {verdict.detail.map((d, i) => (
              <li key={i}>{d}</li>
            ))}
          </ul>
          {verdict.qualifyRanking && (
            <p className="text-xs pt-2 mt-1 border-t border-current/20 opacity-80">
              The ranking below is shown for completeness. It is <strong>not</strong> a recommendation.
            </p>
          )}
        </div>
      </div>
    </div>
  );
};