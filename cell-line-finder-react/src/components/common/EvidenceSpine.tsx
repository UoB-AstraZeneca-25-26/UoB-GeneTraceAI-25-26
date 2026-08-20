import React from 'react';
import { TrackScores } from '../../types';
import { TRACKS, pct } from '../../lib/evidence';

const segAlpha = (s: number) => 0.28 + 0.72 * s;

interface EvidenceSpineProps {
  tracks?: TrackScores;         // full per-layer detail (future /cellline endpoint)
  nLayers?: number;             // fallback: how many layers fired
  driverAlteration?: boolean;   // fallback: mutation track carries a driver call
}

/**
 * A compact evidence fingerprint per cell line.
 *
 * When per-track detail is available, renders the true six-segment spine — one
 * colored bar per evidence layer, opacity encoding strength, hollow when silent.
 *
 * Until then (the live /gene endpoint only returns a layer *count*), it degrades
 * to an honest N-of layers strip: it shows how many layers fired without claiming
 * which, so no signal is invented.
 */
export const EvidenceSpine: React.FC<EvidenceSpineProps> = ({ tracks, nLayers, driverAlteration }) => {
  if (tracks) {
    return (
      <div className="flex gap-[3px] justify-center" role="img" aria-label="Evidence signature">
        {TRACKS.map((tr) => {
          const s = tracks[tr.key]?.score || 0;
          const on = s > 0;
          return (
            <span
              key={tr.key}
              title={`${tr.label}: ${on ? pct(s) + '%' : 'no signal'}`}
              className={`w-4 h-[9px] rounded-[2px] block ${on ? '' : 'ring-1 ring-inset ring-slate-200'}`}
              style={on ? { background: tr.color, opacity: segAlpha(s) } : undefined}
            />
          );
        })}
      </div>
    );
  }

  const n = Math.max(0, nLayers ?? 0);
  const title = `${n} evidence layer${n === 1 ? '' : 's'}${driverAlteration ? ' · driver alteration' : ''}`;
  return (
    <div className="flex gap-[3px] justify-center" title={title} aria-label={title}>
      {n === 0 ? (
        <span className="w-4 h-[9px] rounded-[2px] block ring-1 ring-inset ring-slate-200" />
      ) : (
        Array.from({ length: Math.min(n, 6) }).map((_, i) => (
          <span key={i} className="w-4 h-[9px] rounded-[2px] block bg-indigo-400" />
        ))
      )}
    </div>
  );
};
