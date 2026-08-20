import React, { useEffect, useMemo, useState } from 'react';
import { TrackScores } from '../../types';
import { fetchCellLineDetail, toCellLineDetail } from '../../lib/api';
import { TRACKS, pct } from '../../lib/evidence';
import { Loader2 } from 'lucide-react';

interface OmicsMapProps {
  gene: string;
  lines: { modelId: string; cellLine: string }[]; // top-N to map
  onSelect: (modelId: string) => void;
}

type LineTracks = { modelId: string; cellLine: string; tracks: TrackScores };

/**
 * Sankey-style flow: gene → evidence layers → cell lines. The per-layer signal
 * isn't in the ranking response, so this fetches the detail endpoint for each
 * shown line (in parallel) and builds the map from real track data. Only the
 * layers the pipeline actually reports (mutation / fusion / CNA) appear.
 */
export const OmicsMap: React.FC<OmicsMapProps> = ({ gene, lines, onSelect }) => {
  const [data, setData] = useState<LineTracks[] | null>(null);
  const [error, setError] = useState<string | null>(null);
  const [hover, setHover] = useState<string | null>(null);

  useEffect(() => {
    let cancelled = false;
    const ctrl = new AbortController();
    setData(null);
    setError(null);
    Promise.allSettled(
      lines.map((l) =>
        fetchCellLineDetail(gene, l.modelId, ctrl.signal).then((r) => ({
          modelId: l.modelId,
          cellLine: l.cellLine,
          tracks: toCellLineDetail(r).tracks,
        }))
      )
    )
      .then((settled) => {
        if (cancelled) return;
        const ok = settled
          .filter((s): s is PromiseFulfilledResult<LineTracks> => s.status === 'fulfilled')
          .map((s) => s.value);
        if (ok.length === 0) setError('Could not load layer detail for these lines.');
        setData(ok);
      })
      .catch(() => {
        if (!cancelled) setError('Could not load layer detail.');
      });
    return () => {
      cancelled = true;
      ctrl.abort();
    };
  }, [gene, JSON.stringify(lines.map((l) => l.modelId))]);

  const layout = useMemo(() => {
    if (!data) return null;
    const activeTracks = TRACKS.filter((t) => data.some((d) => (d.tracks[t.key]?.score ?? 0) > 0));
    return { activeTracks, rows: data };
  }, [data]);

  if (error) {
    return (
      <div className="bg-white rounded-xl shadow-sm border border-slate-200 p-8 text-center text-sm text-slate-500">
        {error}
      </div>
    );
  }

  if (!data || !layout) {
    return (
      <div className="bg-white rounded-xl shadow-sm border border-slate-200 p-12 text-center">
        <Loader2 className="w-6 h-6 mx-auto animate-spin text-indigo-600" />
        <p className="text-xs text-slate-500 mt-3">Loading evidence layers for the top lines…</p>
      </div>
    );
  }

  const { activeTracks, rows } = layout;

  if (activeTracks.length === 0) {
    return (
      <div className="bg-white rounded-xl shadow-sm border border-slate-200 p-8 text-center text-sm text-slate-500">
        None of the top lines have alteration-layer signal (mutation / fusion / CNA) for {gene}.
      </div>
    );
  }

  // ----- layout geometry -----
  const W = 660;
  const rowH = 34;
  const H = Math.max(activeTracks.length, rows.length) * rowH + 60;
  const xGene = 70;
  const xTrack = 320;
  const xLine = 560;
  const geneY = H / 2;

  const trackY = (i: number, n: number) => (H - (n - 1) * rowH) / 2 + i * rowH;
  const lineY = (i: number, n: number) => (H - (n - 1) * rowH) / 2 + i * rowH;

  // Total flow per track (for gene->track ribbon width)
  const trackTotals = activeTracks.map((t) =>
    rows.reduce((sum, r) => sum + (r.tracks[t.key]?.score ?? 0), 0)
  );
  const maxTrackTotal = Math.max(...trackTotals, 0.0001);

  const ribbon = (x1: number, y1: number, x2: number, y2: number) => {
    const mx = (x1 + x2) / 2;
    return `M ${x1} ${y1} C ${mx} ${y1}, ${mx} ${y2}, ${x2} ${y2}`;
  };

  return (
    <div className="bg-white rounded-xl shadow-sm border border-slate-200 p-4">
      <p className="text-xs text-slate-500 mb-2 px-1">
        Evidence flow for <strong>{gene}</strong>: the gene connects through each reported alteration
        layer to the top cell lines. Ribbon thickness follows the per-layer signal. Hover a line to
        trace its evidence; click to inspect.
      </p>
      <svg viewBox={`0 0 ${W} ${H}`} className="w-full" role="img" aria-label={`Omics evidence map for ${gene}`}>
        {/* gene -> track ribbons */}
        {activeTracks.map((t, i) => {
          const y = trackY(i, activeTracks.length);
          const w = 2 + (trackTotals[i] / maxTrackTotal) * 10;
          return (
            <path
              key={'gt' + t.key}
              d={ribbon(xGene + 14, geneY, xTrack - 52, y)}
              fill="none"
              stroke={t.color}
              strokeWidth={w}
              opacity={hover ? 0.12 : 0.28}
              strokeLinecap="round"
            />
          );
        })}

        {/* track -> line ribbons */}
        {rows.map((r, li) => {
          const yL = lineY(li, rows.length);
          const dim = hover && hover !== r.modelId;
          return activeTracks.map((t, ti) => {
            const s = r.tracks[t.key]?.score ?? 0;
            if (s <= 0) return null;
            const yT = trackY(ti, activeTracks.length);
            return (
              <path
                key={`tl-${r.modelId}-${t.key}`}
                d={ribbon(xTrack + 52, yT, xLine - 6, yL)}
                fill="none"
                stroke={t.color}
                strokeWidth={2 + s * 9}
                opacity={dim ? 0.06 : 0.45}
                strokeLinecap="round"
              />
            );
          });
        })}

        {/* gene node */}
        <g>
          <circle cx={xGene} cy={geneY} r={16} fill="#4f46e5" />
          <text x={xGene} y={geneY + 1} textAnchor="middle" dominantBaseline="central" className="fill-white" style={{ fontSize: 11, fontWeight: 700 }}>
            {gene.length > 5 ? gene.slice(0, 5) : gene}
          </text>
        </g>

        {/* track nodes */}
        {activeTracks.map((t, i) => {
          const y = trackY(i, activeTracks.length);
          return (
            <g key={'tn' + t.key}>
              <rect x={xTrack - 52} y={y - 11} width={104} height={22} rx={6} fill={t.color} opacity={0.14} />
              <rect x={xTrack - 52} y={y - 11} width={4} height={22} rx={2} fill={t.color} />
              <text x={xTrack} y={y + 1} textAnchor="middle" dominantBaseline="central" style={{ fontSize: 11, fill: '#334155', fontWeight: 600 }}>
                {t.label}
              </text>
            </g>
          );
        })}

        {/* line nodes */}
        {rows.map((r, i) => {
          const y = lineY(i, rows.length);
          const active = hover === r.modelId;
          return (
            <g
              key={'ln' + r.modelId}
              onMouseEnter={() => setHover(r.modelId)}
              onMouseLeave={() => setHover(null)}
              onClick={() => onSelect(r.modelId)}
              style={{ cursor: 'pointer' }}
            >
              <circle cx={xLine} cy={y} r={active ? 7 : 5} fill={active ? '#0f172a' : '#6366f1'} />
              <text x={xLine + 12} y={y + 1} dominantBaseline="central" style={{ fontSize: 11, fill: active ? '#0f172a' : '#475569', fontWeight: active ? 700 : 500 }}>
                {r.cellLine}
              </text>
            </g>
          );
        })}
      </svg>

      <div className="flex flex-wrap gap-3 mt-2 px-1 text-xs text-slate-500">
        {activeTracks.map((t) => (
          <span key={t.key} className="inline-flex items-center gap-1.5">
            <span className="w-2.5 h-2.5 rounded-sm" style={{ background: t.color }} />
            {t.label}
          </span>
        ))}
        {hover && (
          <span className="ml-auto font-mono text-slate-400">
            {rows.find((r) => r.modelId === hover)?.cellLine}:{' '}
            {activeTracks
              .filter((t) => (rows.find((r) => r.modelId === hover)?.tracks[t.key]?.score ?? 0) > 0)
              .map((t) => `${t.label} ${pct(rows.find((r) => r.modelId === hover)!.tracks[t.key]!.score)}%`)
              .join(' · ') || 'no alteration signal'}
          </span>
        )}
      </div>
    </div>
  );
};