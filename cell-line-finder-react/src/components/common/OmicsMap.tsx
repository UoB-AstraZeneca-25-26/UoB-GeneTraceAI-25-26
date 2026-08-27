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
 * shown line (in parallel) and builds the map from real track data. Expression
 * and proteomics are 0-1 levels; mutation / fusion / CNA are categorical.
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
        None of the top lines have measured evidence-layer signal for {gene}.
      </div>
    );
  }

  // ----- layout geometry -----
  const W = 680;
  const rowH = 46;
  const padTop = 54;
  const padBottom = 20;
  const H = Math.max(activeTracks.length, rows.length) * rowH + padTop + padBottom;

  const geneW = Math.max(58, gene.length * 8.5 + 22);
  const xGene = geneW / 2 + 14;
  const xGeneR = xGene + geneW / 2; // gene pill right edge
  const trackW = 116;
  const xTrack = 330;
  const xTrackL = xTrack - trackW / 2;
  const xTrackR = xTrack + trackW / 2;
  const xLine = 592;
  const midY = padTop + ((Math.max(activeTracks.length, rows.length) - 1) * rowH) / 2;

  const colY = (i: number, n: number) => padTop + ((Math.max(activeTracks.length, rows.length) - n) * rowH) / 2 + i * rowH;
  const trackY = (i: number) => colY(i, activeTracks.length);
  const lineY = (i: number) => colY(i, rows.length);

  // Total flow per track (drives the gene→track ribbon width)
  const trackTotals = activeTracks.map((t) =>
    rows.reduce((sum, r) => sum + (r.tracks[t.key]?.score ?? 0), 0)
  );
  const maxTrackTotal = Math.max(...trackTotals, 0.0001);

  const ribbon = (x1: number, y1: number, x2: number, y2: number) => {
    const mx = (x1 + x2) / 2;
    return `M ${x1} ${y1} C ${mx} ${y1}, ${mx} ${y2}, ${x2} ${y2}`;
  };

  const hoveredRow = hover ? rows.find((r) => r.modelId === hover) : null;

  return (
    <div className="bg-white rounded-xl shadow-sm border border-slate-200 p-5">
      <p className="text-xs text-slate-500 mb-1 px-1 leading-relaxed">
        Evidence flow for <strong className="text-slate-700">{gene}</strong>: the gene connects through
        each measured layer to the top cell lines. Ribbon thickness follows the per-layer signal — hover a
        line to trace its evidence, click to inspect.
      </p>

      <svg
        viewBox={`0 0 ${W} ${H}`}
        className="w-full"
        role="img"
        aria-label={`Omics evidence map for ${gene}`}
      >
        {/* column headers */}
        <text x={xGene} y={26} textAnchor="middle" style={{ fontSize: 10, fill: '#94a3b8', fontWeight: 600, letterSpacing: 0.6 }}>
          GENE
        </text>
        <text x={xTrack} y={26} textAnchor="middle" style={{ fontSize: 10, fill: '#94a3b8', fontWeight: 600, letterSpacing: 0.6 }}>
          EVIDENCE LAYERS
        </text>
        <text x={xLine} y={26} textAnchor="middle" style={{ fontSize: 10, fill: '#94a3b8', fontWeight: 600, letterSpacing: 0.6 }}>
          CELL LINES
        </text>

        {/* gene -> track ribbons */}
        {activeTracks.map((t, i) => {
          const y = trackY(i);
          const w = 3 + (trackTotals[i] / maxTrackTotal) * 11;
          return (
            <path
              key={'gt' + t.key}
              d={ribbon(xGeneR, midY, xTrackL, y)}
              fill="none"
              stroke={t.color}
              strokeWidth={w}
              opacity={hover ? 0.1 : 0.26}
              strokeLinecap="round"
            />
          );
        })}

        {/* track -> line ribbons */}
        {rows.map((r) => {
          const yL = lineY(rows.indexOf(r));
          const dim = hover && hover !== r.modelId;
          return activeTracks.map((t, ti) => {
            const s = r.tracks[t.key]?.score ?? 0;
            if (s <= 0) return null;
            const yT = trackY(ti);
            return (
              <path
                key={`tl-${r.modelId}-${t.key}`}
                d={ribbon(xTrackR, yT, xLine - 8, yL)}
                fill="none"
                stroke={t.color}
                strokeWidth={2 + s * 10}
                opacity={dim ? 0.05 : hover === r.modelId ? 0.6 : 0.4}
                strokeLinecap="round"
              />
            );
          });
        })}

        {/* gene node (pill) */}
        <g>
          <rect
            x={xGene - geneW / 2}
            y={midY - 15}
            width={geneW}
            height={30}
            rx={15}
            fill="#4f46e5"
          />
          <text
            x={xGene}
            y={midY + 1}
            textAnchor="middle"
            dominantBaseline="central"
            className="fill-white"
            style={{ fontSize: 12, fontWeight: 700 }}
          >
            {gene}
          </text>
        </g>

        {/* track nodes (pills) */}
        {activeTracks.map((t, i) => {
          const y = trackY(i);
          const strength = pct(trackTotals[i] / Math.max(rows.length, 1));
          return (
            <g key={'tn' + t.key}>
              <rect x={xTrackL} y={y - 14} width={trackW} height={28} rx={8} fill={t.color} opacity={0.12} />
              <rect x={xTrackL} y={y - 14} width={4} height={28} rx={2} fill={t.color} />
              <text x={xTrackL + 14} y={y - 2} dominantBaseline="central" style={{ fontSize: 11, fill: '#334155', fontWeight: 600 }}>
                {t.label}
              </text>
              <text x={xTrackL + 14} y={y + 9} dominantBaseline="central" style={{ fontSize: 8.5, fill: '#94a3b8', fontWeight: 500 }}>
                {t.kind === 'call' ? 'categorical' : `avg ${strength}%`}
              </text>
            </g>
          );
        })}

        {/* line nodes */}
        {rows.map((r, i) => {
          const y = lineY(i);
          const active = hover === r.modelId;
          const nLayers = activeTracks.filter((t) => (r.tracks[t.key]?.score ?? 0) > 0).length;
          return (
            <g
              key={'ln' + r.modelId}
              onMouseEnter={() => setHover(r.modelId)}
              onMouseLeave={() => setHover(null)}
              onClick={() => onSelect(r.modelId)}
              style={{ cursor: 'pointer' }}
            >
              <circle cx={xLine} cy={y} r={active ? 7 : 5} fill={active ? '#0f172a' : '#6366f1'} />
              <text
                x={xLine + 13}
                y={y - 3}
                dominantBaseline="central"
                style={{ fontSize: 11.5, fill: active ? '#0f172a' : '#475569', fontWeight: active ? 700 : 500 }}
              >
                {r.cellLine}
              </text>
              <text x={xLine + 13} y={y + 9} dominantBaseline="central" style={{ fontSize: 8.5, fill: '#94a3b8' }}>
                {nLayers} layer{nLayers === 1 ? '' : 's'}
              </text>
            </g>
          );
        })}
      </svg>

      {/* legend + hover readout */}
      <div className="flex flex-wrap items-center gap-x-4 gap-y-1.5 mt-3 pt-3 border-t border-slate-100 px-1 text-xs text-slate-500">
        {activeTracks.map((t) => (
          <span key={t.key} className="inline-flex items-center gap-1.5">
            <span className="w-2.5 h-2.5 rounded-sm" style={{ background: t.color }} />
            {t.label}
          </span>
        ))}
        {hoveredRow && (
          <span className="ml-auto font-mono text-[11px] text-slate-500">
            <strong className="text-slate-700">{hoveredRow.cellLine}</strong>{' '}
            {activeTracks
              .filter((t) => (hoveredRow.tracks[t.key]?.score ?? 0) > 0)
              .map((t) =>
                t.kind === 'call'
                  ? t.label
                  : `${t.label} ${pct(hoveredRow.tracks[t.key]!.score)}%`
              )
              .join(' · ') || 'no measured signal'}
          </span>
        )}
      </div>
    </div>
  );
};
