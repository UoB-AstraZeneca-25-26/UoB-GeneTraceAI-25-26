import React, { useState } from 'react';
import { SelectivityRanked } from '../../types';

interface SelectivityScatterProps {
  lines: SelectivityRanked[];
  geneHigh: string;
  geneLow: string;
  onSelect: (modelId: string) => void;
}

/* Plots each cell line by its two component scores: x = excluded gene (want LOW),
   y = target gene (want HIGH). Top-left is the selective sweet spot. Uses real
   /exclude data — no per-track detail needed. */
export const SelectivityScatter: React.FC<SelectivityScatterProps> = ({
  lines,
  geneHigh,
  geneLow,
  onSelect,
}) => {
  const W = 620;
  const H = 470;
  const pad = { l: 58, r: 26, t: 26, b: 52 };
  const plotW = W - pad.l - pad.r;
  const plotH = H - pad.t - pad.b;
  const [hover, setHover] = useState<string | null>(null);

  const sx = (v: number) => pad.l + v * plotW; // scoreLow: 0 left, 1 right
  const sy = (v: number) => pad.t + (1 - v) * plotH; // scoreHigh: 1 top, 0 bottom
  const ticks = [0, 0.25, 0.5, 0.75, 1];

  return (
    <div className="bg-white rounded-xl shadow-sm border border-slate-200 p-4">
      <p className="text-xs text-slate-500 mb-2 px-1">
        Each point is a cell line. Up = high <strong>{geneHigh}</strong>; left = low{' '}
        <strong>{geneLow}</strong>. The <span className="text-emerald-600 font-semibold">top-left</span>{' '}
        region is the selective sweet spot.
      </p>
      <svg viewBox={`0 0 ${W} ${H}`} className="w-full" role="img" aria-label={`Selectivity scatter for ${geneHigh} vs ${geneLow}`}>
        {/* sweet-spot shading: high target (top) + low exclusion (left) */}
        <rect x={sx(0)} y={sy(1)} width={plotW / 2} height={plotH / 2} fill="#059669" opacity={0.06} />

        {/* gridlines + axis ticks */}
        {ticks.map((t) => (
          <g key={'gx' + t}>
            <line x1={sx(t)} y1={pad.t} x2={sx(t)} y2={pad.t + plotH} stroke="#eef2f5" strokeWidth={1} />
            <text x={sx(t)} y={pad.t + plotH + 16} textAnchor="middle" style={{ fontSize: 10, fill: '#94a3b8' }}>
              {t}
            </text>
          </g>
        ))}
        {ticks.map((t) => (
          <g key={'gy' + t}>
            <line x1={pad.l} y1={sy(t)} x2={pad.l + plotW} y2={sy(t)} stroke="#eef2f5" strokeWidth={1} />
            <text x={pad.l - 8} y={sy(t) + 3} textAnchor="end" style={{ fontSize: 10, fill: '#94a3b8' }}>
              {t}
            </text>
          </g>
        ))}

        {/* axis labels */}
        <text x={pad.l + plotW / 2} y={H - 8} textAnchor="middle" style={{ fontSize: 12, fill: '#475569', fontWeight: 600 }}>
          {geneLow} score  (want low →)
        </text>
        <text
          x={16}
          y={pad.t + plotH / 2}
          textAnchor="middle"
          transform={`rotate(-90 16 ${pad.t + plotH / 2})`}
          style={{ fontSize: 12, fill: '#475569', fontWeight: 600 }}
        >
          {geneHigh} score  (want high ↑)
        </text>

        {/* points */}
        {lines.map((l) => {
          const active = hover === l.modelId;
          const x = sx(l.scoreLow);
          const y = sy(l.scoreHigh);
          const r = 5 + l.selectivity * 6;
          return (
            <g
              key={l.modelId}
              onMouseEnter={() => setHover(l.modelId)}
              onMouseLeave={() => setHover(null)}
              onClick={() => onSelect(l.modelId)}
              style={{ cursor: 'pointer' }}
            >
              <circle
                cx={x}
                cy={y}
                r={active ? r + 3 : r}
                fill="#4f46e5"
                fillOpacity={0.55 + l.selectivity * 0.4}
                stroke={active ? '#0f172a' : '#fff'}
                strokeWidth={active ? 2 : 1.25}
              />
              {active && (
                <text x={x} y={y - r - 6} textAnchor="middle" style={{ fontSize: 11, fill: '#0f172a', fontWeight: 700 }}>
                  {l.cellLine} · {l.selectivity.toFixed(3)}
                </text>
              )}
            </g>
          );
        })}
      </svg>
    </div>
  );
};
