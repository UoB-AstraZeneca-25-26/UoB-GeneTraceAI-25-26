import React, { useMemo, useState } from 'react';
import { SingleRanked } from '../../types';
import { TIER_HEX, tierToneOf, confTier, pct } from '../../lib/evidence';

interface NetworkGraphProps {
  gene: string;
  lines: SingleRanked[];
  onSelect: (modelId: string) => void;
}

/* Force-ish layout: each node sits at a distance set by its relativeScore
   (which spreads the near-identical raw scores into a readable gradient), with
   a light collision pass so labels don't overlap. Adapted from GeneTraceAI, but
   colored by confidence tier since the live API returns no lineage. */
function useLayout(lines: SingleRanked[], W: number, H: number) {
  return useMemo(() => {
    const cx = W / 2;
    const cy = H / 2;
    const rMin = 74;
    const rMax = Math.min(W, H) / 2 - 54;
    const n = lines.length;

    const nodes = lines.map((l, i) => {
      const closeness = (l.relativeScore ?? 0) / 100; // 1 = best = nearest
      const a = ((i + 0.5) / n) * Math.PI * 2 + ((i % 3) - 1) * 0.18;
      const target = rMin + (1 - closeness) * (rMax - rMin);
      return {
        line: l,
        x: cx + Math.cos(a) * target,
        y: cy + Math.sin(a) * target,
        target,
      };
    });

    for (let it = 0; it < 200; it++) {
      for (let i = 0; i < nodes.length; i++) {
        for (let j = i + 1; j < nodes.length; j++) {
          let dx = nodes[j].x - nodes[i].x;
          let dy = nodes[j].y - nodes[i].y;
          const d = Math.hypot(dx, dy) || 0.01;
          const min = 60;
          if (d < min) {
            const push = ((min - d) / d) * 0.5;
            dx *= push;
            dy *= push;
            nodes[i].x -= dx;
            nodes[i].y -= dy;
            nodes[j].x += dx;
            nodes[j].y += dy;
          }
        }
      }
      for (const nd of nodes) {
        const dx = nd.x - cx;
        const dy = nd.y - cy;
        const d = Math.hypot(dx, dy) || 0.01;
        const k = ((nd.target - d) / d) * 0.1;
        nd.x += dx * k;
        nd.y += dy * k;
        nd.x = Math.max(36, Math.min(W - 36, nd.x));
        nd.y = Math.max(30, Math.min(H - 30, nd.y));
      }
    }
    return nodes;
  }, [lines, W, H]);
}

export const NetworkGraph: React.FC<NetworkGraphProps> = ({ gene, lines, onSelect }) => {
  const W = 620;
  const H = 470;
  const cx = W / 2;
  const cy = H / 2;
  const [hover, setHover] = useState<string | null>(null);
  const nodes = useLayout(lines, W, H);

  const tones = Array.from(new Set(lines.map((l) => tierToneOf(l.tier))));

  return (
    <div className="bg-white rounded-xl shadow-sm border border-slate-200 p-4">
      <p className="text-xs text-slate-500 mb-2 px-1">
        Each cell line sits at a distance set by its rank for{' '}
        <strong>{gene}</strong> — nearest the centre is highest-ranked. Node size follows the same
        order; colour marks confidence tier. Distance is a ranking, not a probability.
      </p>
      <svg viewBox={`0 0 ${W} ${H}`} className="w-full" role="img" aria-label={`Ranking network for ${gene}`}>
        {nodes.map(({ line, x, y }) => {
          const dim = hover && hover !== line.modelId;
          const w = (line.relativeScore ?? 0) / 100;
          return (
            <line
              key={'e' + line.modelId}
              x1={cx}
              y1={cy}
              x2={x}
              y2={y}
              stroke="#6366f1"
              strokeWidth={1 + w * 3.4}
              strokeLinecap="round"
              opacity={dim ? 0.06 : 0.14 + w * 0.45}
            />
          );
        })}

        <g>
          <circle cx={cx} cy={cy} r={32} fill="#4f46e5" />
          <text x={cx} y={cy - 3} textAnchor="middle" className="fill-white" style={{ fontSize: 14, fontWeight: 700 }}>
            {gene}
          </text>
          <text x={cx} y={cy + 12} textAnchor="middle" className="fill-white" style={{ fontSize: 9, opacity: 0.8 }}>
            gene
          </text>
        </g>

        {nodes.map(({ line, x, y }) => {
          const active = hover === line.modelId;
          const dim = hover && !active;
          const r = 8 + ((line.relativeScore ?? 0) / 100) * 12;
          const fill = TIER_HEX[tierToneOf(line.tier)];
          return (
            <g
              key={line.modelId}
              opacity={dim ? 0.4 : 1}
              onMouseEnter={() => setHover(line.modelId)}
              onMouseLeave={() => setHover(null)}
              onClick={() => onSelect(line.modelId)}
              style={{ cursor: 'pointer' }}
            >
              <circle cx={x} cy={y} r={active ? r + 3 : r} fill={fill} stroke={active ? '#0f172a' : '#fff'} strokeWidth={active ? 2 : 1.5} />
              {active && (
                <text x={x} y={y + 1} textAnchor="middle" dominantBaseline="central" className="fill-white" style={{ fontSize: 10, fontWeight: 700 }}>
                  {pct((line.relativeScore ?? 0) / 100)}
                </text>
              )}
              <text x={x} y={y + r + 12} textAnchor="middle" style={{ fontSize: 10, fill: '#475569' }}>
                {line.cellLine}
              </text>
            </g>
          );
        })}
      </svg>

      <div className="flex flex-wrap gap-3 mt-2 px-1 text-xs text-slate-500">
        {tones.map((tone) => (
          <span key={tone} className="inline-flex items-center gap-1.5">
            <span className="w-2.5 h-2.5 rounded-full" style={{ background: TIER_HEX[tone] }} />
            {confTier(tone === 'high' ? 'HIGH' : tone === 'medium' ? 'MEDIUM' : tone === 'low' ? 'LOW' : '').label}
          </span>
        ))}
      </div>
    </div>
  );
};