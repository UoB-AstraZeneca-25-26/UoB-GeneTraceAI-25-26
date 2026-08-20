import React, { useMemo, useState } from 'react';
import { SingleRanked } from '../../types';
import { lineageColor, pct } from '../../lib/evidence';

interface NetworkGraphProps {
  gene: string;
  lines: SingleRanked[];
  onSelect: (modelId: string) => void;
}

<<<<<<< Updated upstream
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
=======
type Node = { line: SingleRanked; x: number; y: number; a: number; r: number };

/* Clean radial layout: cell lines are spread evenly by angle around the gene so
   the canvas stays balanced, and each node's distance from the centre follows
   its rank (nearest = highest-ranked). A light collision pass only nudges nodes
   that would otherwise touch — no force simulation, so the shape is stable. */
function useLayout(lines: SingleRanked[], W: number, H: number): Node[] {
  return useMemo(() => {
    const cx = W / 2;
    const cy = H / 2;
    const rMin = 96;
    const rMax = Math.min(W, H) / 2 - 58;
    const n = Math.max(lines.length, 1);
    const goldenOffset = -Math.PI / 2; // start at the top

    const nodes: Node[] = lines.map((l, i) => {
      // rank order = array order; spread radius across the rank range
      const t = n === 1 ? 0 : i / (n - 1);
      const r = rMin + t * (rMax - rMin);
      const a = goldenOffset + (i / n) * Math.PI * 2;
      return { line: l, a, r, x: cx + Math.cos(a) * r, y: cy + Math.sin(a) * r };
    });

    // Gentle de-overlap: only separate nodes closer than the min gap.
    for (let it = 0; it < 60; it++) {
>>>>>>> Stashed changes
      for (let i = 0; i < nodes.length; i++) {
        for (let j = i + 1; j < nodes.length; j++) {
          let dx = nodes[j].x - nodes[i].x;
          let dy = nodes[j].y - nodes[i].y;
          const d = Math.hypot(dx, dy) || 0.01;
<<<<<<< Updated upstream
          const min = 60;
=======
          const min = 74;
>>>>>>> Stashed changes
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
<<<<<<< Updated upstream
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
=======
>>>>>>> Stashed changes
    }
    return nodes;
  }, [lines, W, H]);
}

export const NetworkGraph: React.FC<NetworkGraphProps> = ({ gene, lines, onSelect }) => {
  const W = 620;
<<<<<<< Updated upstream
  const H = 470;
=======
  const H = 440;
>>>>>>> Stashed changes
  const cx = W / 2;
  const cy = H / 2;
  const [hover, setHover] = useState<string | null>(null);
  const nodes = useLayout(lines, W, H);

  const lineages = Array.from(new Set(lines.map((l) => l.lineage)));

  return (
<<<<<<< Updated upstream
    <div className="bg-white rounded-xl shadow-sm border border-slate-200 p-4">
      <p className="text-xs text-slate-500 mb-2 px-1">
        Each cell line sits at a distance set by its rank for{' '}
        <strong>{gene}</strong> — nearest the centre is highest-ranked. Node size follows the same
        order; colour marks tissue lineage. Distance is a ranking, not a probability.
      </p>
      <svg viewBox={`0 0 ${W} ${H}`} className="w-full" role="img" aria-label={`Ranking network for ${gene}`}>
=======
    <div className="bg-white rounded-xl shadow-sm border border-slate-200 p-5">
      <p className="text-xs text-slate-500 mb-1 px-1 leading-relaxed">
        Each cell line sits at a distance set by its rank for <strong className="text-slate-700">{gene}</strong> —
        nearest the centre is highest-ranked. Node size follows the same order; colour marks tissue lineage.
        Distance is a ranking, not a probability.
      </p>
      <svg viewBox={`0 0 ${W} ${H}`} className="w-full" role="img" aria-label={`Ranking network for ${gene}`}>
        {/* faint guide rings */}
        {[0.45, 0.72, 1].map((f) => (
          <circle
            key={'ring' + f}
            cx={cx}
            cy={cy}
            r={(Math.min(W, H) / 2 - 58) * f}
            fill="none"
            stroke="#eef2f7"
            strokeWidth={1}
          />
        ))}

        {/* edges gene -> line */}
>>>>>>> Stashed changes
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
<<<<<<< Updated upstream
              opacity={dim ? 0.06 : 0.14 + w * 0.45}
=======
              opacity={dim ? 0.05 : 0.12 + w * 0.4}
>>>>>>> Stashed changes
            />
          );
        })}

<<<<<<< Updated upstream
        <g>
          <circle cx={cx} cy={cy} r={32} fill="#4f46e5" />
          <text x={cx} y={cy - 3} textAnchor="middle" className="fill-white" style={{ fontSize: 14, fontWeight: 700 }}>
            {gene}
          </text>
          <text x={cx} y={cy + 12} textAnchor="middle" className="fill-white" style={{ fontSize: 9, opacity: 0.8 }}>
=======
        {/* gene node */}
        <g>
          <circle cx={cx} cy={cy} r={34} fill="#4f46e5" />
          <text x={cx} y={cy - 3} textAnchor="middle" className="fill-white" style={{ fontSize: 14, fontWeight: 700 }}>
            {gene}
          </text>
          <text x={cx} y={cy + 13} textAnchor="middle" className="fill-white" style={{ fontSize: 9, opacity: 0.75 }}>
>>>>>>> Stashed changes
            gene
          </text>
        </g>

<<<<<<< Updated upstream
        {nodes.map(({ line, x, y }) => {
          const active = hover === line.modelId;
          const dim = hover && !active;
          const r = 8 + ((line.relativeScore ?? 0) / 100) * 12;
          const fill = lineageColor(line.lineage);
          return (
            <g
              key={line.modelId}
              opacity={dim ? 0.4 : 1}
=======
        {/* line nodes */}
        {nodes.map(({ line, x, y, a }) => {
          const active = hover === line.modelId;
          const dim = hover && !active;
          const r = 9 + ((line.relativeScore ?? 0) / 100) * 11;
          const fill = lineageColor(line.lineage);
          // push the label radially outward so it never sits on an edge
          const lx = x + Math.cos(a) * (r + 10);
          const ly = y + Math.sin(a) * (r + 10);
          const anchor = Math.cos(a) > 0.25 ? 'start' : Math.cos(a) < -0.25 ? 'end' : 'middle';
          const dy = Math.sin(a) > 0.35 ? 10 : Math.sin(a) < -0.35 ? -4 : 4;
          return (
            <g
              key={line.modelId}
              opacity={dim ? 0.35 : 1}
>>>>>>> Stashed changes
              onMouseEnter={() => setHover(line.modelId)}
              onMouseLeave={() => setHover(null)}
              onClick={() => onSelect(line.modelId)}
              style={{ cursor: 'pointer' }}
            >
<<<<<<< Updated upstream
              <circle cx={x} cy={y} r={active ? r + 3 : r} fill={fill} stroke={active ? '#0f172a' : '#fff'} strokeWidth={active ? 2 : 1.5} />
=======
              <circle
                cx={x}
                cy={y}
                r={active ? r + 3 : r}
                fill={fill}
                stroke={active ? '#0f172a' : '#fff'}
                strokeWidth={active ? 2.5 : 1.5}
              />
>>>>>>> Stashed changes
              {active && (
                <text x={x} y={y + 1} textAnchor="middle" dominantBaseline="central" className="fill-white" style={{ fontSize: 10, fontWeight: 700 }}>
                  {pct((line.relativeScore ?? 0) / 100)}
                </text>
              )}
<<<<<<< Updated upstream
              <text x={x} y={y + r + 12} textAnchor="middle" style={{ fontSize: 10, fill: '#475569' }}>
=======
              <text
                x={lx}
                y={ly + dy}
                textAnchor={anchor}
                dominantBaseline="central"
                style={{ fontSize: 10.5, fill: active ? '#0f172a' : '#64748b', fontWeight: active ? 700 : 500 }}
              >
>>>>>>> Stashed changes
                {line.cellLine}
              </text>
            </g>
          );
        })}
      </svg>

<<<<<<< Updated upstream
      <div className="flex flex-wrap gap-x-3 gap-y-1 mt-2 px-1 text-xs text-slate-500">
=======
      <div className="flex flex-wrap gap-x-3 gap-y-1 mt-3 pt-3 border-t border-slate-100 px-1 text-xs text-slate-500">
>>>>>>> Stashed changes
        {lineages.map((lin) => (
          <span key={lin} className="inline-flex items-center gap-1.5 capitalize">
            <span className="w-2.5 h-2.5 rounded-full" style={{ background: lineageColor(lin) }} />
            {lin}
          </span>
        ))}
      </div>
    </div>
  );
};
