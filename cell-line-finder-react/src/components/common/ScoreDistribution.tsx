import React, { useMemo, useState } from 'react';
import { RankedCellLine } from '../../types';

interface ScoreDistributionProps {
  results: RankedCellLine[];
  total?: number; // full ranked population size, for context
  onSelect: (modelId: string) => void;
}

const metricOf = (r: RankedCellLine) =>
  r.mode === 'selectivity' ? r.selectivity : r.mode === 'multi' ? r.jointScore : r.score;

const metricLabel = (r?: RankedCellLine) =>
  r?.mode === 'selectivity' ? 'selectivity' : r?.mode === 'multi' ? 'joint score' : 'score';

type Placed = { r: RankedCellLine; v: number; x: number; y: number; level: number };

/**
 * A strip/beeswarm plot of the shown lines' scores. Its job is to make the
 * clustering visible: when the top lines are near-tied, you see one dense clump
 * instead of a false ranking. Two scales — "Full (0–1)" shows how high and how
 * bunched the set is against the whole range; "Zoomed" spreads the shown range
 * so you can still read the order.
 */
export const ScoreDistribution: React.FC<ScoreDistributionProps> = ({ results, total, onSelect }) => {
  const [full, setFull] = useState(true);
  const [hover, setHover] = useState<string | null>(null);

  const stats = useMemo(() => {
    const vals = results.map(metricOf);
    const max = Math.max(...vals);
    const min = Math.min(...vals);
    const sorted = [...vals].sort((a, b) => a - b);
    const median = sorted[Math.floor(sorted.length / 2)] ?? 0;
    return { max, min, median, range: max - min };
  }, [results]);

  // Geometry
  const W = 680;
  const H = 220;
  const padL = 24;
  const padR = 24;
  const axisY = 150;
  const plotW = W - padL - padR;

  // Domain: full scale pins to [0,1]; zoomed pads the shown range a touch.
  const [d0, d1] = full
    ? [0, 1]
    : (() => {
        const pad = Math.max(stats.range * 0.15, 0.002);
        return [Math.max(0, stats.min - pad), Math.min(1, stats.max + pad)];
      })();
  const span = d1 - d0 || 1;
  const scaleX = (v: number) => padL + ((v - d0) / span) * plotW;

  // Beeswarm-lite: walk left→right, stack a dot upward when it would collide.
  const placed = useMemo<Placed[]>(() => {
    const R = 6.5;
    const step = R * 1.7;
    const topLimit = 34; // don't let stacks run off the top of the canvas
    const sorted = [...results].sort((a, b) => metricOf(a) - metricOf(b));
    const done: Placed[] = [];
    for (const r of sorted) {
      const v = metricOf(r);
      const x = scaleX(v);
      let level = 0;
      while (done.some((p) => Math.abs(p.x - x) < R * 2 && p.level === level)) level++;
      const y = Math.max(topLimit, axisY - level * step);
      done.push({ r, v, x, y, level });
    }
    return done;
  }, [results, full, d0, d1]);

  // Axis ticks
  const ticks = full ? [0, 0.25, 0.5, 0.75, 1] : [d0, d0 + span / 2, d1];

  const tight = results.length > 1 && stats.max > 0.9 && stats.range < 0.02;
  const topPct = total ? Math.max(0.1, (results.length / total) * 100) : null;

  return (
    <div className="bg-white rounded-xl shadow-sm border border-slate-200 p-5">
      <div className="flex items-start justify-between gap-3 flex-wrap">
        <div>
          <h3 className="text-base font-bold text-slate-900">Score distribution</h3>
          <p className="text-xs text-slate-500 mt-0.5 max-w-md leading-relaxed">
            Where your {results.length} shown lines fall on the {metricLabel(results[0])} scale. Each dot is a
            line — hover to see it, click to inspect.
          </p>
        </div>
        <div className="inline-flex bg-slate-100 rounded-lg p-0.5 text-xs font-medium">
          <button
            onClick={() => setFull(true)}
            className={`px-2.5 py-1 rounded-md transition ${full ? 'bg-white text-slate-900 shadow-sm' : 'text-slate-500 hover:text-slate-700'}`}
          >
            Full (0–1)
          </button>
          <button
            onClick={() => setFull(false)}
            className={`px-2.5 py-1 rounded-md transition ${!full ? 'bg-white text-slate-900 shadow-sm' : 'text-slate-500 hover:text-slate-700'}`}
          >
            Zoomed
          </button>
        </div>
      </div>

      <svg viewBox={`0 0 ${W} ${H}`} className="w-full mt-2" role="img" aria-label="Distribution of shown scores">
        {/* axis */}
        <line x1={padL} y1={axisY} x2={W - padR} y2={axisY} stroke="#e2e8f0" strokeWidth={1} />
        {ticks.map((t) => (
          <g key={t}>
            <line x1={scaleX(t)} y1={axisY} x2={scaleX(t)} y2={axisY + 5} stroke="#cbd5e1" strokeWidth={1} />
            <text x={scaleX(t)} y={axisY + 17} textAnchor="middle" style={{ fontSize: 10, fill: '#94a3b8' }} className="font-mono">
              {t.toFixed(full ? 2 : 4)}
            </text>
          </g>
        ))}

        {/* dots */}
        {placed.map(({ r, v, x, y }) => {
          const active = hover === r.modelId;
          const isTop = r.rank === 1;
          return (
            <g
              key={r.modelId}
              onMouseEnter={() => setHover(r.modelId)}
              onMouseLeave={() => setHover(null)}
              onClick={() => onSelect(r.modelId)}
              style={{ cursor: 'pointer' }}
            >
              <circle
                cx={x}
                cy={y}
                r={active ? 7 : 5.5}
                fill={isTop ? '#4f46e5' : active ? '#4338ca' : '#a5b4fc'}
                stroke={isTop ? '#312e81' : '#fff'}
                strokeWidth={isTop ? 1.5 : 1}
              />
              {active && (
                <>
                  <rect
                    x={Math.min(Math.max(x - 60, 4), W - 124)}
                    y={y - 34}
                    width={120}
                    height={24}
                    rx={4}
                    fill="#0f172a"
                    opacity={0.92}
                  />
                  <text
                    x={Math.min(Math.max(x, 64), W - 64)}
                    y={y - 18}
                    textAnchor="middle"
                    style={{ fontSize: 10, fill: '#fff' }}
                  >
                    #{r.rank} {r.cellLine} · {v.toFixed(4)}
                  </text>
                </>
              )}
            </g>
          );
        })}
      </svg>

      {/* readout */}
      <div className="flex flex-wrap items-center gap-x-5 gap-y-1 text-xs text-slate-500 border-t border-slate-100 pt-3 mt-1">
        <span>Top <span className="font-mono text-slate-700">{stats.max.toFixed(4)}</span></span>
        <span>Median <span className="font-mono text-slate-700">{stats.median.toFixed(4)}</span></span>
        <span>Spread <span className="font-mono text-slate-700">{(stats.range * 100).toFixed(2)}</span> pts</span>
        {topPct && (
          <span>
            Shown lines are the top <span className="font-mono text-slate-700">{results.length}</span> of{' '}
            <span className="font-mono text-slate-700">{total!.toLocaleString()}</span> (top {topPct.toFixed(1)}%)
          </span>
        )}
      </div>

      {tight && (
        <p className="text-[12px] text-amber-700 bg-amber-50 border border-amber-200 rounded-lg p-3 mt-3 leading-relaxed">
          These lines are near-indistinguishable — they span only {(stats.range * 100).toFixed(2)} points at the
          top of the scale. Treat them as a comparable set and let the evidence layers, not the rank order, break the tie.
        </p>
      )}
    </div>
  );
};
