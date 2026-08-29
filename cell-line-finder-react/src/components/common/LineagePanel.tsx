import React, { useMemo, useState } from 'react';
import { RankedCellLine } from '../../types';
import { ChevronRight, Layers } from 'lucide-react';

interface LineagePanelProps {
  gene: string;
  lines: RankedCellLine[];
  onSelect: (modelId: string) => void;
}

const metricOf = (r: RankedCellLine) =>
  r.mode === 'selectivity' ? r.selectivity : r.mode === 'multi' ? r.jointScore : r.score;

type Group = {
  lineage: string;
  lines: RankedCellLine[];
  bestRank: number;
  topScore: number;
  meanScore: number;
};

export const LineagePanel: React.FC<LineagePanelProps> = ({ gene, lines, onSelect }) => {
  const groups = useMemo<Group[]>(() => {
    const byLineage = new Map<string, RankedCellLine[]>();
    for (const l of lines) {
      const key = l.lineage || 'Unknown';
      if (!byLineage.has(key)) byLineage.set(key, []);
      byLineage.get(key)!.push(l);
    }
    return Array.from(byLineage.entries())
      .map(([lineage, ls]) => {
        const scores = ls.map(metricOf);
        return {
          lineage,
          lines: [...ls].sort((a, b) => a.rank - b.rank),
          bestRank: Math.min(...ls.map((l) => l.rank)),
          topScore: Math.max(...scores),
          meanScore: scores.reduce((s, v) => s + v, 0) / scores.length,
        };
      })
      // Rank lineages by their best-placed line, then by depth (how many they have).
      .sort((a, b) => a.bestRank - b.bestRank || b.lines.length - a.lines.length);
  }, [lines]);

  const [filter, setFilter] = useState<string>('__all__');

  if (groups.length === 0) {
    return (
      <div className="bg-white rounded-xl shadow-sm border border-slate-200 p-8 text-center text-sm text-slate-500">
        No lineage information available for these results.
      </div>
    );
  }

  const shown = filter === '__all__' ? groups : groups.filter((g) => g.lineage === filter);

  return (
    <div className="space-y-4">
      {/* Lineage overview */}
      <div>
        <h3 className="text-lg font-bold text-slate-900">Lineage breakdown for {gene}</h3>
        <p className="text-xs text-slate-500 mt-1">
          Your top {lines.length} lines grouped by tissue lineage, ordered by their best-placed line. Use the
          filters below to focus on a single lineage.
        </p>
      </div>

      {/* Lineage filter */}
      <div className="flex flex-wrap gap-2">
        <FilterPill label={`All (${lines.length})`} active={filter === '__all__'} onClick={() => setFilter('__all__')} />
        {groups.map((g) => (
          <FilterPill
            key={g.lineage}
            label={`${g.lineage} (${g.lines.length})`}
            active={filter === g.lineage}
            onClick={() => setFilter(g.lineage)}
          />
        ))}
      </div>

      {/* Grouped / filtered list */}
      <div className="space-y-4">
        {shown.map((g) => (
          <div key={g.lineage} className="bg-white rounded-xl shadow-sm border border-slate-200 overflow-hidden">
            <div className="px-4 py-3 border-b border-slate-100 flex items-center justify-between">
              <div className="flex items-center gap-2">
                <Layers className="w-4 h-4 text-slate-400" />
                <span className="font-bold text-slate-800 capitalize">{g.lineage}</span>
              </div>
            </div>
            <ul className="divide-y divide-slate-100">
              {g.lines.map((l) => (
                <li key={l.modelId}>
                  <button
                    onClick={() => onSelect(l.modelId)}
                    className="w-full flex items-center justify-between px-4 py-2.5 hover:bg-slate-50 transition text-left"
                  >
                    <div className="flex items-center gap-3">
                      <span className="text-xs font-semibold text-slate-400 w-8">#{l.rank}</span>
                      <div>
                        <span className="font-bold font-mono text-slate-900">{l.modelId}</span>
                        <span className="text-[11px] text-slate-400 ml-2">{l.cellLine}</span>
                      </div>
                    </div>
                    <div className="flex items-center gap-3">
                      <span className="text-xs font-mono text-slate-500">{metricOf(l).toFixed(4)}</span>
                      <ChevronRight className="w-4 h-4 text-slate-300" />
                    </div>
                  </button>
                </li>
              ))}
            </ul>
          </div>
        ))}
      </div>
    </div>
  );
};

const FilterPill: React.FC<{ label: string; active: boolean; onClick: () => void }> = ({
  label,
  active,
  onClick,
}) => (
  <button
    onClick={onClick}
    className={`inline-flex items-center gap-1 text-xs font-medium px-3 py-1.5 rounded-full border transition capitalize ${
      active
        ? 'bg-mulberry-600 text-white border-mulberry-600'
        : 'bg-white text-slate-600 border-slate-200 hover:border-mulberry-300 hover:text-mulberry-700'
    }`}
  >
    {label}
  </button>
);
