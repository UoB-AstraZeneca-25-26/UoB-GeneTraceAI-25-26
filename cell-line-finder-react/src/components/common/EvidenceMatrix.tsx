import React, { useEffect, useState } from 'react';
import { CellLineDetail } from '../../types';
import { fetchCellLineDetail, toCellLineDetail } from '../../lib/api';
import { levelBand, TIER_PILL_CLASS } from '../../lib/evidence';
import { Loader2, Check, Minus } from 'lucide-react';

interface EvidenceMatrixProps {
  genes: string[]; // one gene (single) or several (multi / selectivity: target + excluded)
  lines: { modelId: string; cellLine: string; lineage?: string; rank?: number }[];
  onSelect: (modelId: string) => void;
}

const LAYERS = ['Expression', 'Proteomics', 'Mutation', 'Fusion', 'Copy number'];
const key = (gene: string, modelId: string) => `${gene}::${modelId}`;

/**
 * A grid of cell lines (rows) × evidence layers (columns). For a single-gene
 * query it shows that gene's layers; for multi-gene and selectivity queries it
 * repeats the layer block per gene (grouped headers), so you can see — on one
 * screen — that a line is high in the target AND low in the excluded gene.
 * Every cell is categorical: Low/Moderate/High for levels, present/absent for
 * alterations, so it stays scannable as rows and genes grow.
 */
export const EvidenceMatrix: React.FC<EvidenceMatrixProps> = ({ genes, lines, onSelect }) => {
  const [details, setDetails] = useState<Record<string, CellLineDetail>>({});
  const [error, setError] = useState<string | null>(null);

  const geneKey = genes.join(',');
  const lineKey = lines.map((l) => l.modelId).join(',');

  useEffect(() => {
    setDetails({});
    setError(null);
    const ctrl = new AbortController();
    let cancelled = false;
    const pairs = genes.flatMap((g) => lines.map((l) => ({ g, id: l.modelId })));
    Promise.allSettled(
      pairs.map((p) => fetchCellLineDetail(p.g, p.id, ctrl.signal).then((r) => [key(p.g, p.id), toCellLineDetail(r)] as const))
    ).then((settled) => {
      if (cancelled) return;
      const map: Record<string, CellLineDetail> = {};
      settled.forEach((s) => {
        if (s.status === 'fulfilled') map[s.value[0]] = s.value[1];
      });
      if (Object.keys(map).length === 0) setError('Could not load evidence for these lines.');
      setDetails(map);
    });
    return () => {
      cancelled = true;
      ctrl.abort();
    };
    // eslint-disable-next-line react-hooks/exhaustive-deps
  }, [geneKey, lineKey]);

  if (error) {
    return (
      <div className="bg-white rounded-xl shadow-sm border border-slate-200 p-8 text-center text-sm text-slate-500">
        {error}
      </div>
    );
  }

  const total = genes.length * lines.length;
  const loaded = Object.keys(details).length;
  if (loaded === 0) {
    return (
      <div className="bg-white rounded-xl shadow-sm border border-slate-200 p-12 text-center">
        <Loader2 className="w-6 h-6 mx-auto animate-spin text-mulberry-600" />
        <p className="text-xs text-slate-500 mt-3">Loading evidence layers for the top lines…</p>
      </div>
    );
  }

  const multi = genes.length > 1;

  return (
    <div className="bg-white rounded-xl shadow-sm border border-slate-200 p-5">
      <div className="flex items-baseline justify-between gap-3 flex-wrap">
        <div>
          <h3 className="text-base font-bold text-slate-900">
            Evidence grid{multi ? '' : ` for ${genes[0]}`}
          </h3>
          <p className="text-xs text-slate-500 mt-0.5">
            {multi
              ? `Each row is a cell line; layers are shown for ${genes.join(', ')}. Click a line to inspect it.`
              : 'Each row is a cell line; each column an evidence layer. Click a line to inspect it.'}
          </p>
        </div>
        {loaded < total && (
          <span className="text-[11px] text-slate-400 inline-flex items-center gap-1">
            <Loader2 className="w-3 h-3 animate-spin" /> {loaded}/{total}
          </span>
        )}
      </div>

      <div className="overflow-x-auto mt-3">
        <table className="text-sm border-separate border-spacing-0">
          <thead>
            {multi && (
              <tr>
                <th className="sticky left-0 bg-white z-10 border-b border-slate-200" />
                {genes.map((g, gi) => (
                  <th
                    key={g}
                    colSpan={LAYERS.length}
                    className={`text-center text-xs font-bold text-mulberry-700 p-2 border-b border-slate-200 whitespace-nowrap ${
                      gi > 0 ? 'border-l border-slate-200' : ''
                    }`}
                  >
                    {g}
                  </th>
                ))}
              </tr>
            )}
            <tr>
              <th className="sticky left-0 bg-white z-10 text-left text-[11px] uppercase tracking-wide text-slate-400 font-semibold p-2 border-b border-slate-200">
                Cell line
              </th>
              {genes.map((g, gi) =>
                LAYERS.map((c, ci) => (
                  <th
                    key={g + c}
                    className={`text-left text-[11px] uppercase tracking-wide text-slate-400 font-semibold p-2 border-b border-slate-200 whitespace-nowrap ${
                      gi > 0 && ci === 0 ? 'border-l border-slate-200' : ''
                    }`}
                  >
                    {c}
                  </th>
                ))
              )}
            </tr>
          </thead>
          <tbody>
            {lines.map((l) => (
              <tr key={l.modelId} className="hover:bg-slate-50/70 transition">
                <td className="sticky left-0 bg-white z-10 p-2 border-b border-slate-100">
                  <button onClick={() => onSelect(l.modelId)} className="text-left group">
                    <span className="flex items-center gap-2">
                      {l.rank != null && <span className="text-[11px] font-mono text-slate-400">#{l.rank}</span>}
                      <span className="font-semibold font-mono text-slate-800 group-hover:text-mulberry-700 transition">
                        {l.modelId}
                      </span>
                    </span>
                    {l.cellLine && (
                      <span className="block text-[10px] text-slate-400">
                        {l.cellLine}
                        {l.lineage ? ` · ${l.lineage}` : ''}
                      </span>
                    )}
                  </button>
                </td>
                {genes.map((g, gi) => {
                  const d = details[key(g, l.modelId)];
                  if (!d) {
                    return (
                      <td
                        key={g}
                        colSpan={LAYERS.length}
                        className={`p-2 border-b border-slate-100 text-xs text-slate-300 ${gi > 0 ? 'border-l border-slate-200' : ''}`}
                      >
                        <Loader2 className="w-3.5 h-3.5 animate-spin inline" />
                      </td>
                    );
                  }
                  return (
                    <React.Fragment key={g}>
                      <BandCell level={d.expressionLevel} edge={gi > 0} />
                      <BandCell level={d.proteomicsLevel} />
                      <MarkCell yes={(d.pMutation ?? 0) > 0} tag={d.driverAlteration ? 'driver' : undefined} />
                      <MarkCell yes={(d.pFusion ?? 0) > 0} />
                      <MarkCell yes={d.hasCnaAlteration} />
                    </React.Fragment>
                  );
                })}
              </tr>
            ))}
          </tbody>
        </table>
      </div>

      {/* Legend */}
      <div className="flex flex-wrap items-center gap-x-4 gap-y-1.5 mt-3 pt-3 border-t border-slate-100 text-[11px] text-slate-500">
        <span className="font-medium text-slate-400 uppercase tracking-wide text-[10px]">Levels</span>
        {(['Low', 'Moderate', 'High'] as const).map((label) => {
          const tone = label === 'Low' ? 'low' : label === 'Moderate' ? 'medium' : 'high';
          return (
            <span key={label} className={`text-[11px] font-bold px-2 py-0.5 rounded-full border ${TIER_PILL_CLASS[tone]}`}>
              {label}
            </span>
          );
        })}
        <span className="ml-2 inline-flex items-center gap-1 text-emerald-700 font-medium">
          <Check className="w-3.5 h-3.5" /> present
        </span>
        <span className="inline-flex items-center gap-1 text-slate-400 font-medium">
          <Minus className="w-3.5 h-3.5" /> absent
        </span>
        <span className="text-slate-400">· levels are relative to other cell lines for the gene</span>
      </div>
    </div>
  );
};

const BandCell: React.FC<{ level: number | null; edge?: boolean }> = ({ level, edge }) => {
  const border = `p-2 border-b border-slate-100 ${edge ? 'border-l border-slate-200' : ''}`;
  if (level == null) {
    return (
      <td className={border}>
        <span className="text-xs text-slate-300">n/a</span>
      </td>
    );
  }
  const band = levelBand(level);
  return (
    <td className={border}>
      <span className={`inline-flex text-[11px] font-bold px-2.5 py-0.5 rounded-full border ${TIER_PILL_CLASS[band.tone]}`}>
        {band.label}
      </span>
    </td>
  );
};

const MarkCell: React.FC<{ yes: boolean; tag?: string }> = ({ yes, tag }) => (
  <td className="p-2 border-b border-slate-100">
    {yes ? (
      <span className="inline-flex items-center gap-1 text-xs font-semibold text-emerald-700">
        <Check className="w-3.5 h-3.5" />
        {tag && (
          <span className="text-[9px] font-bold text-purple-700 bg-purple-50 border border-purple-200 px-1 py-0.5 rounded">
            {tag}
          </span>
        )}
      </span>
    ) : (
      <span className="inline-flex items-center text-slate-300">
        <Minus className="w-3.5 h-3.5" />
      </span>
    )}
  </td>
);
