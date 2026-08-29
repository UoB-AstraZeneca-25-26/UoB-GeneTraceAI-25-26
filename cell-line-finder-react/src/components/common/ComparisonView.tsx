import React, { useEffect, useMemo, useState } from 'react';
import { CellLineDetail, RankedCellLine } from '../../types';
import { fetchCellLineDetail, toCellLineDetail } from '../../lib/api';
import { confTier, levelBand, TIER_PILL_CLASS } from '../../lib/evidence';
import { Loader2, Check, Minus } from 'lucide-react';

interface ComparisonViewProps {
  gene: string;
  lines: RankedCellLine[];
  onSelect: (modelId: string) => void;
}

const MAX_COMPARE = 4;
const MIN_COMPARE = 2;

/**
 * Side-by-side matrix for a handful of lines. When the top candidates are tied
 * on score, the decision is "which of these do I pick?" — so this fetches each
 * line's evidence detail and lines them up across every layer, marking the
 * best cell in each measured row so trade-offs are obvious at a glance.
 */
export const ComparisonView: React.FC<ComparisonViewProps> = ({ gene, lines, onSelect }) => {
  const [selected, setSelected] = useState<string[]>(() => lines.slice(0, 3).map((l) => l.modelId));
  const [details, setDetails] = useState<Record<string, CellLineDetail>>({});
  const [pending, setPending] = useState(false);

  const byId = useMemo(() => new Map(lines.map((l) => [l.modelId, l])), [lines]);

  // Reset the cache when the gene changes — cached detail is gene-specific.
  useEffect(() => {
    setDetails({});
  }, [gene]);

  // Fetch detail for any selected line we don't already have.
  useEffect(() => {
    const missing = selected.filter((id) => !details[id]);
    if (missing.length === 0) return;
    const ctrl = new AbortController();
    setPending(true);
    Promise.allSettled(
      missing.map((id) =>
        fetchCellLineDetail(gene, id, ctrl.signal).then((r) => [id, toCellLineDetail(r)] as const)
      )
    )
      .then((settled) => {
        const add: Record<string, CellLineDetail> = {};
        settled.forEach((s) => {
          if (s.status === 'fulfilled') add[s.value[0]] = s.value[1];
        });
        if (Object.keys(add).length) setDetails((prev) => ({ ...prev, ...add }));
      })
      .finally(() => setPending(false));
    return () => ctrl.abort();
    // eslint-disable-next-line react-hooks/exhaustive-deps
  }, [selected.join(','), gene]);

  const toggle = (id: string) => {
    setSelected((prev) => {
      if (prev.includes(id)) {
        return prev.length > MIN_COMPARE ? prev.filter((x) => x !== id) : prev;
      }
      return prev.length < MAX_COMPARE ? [...prev, id] : prev;
    });
  };

  // Columns in rank order.
  const cols = useMemo(
    () =>
      selected
        .map((id) => byId.get(id))
        .filter((l): l is RankedCellLine => Boolean(l))
        .sort((a, b) => a.rank - b.rank),
    [selected, byId]
  );

  const ready = cols.every((l) => details[l.modelId]);

  // Best value per measured row (for highlighting), computed over loaded details.
  const best = useMemo(() => {
    const ds = cols.map((l) => details[l.modelId]).filter(Boolean) as CellLineDetail[];
    const maxOf = (pick: (d: CellLineDetail) => number | null) => {
      const vals = ds.map(pick).filter((v): v is number => v != null);
      return vals.length ? Math.max(...vals) : null;
    };
    return {
      score: maxOf((d) => d.score),
      lin: maxOf((d) => d.lineageScore),
      expr: maxOf((d) => d.expressionLevel),
      prot: maxOf((d) => d.proteomicsLevel),
    };
  }, [cols, details]);

  return (
    <div className="bg-white rounded-xl shadow-sm border border-slate-200 p-5 space-y-4">
      <div>
        <h3 className="text-base font-bold text-slate-900">Compare cell lines</h3>
        <p className="text-xs text-slate-500 mt-0.5 leading-relaxed">
          Pick {MIN_COMPARE}–{MAX_COMPARE} lines to line up across every evidence layer for{' '}
          <strong className="text-slate-700">{gene}</strong>. The best cell in each measured row is ringed.
        </p>
      </div>

      {/* Line picker */}
      <div className="flex flex-wrap gap-1.5">
        {lines.map((l) => {
          const on = selected.includes(l.modelId);
          const atCap = !on && selected.length >= MAX_COMPARE;
          return (
            <button
              key={l.modelId}
              onClick={() => toggle(l.modelId)}
              disabled={atCap}
              title={atCap ? `Remove one first (max ${MAX_COMPARE})` : undefined}
              className={`inline-flex items-center gap-1 text-xs font-medium px-2.5 py-1.5 rounded-full border transition ${
                on
                  ? 'bg-indigo-600 text-white border-indigo-600'
                  : atCap
                  ? 'bg-slate-50 text-slate-300 border-slate-200 cursor-not-allowed'
                  : 'bg-white text-slate-600 border-slate-200 hover:border-indigo-300 hover:text-indigo-700'
              }`}
            >
              <span className="font-mono opacity-70">#{l.rank}</span>
              <span className="font-mono">{l.modelId}</span>
            </button>
          );
        })}
      </div>

      {!ready ? (
        <div className="py-16 text-center">
          <Loader2 className="w-6 h-6 mx-auto animate-spin text-indigo-600" />
          <p className="text-xs text-slate-500 mt-3">Loading evidence for the selected lines…</p>
        </div>
      ) : (
        <div className="overflow-x-auto">
          <table className="w-full text-sm border-separate border-spacing-0">
            <thead>
              <tr>
                <th className="sticky left-0 bg-white z-10 text-left align-bottom p-2 w-40" />
                {cols.map((l) => {
                  const d = details[l.modelId];
                  const t = confTier(d.tier);
                  return (
                    <th key={l.modelId} className="p-2 text-left align-bottom min-w-[130px]">
                      <button onClick={() => onSelect(l.modelId)} className="group text-left">
                        <span className="block font-bold font-mono text-slate-900 group-hover:text-indigo-700 transition">
                          {l.modelId}
                        </span>
                        <span className="block text-[10px] text-slate-400">{l.cellLine}</span>
                        <span className="mt-1 flex items-center gap-1.5 flex-wrap">
                          <span className="text-[10px] text-slate-500 font-mono">#{l.rank}</span>
                          <span className={`text-[10px] font-semibold px-1.5 py-0.5 rounded-full border ${TIER_PILL_CLASS[t.tone]}`}>
                            {t.label}
                          </span>
                        </span>
                        <span className="block text-[10px] text-slate-400 capitalize mt-0.5">{l.lineage}</span>
                      </button>
                    </th>
                  );
                })}
              </tr>
            </thead>
            <tbody>
              <SectionRow label="Ranking" span={cols.length + 1} />
              <Row label="Global rank">
                {cols.map((l) => (
                  <Cell key={l.modelId}>
                    <span className="font-mono text-slate-600">#{l.rank}</span>
                  </Cell>
                ))}
              </Row>
              <Row label={`Global score (${gene})`}>
                {cols.map((l) => {
                  const d = details[l.modelId];
                  const isBest = best.score != null && d.score === best.score && cols.length > 1;
                  return (
                    <Cell key={l.modelId} best={isBest}>
                      <span className="font-mono font-semibold text-slate-800">{d.score.toFixed(4)}</span>
                    </Cell>
                  );
                })}
              </Row>
              <Row label="In-lineage score">
                {cols.map((l) => {
                  const d = details[l.modelId];
                  const isBest = best.lin != null && d.lineageScore === best.lin && cols.length > 1;
                  return (
                    <Cell key={l.modelId} best={isBest}>
                      {d.lineageScore != null ? (
                        <span className="font-mono font-semibold text-slate-800">
                          {d.lineageScore.toFixed(3)}
                          {d.lineageRank != null && d.lineageTotal != null && (
                            <span className="text-[10px] text-slate-400 ml-1">#{d.lineageRank}/{d.lineageTotal}</span>
                          )}
                        </span>
                      ) : (
                        <span className="text-xs text-slate-300">—</span>
                      )}
                    </Cell>
                  );
                })}
              </Row>

              <SectionRow label="Measured levels" span={cols.length + 1} />
              <Row label="Expression (RNA)">
                {cols.map((l) => (
                  <LevelCell
                    key={l.modelId}
                    level={details[l.modelId].expressionLevel}
                    best={best.expr != null && details[l.modelId].expressionLevel === best.expr && cols.length > 1}
                  />
                ))}
              </Row>
              <Row label="Proteomics (protein)">
                {cols.map((l) => (
                  <LevelCell
                    key={l.modelId}
                    level={details[l.modelId].proteomicsLevel}
                    best={best.prot != null && details[l.modelId].proteomicsLevel === best.prot && cols.length > 1}
                  />
                ))}
              </Row>

              <SectionRow label="Alterations" span={cols.length + 1} />
              <Row label="Mutation">
                {cols.map((l) => {
                  const d = details[l.modelId];
                  return (
                    <Cell key={l.modelId}>
                      <YesNo yes={d.pMutation != null} extra={d.driverAlteration ? 'driver' : undefined} />
                    </Cell>
                  );
                })}
              </Row>
              <Row label="Fusion">
                {cols.map((l) => (
                  <Cell key={l.modelId}>
                    <YesNo yes={details[l.modelId].pFusion != null} />
                  </Cell>
                ))}
              </Row>
              <Row label="Copy number">
                {cols.map((l) => (
                  <Cell key={l.modelId}>
                    <YesNo yes={details[l.modelId].hasCnaAlteration} />
                  </Cell>
                ))}
              </Row>
            </tbody>
          </table>

          <div className="flex items-center justify-between mt-3 pt-3 border-t border-slate-100">
            <p className="text-[11px] text-slate-400">
              Levels are Low / Moderate / High relative to other cell lines for the gene. Alterations are present / absent.
            </p>
            {pending && <span className="text-[11px] text-slate-400 inline-flex items-center gap-1"><Loader2 className="w-3 h-3 animate-spin" /> updating</span>}
          </div>
        </div>
      )}
    </div>
  );
};

const SectionRow: React.FC<{ label: string; span: number }> = ({ label, span }) => (
  <tr>
    <td
      colSpan={span}
      className="sticky left-0 bg-slate-50 text-[10px] font-semibold uppercase tracking-wide text-slate-400 px-2 py-1.5 border-y border-slate-100"
    >
      {label}
    </td>
  </tr>
);

const Row: React.FC<{ label: string; children: React.ReactNode }> = ({ label, children }) => (
  <tr>
    <td className="sticky left-0 bg-white z-10 text-xs font-medium text-slate-500 p-2 border-b border-slate-100 align-middle">
      {label}
    </td>
    {children}
  </tr>
);

const Cell: React.FC<{ best?: boolean; children: React.ReactNode }> = ({ best, children }) => (
  <td className="p-2 border-b border-slate-100 align-middle">
    <div className={`inline-flex items-center ${best ? 'ring-2 ring-emerald-300 rounded-md px-1.5 py-0.5 bg-emerald-50/40' : ''}`}>
      {children}
    </div>
  </td>
);

const LevelCell: React.FC<{ level: number | null; best?: boolean }> = ({ level, best }) => {
  if (level == null) {
    return (
      <td className="p-2 border-b border-slate-100 align-middle">
        <span className="text-xs text-slate-300">n/a</span>
      </td>
    );
  }
  const band = levelBand(level);
  return (
    <td className="p-2 border-b border-slate-100 align-middle">
      <span
        className={`inline-flex items-center text-[11px] font-bold px-2.5 py-0.5 rounded-full border ${TIER_PILL_CLASS[band.tone]} ${
          best ? 'ring-2 ring-emerald-300' : ''
        }`}
      >
        {band.label}
      </span>
    </td>
  );
};

const YesNo: React.FC<{ yes: boolean; extra?: string }> = ({ yes, extra }) =>
  yes ? (
    <span className="inline-flex items-center gap-1 text-xs font-semibold text-emerald-700">
      <Check className="w-3.5 h-3.5" /> Yes
      {extra && (
        <span className="ml-1 text-[9px] font-bold text-purple-700 bg-purple-50 border border-purple-200 px-1 py-0.5 rounded">
          {extra}
        </span>
      )}
    </span>
  ) : (
    <span className="inline-flex items-center gap-1 text-xs font-medium text-slate-400">
      <Minus className="w-3.5 h-3.5" /> No
    </span>
  );
