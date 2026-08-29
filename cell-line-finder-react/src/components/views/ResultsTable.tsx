import React, { useState, useEffect } from 'react';
import { GeneAliasInfo, InspectTarget, QueryParams, RankedCellLine, ResultMeta } from '../../types';
import { EvidenceMatrix } from '../common/EvidenceMatrix';
import { LineagePanel } from '../common/LineagePanel';
import { ChevronRight, AlertCircle, Loader2, Table2, LayoutGrid, Layers } from 'lucide-react';
import { AGENT_API_URL } from '../../config';

interface ResultsTableProps {
  queryParams: QueryParams;
  results: RankedCellLine[];
  meta: ResultMeta | null;
  loading: boolean;
  error: string | null;
  onRetry: () => void;
  onInspect: (t: InspectTarget) => void;
}

export const ResultsTable: React.FC<ResultsTableProps> = ({
  queryParams,
  results,
  meta,
  loading,
  error,
  onRetry,
  onInspect,
}) => {
  const [view, setView] = useState<string>('table');
  const [aliasInfo, setAliasInfo] = useState<Record<string, GeneAliasInfo>>({});
  const [aliasFailed, setAliasFailed] = useState<Record<string, boolean>>({});

  // Gene full-name / alias enrichment — purely additive, never blocks or alters the ranking.
  const aliasGenes = React.useMemo(() => {
    const genes =
      meta?.mode === 'multi'
        ? meta.genes ?? []
        : [meta?.primaryGene ?? queryParams.targets[0], ...(meta?.excludedGenes ?? [])];
    return Array.from(new Set(genes.filter((g): g is string => Boolean(g))));
  }, [meta, queryParams.targets]);

  useEffect(() => {
    aliasGenes.forEach((gene) => {
      fetch(`${AGENT_API_URL}/v1/agent/query`, {
        method: 'POST',
        headers: { 'Content-Type': 'application/json' },
        body: JSON.stringify({ query: `What are the aliases for ${gene}` }),
      })
        .then((res) => res.json())
        .then((data) => {
          if (data?.data?.found) {
            setAliasInfo((prev) => ({ ...prev, [gene]: data.data }));
          } else {
            setAliasFailed((prev) => ({ ...prev, [gene]: true }));
          }
        })
        .catch((err) => {
          // Alias info is enrichment only — fail quietly, but stop showing
          // the loading shimmer so it doesn't spin forever.
          console.error(`Alias fetch failed for ${gene}:`, err);
          setAliasFailed((prev) => ({ ...prev, [gene]: true }));
        });
    });
    // eslint-disable-next-line react-hooks/exhaustive-deps
  }, [aliasGenes.join(',')]);

  // Lineage filter for the table. Declared here (before any early return) so hook
  // order stays constant. When a lineage is picked, the table shows only that
  // lineage's lines, re-ranked by their in-lineage score (highest first).
  const [lineageFilter, setLineageFilter] = useState<string>('all');
  const resultSig = results.map((r) => r.modelId).join(',');
  useEffect(() => {
    setLineageFilter('all');
  }, [resultSig]);

  const lineageCounts = React.useMemo(() => {
    const m = new Map<string, number>();
    results.forEach((r) => m.set(r.lineage, (m.get(r.lineage) ?? 0) + 1));
    return Array.from(m.entries()).sort((a, b) => b[1] - a[1] || a[0].localeCompare(b[0]));
  }, [resultSig]);

  const tableRows = React.useMemo(() => {
    if (lineageFilter === 'all') return results.map((r) => ({ r, displayRank: r.rank, filtered: false }));
    return results
      .filter((r) => r.lineage === lineageFilter)
      .sort((a, b) => b.lineageScore - a.lineageScore)
      .map((r, i) => ({ r, displayRank: i + 1, filtered: true }));
  }, [resultSig, lineageFilter]);

  const pillClass = (active: boolean) =>
    `inline-flex items-center gap-1 text-xs font-medium px-2.5 py-1 rounded-full border transition capitalize ${
      active
        ? 'bg-mulberry-600 text-white border-mulberry-600'
        : 'bg-white text-slate-600 border-slate-200 hover:border-mulberry-300 hover:text-mulberry-700'
    }`;

  if (loading) {
    return (
      <div className="max-w-4xl mx-auto text-center py-24 space-y-4">
        <Loader2 className="w-8 h-8 mx-auto animate-spin text-mulberry-600" />
        <p className="text-sm text-slate-500">
          Querying live ranking for <strong>{queryParams.targets[0]}</strong>
          {queryParams.exclusions.length > 0 && (
            <> vs <strong>{queryParams.exclusions.join(', ')}</strong></>
          )}…
        </p>
      </div>
    );
  }

  if (error) {
    return (
      <div className="max-w-4xl mx-auto text-center py-16 space-y-4">
        <div className="inline-flex p-3 bg-rose-50 text-rose-600 rounded-full">
          <AlertCircle className="w-8 h-8" />
        </div>
        <h2 className="text-xl font-bold text-slate-800">Couldn't load results</h2>
        <p className="text-sm text-slate-500 max-w-md mx-auto">{error}</p>
        <button
          onClick={onRetry}
          className="mt-2 bg-slate-900 hover:bg-slate-800 text-white font-medium px-4 py-2 rounded-lg text-sm transition"
        >
          Retry
        </button>
      </div>
    );
  }

  if (results.length === 0) {
    return (
      <div className="max-w-3xl mx-auto py-16 text-center space-y-3">
        <div className="inline-flex p-3 bg-amber-50 text-amber-600 rounded-full">
          <AlertCircle className="w-8 h-8" />
        </div>
        <h2 className="text-xl font-bold text-slate-800">No Matching Cell Lines</h2>
        <p className="text-sm text-slate-500 max-w-md mx-auto">
          The ranking service returned no cell lines for this query.
        </p>
      </div>
    );
  }

  const mode = results[0].mode;

  // Views available per mode. 'omics' (grid) needs per-line detail fetches.
  const VIEW_META: Record<string, { label: string; icon: React.ElementType }> = {
    table: { label: 'Table', icon: Table2 },
    omics: { label: 'Evidence grid', icon: LayoutGrid },
    lineages: { label: 'Lineages', icon: Layers },
  };
  const views = ['table', 'omics', 'lineages'];
  const activeView = views.includes(view) ? view : 'table';

  // Which single gene the detail endpoint is queried for: the primary target.
  // (multi mode's primaryGene is a joined label, so prefer the first real gene.)
  const inspectGene = meta?.genes?.[0] ?? meta?.primaryGene ?? queryParams.targets[0];

  // Genes whose evidence the grid should show: single -> [gene]; multi -> all
  // targets; selectivity -> target + excluded genes (so you see high-vs-low).
  const gridGenes =
    mode === 'multi'
      ? meta?.genes ?? queryParams.targets
      : mode === 'selectivity'
      ? [meta?.primaryGene ?? queryParams.targets[0], ...(meta?.excludedGenes ?? queryParams.exclusions)]
      : [inspectGene];

  const inspect = (modelId: string) => {
    const row = results.find((r) => r.modelId === modelId);
    onInspect({ gene: inspectGene, modelId, cellLine: row?.cellLine ?? modelId, ranked: row });
  };

  return (
    <div className="space-y-6 max-w-6xl mx-auto">
      {/* Header */}
      <div>
        <h2 className="text-2xl font-bold text-slate-900 font-display tracking-tight">Ranked Recommendations</h2>
        <p className="text-sm text-slate-500 mt-1">
          {meta?.mode === 'selectivity' ? (
            <>
              <strong>Selective for {meta.primaryGene}</strong> against{' '}
              <strong>{meta.excludedGenes?.join(', ')}</strong>
            </>
          ) : meta?.mode === 'multi' ? (
            <>
              <strong>Joint ranking:</strong> {meta.genes?.join(' + ')}
            </>
          ) : (
            <>
              <strong>Target:</strong> {meta?.primaryGene ?? queryParams.targets[0]}
              {meta?.ensg && <span className="text-slate-400"> ({meta.ensg})</span>}
            </>
          )}
          {meta && (
            <> {' '}| Top {results.length} of {meta.total.toLocaleString()} {meta.mode === 'multi' ? 'passing' : 'ranked'}</>
          )}
        </p>

        {/* Gene alias info — agent enrichment, renders below the ranking subtitle */}
        {aliasGenes.some((g) => aliasInfo[g] || !aliasFailed[g]) && (
          <div className="flex flex-col sm:flex-row sm:flex-wrap gap-2 sm:gap-x-6 sm:gap-y-2 mt-2 text-[13px] text-slate-500">
            {aliasGenes.map((gene) => {
              if (aliasInfo[gene]) {
                return (
                  <div key={gene} className="flex items-center flex-wrap gap-1.5">
                    <span className="font-semibold text-slate-800">{gene}</span>
                    <span className="text-slate-300">—</span>
                    <span>{aliasInfo[gene].full_name}</span>
                    {aliasInfo[gene].synonyms && aliasInfo[gene].synonyms!.length > 0 && (
                      <span className="px-2 py-0.5 bg-slate-100 rounded text-xs">
                        aka {aliasInfo[gene].synonyms!.slice(0, 3).join(', ')}
                      </span>
                    )}
                  </div>
                );
              }
              if (aliasFailed[gene]) return null; // enrichment only — simply omit it
              return (
                <div
                  key={gene}
                  className="h-4 w-full sm:w-72 rounded bg-gradient-to-r from-slate-100 via-slate-200 to-slate-100 animate-shimmer"
                />
              );
            })}
          </div>
        )}
      </div>

      {/* View switcher */}
      <div className="flex items-center justify-between gap-3 flex-wrap">
        {views.length > 1 ? (
          <div className="inline-flex bg-slate-100 rounded-lg p-0.5">
            {views.map((v) => {
              const Icon = VIEW_META[v].icon;
              return (
                <button
                  key={v}
                  onClick={() => setView(v)}
                  className={`inline-flex items-center gap-1.5 px-3 py-1.5 rounded-md text-sm font-medium transition ${
                    activeView === v ? 'bg-white text-slate-900 shadow-sm' : 'text-slate-500 hover:text-slate-700'
                  }`}
                >
                  <Icon className="w-4 h-4" /> {VIEW_META[v].label}
                </button>
              );
            })}
          </div>
        ) : (
          <span />
        )}
      </div>

      <div className="grid grid-cols-1 lg:grid-cols-3 gap-6">
        {/* Left panel: table or a visualization */}
        {activeView !== 'table' ? (
          <div className="lg:col-span-2">
            {activeView === 'omics' && (
              <EvidenceMatrix
                genes={gridGenes}
                lines={results.map((r) => ({
                  modelId: r.modelId,
                  cellLine: r.cellLine,
                  lineage: r.lineage,
                  rank: r.rank,
                }))}
                onSelect={inspect}
              />
            )}
            {activeView === 'lineages' && (
              <LineagePanel gene={inspectGene} lines={results} onSelect={inspect} />
            )}
          </div>
        ) : (
        <div className="lg:col-span-2 bg-white rounded-xl shadow-sm border border-slate-200 overflow-hidden">
          {/* Lineage filter */}
          <div className="flex items-center gap-2 flex-wrap px-4 py-3 border-b border-slate-100 bg-slate-50/60">
            <span className="text-xs font-medium text-slate-500 flex items-center gap-1">
              <Layers className="w-3.5 h-3.5 text-slate-400" /> Lineage
            </span>
            <button onClick={() => setLineageFilter('all')} className={pillClass(lineageFilter === 'all')}>
              All <span className="opacity-60">{results.length}</span>
            </button>
            {lineageCounts.map(([lin, count]) => (
              <button key={lin} onClick={() => setLineageFilter(lin)} className={pillClass(lineageFilter === lin)}>
                {lin} <span className="opacity-60">{count}</span>
              </button>
            ))}
            {lineageFilter !== 'all' && (
              <span className="text-[11px] text-mulberry-700 ml-1">↕ re-ranked by in-lineage score</span>
            )}
          </div>
          <table className="w-full text-left text-sm text-slate-700">
            <thead className="bg-slate-50 border-b border-slate-200 text-xs font-semibold text-slate-600 uppercase">
              <tr>
                <th className="px-4 py-3">Rank</th>
                <th className="px-4 py-3">Cell Line</th>
                {mode === 'single' ? (
                  <>
                    <th className="px-4 py-3">Lineage</th>
                    <th className="px-4 py-3">Score <span className="normal-case font-normal text-slate-400">(global / in-lineage)</span></th>
                  </>
                ) : mode === 'selectivity' ? (
                  <>
                    <th className="px-4 py-3">Component scores</th>
                    <th className="px-4 py-3">Selectivity <span className="normal-case font-normal text-slate-400">(global / in-lineage)</span></th>
                  </>
                ) : (
                  <>
                    <th className="px-4 py-3">Per-gene scores</th>
                    <th className="px-4 py-3">Joint <span className="normal-case font-normal text-slate-400">(global / in-lineage)</span></th>
                  </>
                )}
                <th className="px-4 py-3 text-right">Action</th>
              </tr>
            </thead>
            <tbody className="divide-y divide-slate-100">
              {tableRows.map(({ r, displayRank, filtered }) => (
                <tr key={r.modelId} className="hover:bg-slate-50/75 transition">
                  <td className="px-4 py-3 font-semibold text-slate-500">
                    #{displayRank}
                    {filtered && <span className="block text-[10px] font-normal text-slate-400">global #{r.rank}</span>}
                  </td>
                  <td className="px-4 py-3">
                    <div className="flex items-center gap-2">
                      <span className="font-bold font-mono text-slate-900">{r.modelId}</span>
                      {r.mode === 'multi' && r.genes.some((g) => r.geneScores[g] === null) && (
                        <span
                          title="Missing a score for at least one requested gene — joint score is based on partial evidence"
                          className="inline-flex items-center text-[10px] font-bold text-amber-700 bg-amber-50 border border-amber-200 px-1.5 py-0.5 rounded"
                        >
                          PARTIAL
                        </span>
                      )}
                    </div>
                    <span className="text-[11px] text-slate-400">{r.cellLine}</span>
                  </td>

                  {r.mode === 'single' ? (
                    <>
                      <td className="px-4 py-3 capitalize text-slate-600">{r.lineage}</td>
                      <td className="px-4 py-3">
                        <ScorePair global={r.score} lineage={r.lineageScore} />
                      </td>
                    </>
                  ) : r.mode === 'selectivity' ? (
                    <>
                      <td className="px-4 py-3">
                        <div className="flex flex-wrap items-center gap-x-3 gap-y-1 text-xs font-mono">
                          <span className="text-emerald-600" title={`${r.geneHigh} score`}>
                            {r.geneHigh} {r.scoreHigh.toFixed(3)}
                          </span>
                          {r.excludedGenes.map((g) => (
                            <span key={g} className="text-rose-500" title={`${g} score`}>
                              {g} {r.exclusionScores[g].toFixed(3)}
                            </span>
                          ))}
                        </div>
                      </td>
                      <td className="px-4 py-3">
                        <ScorePair global={r.selectivity} lineage={r.lineageScore} />
                      </td>
                    </>
                  ) : (
                    <>
                      <td className="px-4 py-3">
                        <div className="flex flex-wrap items-center gap-x-3 gap-y-1 text-xs font-mono">
                          {r.genes.map((g) => {
                            const v = r.geneScores[g];
                            const limits = r.limitingGene === g;
                            return (
                              <span
                                key={g}
                                className={v === null ? 'text-slate-300' : limits ? 'text-amber-700 font-semibold' : 'text-slate-700'}
                                title={limits ? `${g} limits the joint score` : `${g} score`}
                              >
                                {g} {v === null ? '—' : v.toFixed(3)}
                                {limits && <span className="ml-0.5 text-[9px]">▼</span>}
                              </span>
                            );
                          })}
                        </div>
                      </td>
                      <td className="px-4 py-3">
                        <ScorePair global={r.jointScore} lineage={r.lineageScore} />
                      </td>
                    </>
                  )}

                  <td className="px-4 py-3 text-right">
                    <button
                      onClick={() => inspect(r.modelId)}
                      className="inline-flex items-center text-xs font-semibold text-mulberry-600 hover:text-mulberry-800"
                    >
                      <span>Inspect</span>
                      <ChevronRight className="w-3.5 h-3.5 ml-0.5" />
                    </button>
                  </td>
                </tr>
              ))}
            </tbody>
          </table>
        </div>
        )}

        {/* Contextual explainer for the visual on the left */}
        <div className="space-y-4">
          <ViewExplainer view={activeView} count={results.length} />
        </div>
      </div>
    </div>
  );
};

const ScorePair: React.FC<{ global: number; lineage: number }> = ({ global, lineage }) => (
  <div>
    <div className="text-sm font-mono font-semibold text-slate-800">{global.toFixed(4)}</div>
    <div
      className="text-[11px] font-mono text-slate-400"
      title="Within-lineage score: how this line ranks among shown lines of the same tissue lineage (0–1)"
    >
      lin {lineage.toFixed(3)}
    </div>
  </div>
);

const VIEW_HELP: Record<
  string,
  { title: string; blurb: string; tips: string[] }
> = {
  table: {
    title: 'Ranked results',
    blurb:
      'Every cell line that passed, ordered by score. Each row shows two scores: global (across all lines) and in-lineage (among shown lines of the same tissue). Cell lines are shown by their ACH id, with the common name beneath.',
    tips: ['Hit Inspect on any line to open its full evidence profile.', 'Switch to the Evidence grid to compare layers across lines.'],
  },
  omics: {
    title: 'Evidence grid',
    blurb:
      'A table of cell lines (rows) against evidence layers (columns). For selectivity and multi-gene queries the layers are shown per gene, so you can see a line is high in the target and low in the excluded gene at once.',
    tips: ['Read down a column to compare one layer across lines.', 'Click a line to open its profile.'],
  },
  lineages: {
    title: 'Lineages',
    blurb:
      'Groups the results by tissue lineage so you can compare candidates within and across tissues, ordered by each lineage’s best-placed line.',
    tips: ['Use the filter pills to focus on a single lineage.', 'Click any line to inspect it.'],
  },
};

const ViewExplainer: React.FC<{ view: string; count: number }> = ({ view, count }) => {
  const help = VIEW_HELP[view] ?? VIEW_HELP.table;
  return (
    <div className="bg-white rounded-xl shadow-sm border border-slate-200 p-5">
      <span className="text-[11px] font-semibold uppercase tracking-wider text-mulberry-600">About this view</span>
      <h3 className="text-lg font-bold text-slate-900 mt-1">{help.title}</h3>
      <p className="text-sm text-slate-600 mt-2 leading-relaxed">{help.blurb}</p>
      <ul className="mt-3 space-y-1.5">
        {help.tips.map((t, i) => (
          <li key={i} className="flex gap-2 text-xs text-slate-500 leading-relaxed">
            <span className="text-mulberry-400 shrink-0">→</span>
            <span>{t}</span>
          </li>
        ))}
      </ul>
      <p className="text-[11px] text-slate-400 mt-4 pt-3 border-t border-slate-100">
        Showing your top {count} {count === 1 ? 'line' : 'lines'}.
      </p>
    </div>
  );
};
