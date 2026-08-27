import React, { useState, useEffect } from 'react';
import { GeneAliasInfo, InspectTarget, QueryParams, RankedCellLine, ResultMeta } from '../../types';
import { NetworkGraph } from '../common/NetworkGraph';
import { OmicsMap } from '../common/OmicsMap';
import { LineagePanel } from '../common/LineagePanel';
import { Trophy, ChevronRight, AlertCircle, Loader2, Table2, Network, GitFork, Layers } from 'lucide-react';
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
  const [selectedInDropdown, setSelectedInDropdown] = useState<string>('');
  const [view, setView] = useState<string>('table');
  const [aliasInfo, setAliasInfo] = useState<Record<string, GeneAliasInfo>>({});
  const [aliasFailed, setAliasFailed] = useState<Record<string, boolean>>({});

  useEffect(() => {
    if (results.length > 0) setSelectedInDropdown(results[0].modelId);
  }, [results]);

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

  if (loading) {
    return (
      <div className="max-w-4xl mx-auto text-center py-24 space-y-4">
        <Loader2 className="w-8 h-8 mx-auto animate-spin text-indigo-600" />
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
  const topPick = results[0];

  // Views available per mode. 'omics' (Sankey) needs per-line detail fetches, so
  // it's single/multi only; 'network' (radial rank graph) works for every mode
  // that has a comparable per-line score; 'lineages' works for every mode.
  const VIEW_META: Record<string, { label: string; icon: React.ElementType }> = {
    table: { label: 'Table', icon: Table2 },
    network: { label: 'Network', icon: Network },
    omics: { label: 'Omics map', icon: GitFork },
    lineages: { label: 'Lineages', icon: Layers },
  };
  const views =
    mode === 'single' || mode === 'multi'
      ? ['table', 'network', 'omics', 'lineages']
      : mode === 'selectivity'
      ? ['table', 'network', 'lineages']
      : ['table', 'lineages'];
  const activeView = views.includes(view) ? view : 'table';

  // Which single gene the detail endpoint is queried for: the primary target.
  // (multi mode's primaryGene is a joined label, so prefer the first real gene.)
  const inspectGene = meta?.genes?.[0] ?? meta?.primaryGene ?? queryParams.targets[0];
  const inspect = (modelId: string) => {
    const row = results.find((r) => r.modelId === modelId);
    onInspect({ gene: inspectGene, modelId, cellLine: row?.cellLine ?? modelId });
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
            {activeView === 'network' && (
              <NetworkGraph
                gene={meta?.primaryGene ?? queryParams.targets.join(' + ')}
                lines={results}
                onSelect={inspect}
              />
            )}
            {activeView === 'omics' && (
              <OmicsMap
                gene={inspectGene}
                lines={results.slice(0, 6).map((r) => ({ modelId: r.modelId, cellLine: r.cellLine }))}
                onSelect={inspect}
              />
            )}
            {activeView === 'lineages' && (
              <LineagePanel gene={inspectGene} lines={results} onSelect={inspect} />
            )}
          </div>
        ) : (
        <div className="lg:col-span-2 bg-white rounded-xl shadow-sm border border-slate-200 overflow-hidden">
          <table className="w-full text-left text-sm text-slate-700">
            <thead className="bg-slate-50 border-b border-slate-200 text-xs font-semibold text-slate-600 uppercase">
              <tr>
                <th className="px-4 py-3">Rank</th>
                <th className="px-4 py-3">Cell Line</th>
                {mode === 'single' ? (
                  <>
                    <th className="px-4 py-3">Lineage</th>
                    <th className="px-4 py-3">Score</th>
                  </>
                ) : mode === 'selectivity' ? (
                  <>
                    <th className="px-4 py-3">Component scores</th>
                    <th className="px-4 py-3">Selectivity</th>
                  </>
                ) : (
                  <>
                    <th className="px-4 py-3">Per-gene scores</th>
                    <th className="px-4 py-3">Joint</th>
                  </>
                )}
                <th className="px-4 py-3 text-right">Action</th>
              </tr>
            </thead>
            <tbody className="divide-y divide-slate-100">
              {results.map((r) => (
                <tr key={r.modelId} className="hover:bg-slate-50/75 transition">
                  <td className="px-4 py-3 font-semibold text-slate-500">#{r.rank}</td>
                  <td className="px-4 py-3">
                    <div className="flex items-center gap-2">
                      <span className="font-bold text-slate-900">{r.cellLine}</span>
                      {r.mode === 'multi' && r.genes.some((g) => r.geneScores[g] === null) && (
                        <span
                          title="Missing a score for at least one requested gene — joint score is based on partial evidence"
                          className="inline-flex items-center text-[10px] font-bold text-amber-700 bg-amber-50 border border-amber-200 px-1.5 py-0.5 rounded"
                        >
                          PARTIAL
                        </span>
                      )}
                    </div>
                    <span className="text-[11px] text-slate-400 font-mono">{r.modelId}</span>
                  </td>

                  {r.mode === 'single' ? (
                    <>
                      <td className="px-4 py-3 capitalize text-slate-600">{r.lineage}</td>
                      <td className="px-4 py-3">
                        <ScoreCell relative={r.relativeScore} label={r.score.toFixed(4)} />
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
                        <ScoreCell relative={r.relativeScore} label={r.selectivity.toFixed(4)} />
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
                        <ScoreCell relative={r.relativeScore} label={r.jointScore.toFixed(4)} />
                      </td>
                    </>
                  )}

                  <td className="px-4 py-3 text-right">
                    <button
                      onClick={() => inspect(r.modelId)}
                      className="inline-flex items-center text-xs font-semibold text-indigo-600 hover:text-indigo-800"
                    >
                      <span>Inspect</span>
                      <ChevronRight className="w-3.5 h-3.5 ml-0.5" />
                    </button>
                  </td>
                </tr>
              ))}
            </tbody>
          </table>
          <p className="px-4 py-2.5 border-t border-slate-100 text-[11px] text-slate-400">
            The number is the score. The marker shows each line’s position within the shown set (left = lowest,
            right = highest) — lines in your top {results.length} all score highly, so treat these as comparable candidates.
          </p>
        </div>
        )}

        {/* Top pick + quick dive */}
        <div className="space-y-4">
          <div className="bg-gradient-to-br from-indigo-50 to-indigo-100/50 border border-indigo-100 rounded-xl p-5 space-y-4">
            <div className="flex items-center space-x-2 text-indigo-900 font-bold">
              <Trophy className="w-5 h-5 text-indigo-600" />
              <span>Top Recommendation</span>
            </div>
            <div>
              <h3 className="text-2xl font-black text-slate-900 font-display">{topPick.cellLine}</h3>
              <p className="text-xs text-slate-600 font-mono">{topPick.modelId}</p>
            </div>
            <div className="grid grid-cols-2 gap-2 pt-2 border-t border-indigo-200/60">
              {topPick.mode === 'single' ? (
                <>
                  <Metric label="Lineage" value={topPick.lineage} />
                  <Metric label="Score" value={topPick.score.toFixed(4)} />
                </>
              ) : topPick.mode === 'selectivity' ? (
                <>
                  <Metric label="Selectivity" value={topPick.selectivity.toFixed(4)} />
                  <Metric
                    label={`${topPick.geneHigh}↑ / ${topPick.excludedGenes.join(', ')}↓`}
                    value={`${topPick.scoreHigh.toFixed(2)} / ${topPick.excludedGenes
                      .map((g) => topPick.exclusionScores[g].toFixed(2))
                      .join(', ')}`}
                  />
                </>
              ) : (
                <>
                  <Metric label="Joint score" value={topPick.jointScore.toFixed(4)} />
                  <Metric
                    label="Genes with evidence"
                    value={`${topPick.genes.filter((g) => topPick.geneScores[g] !== null).length} / ${topPick.genes.length}`}
                  />
                </>
              )}
            </div>
          </div>

          <div className="bg-white rounded-xl shadow-sm border border-slate-200 p-5 space-y-3">
            <h4 className="text-xs font-semibold uppercase tracking-wider text-slate-500">Quick Deep Dive</h4>
            <select
              value={selectedInDropdown}
              onChange={(e) => setSelectedInDropdown(e.target.value)}
              className="w-full bg-slate-50 border border-slate-300 text-slate-800 rounded-lg p-2 text-sm focus:ring-2 focus:ring-indigo-500"
            >
              {results.map((r) => (
                <option key={r.modelId} value={r.modelId}>
                  {r.cellLine}
                </option>
              ))}
            </select>
            <button
              onClick={() => inspect(selectedInDropdown)}
              className="w-full bg-slate-900 hover:bg-slate-800 text-white font-medium py-2 rounded-lg text-sm transition"
            >
              Inspect Profile
            </button>
          </div>
        </div>
      </div>
    </div>
  );
};

/**
 * The score is the value that matters, so it leads. The thin track below is a
 * *position* marker — where this line sits within the shown set (left = lowest
 * shown, right = highest) — not a fill, so a high score never looks "empty"
 * just because other shown lines score higher. Lines in the top N all score
 * highly; the marker only conveys their ordering.
 */
const ScoreCell: React.FC<{ relative: number; label: string }> = ({ relative, label }) => {
  const pos = Math.max(0, Math.min(100, relative));
  return (
    <div className="w-28">
      <span className="text-sm font-mono font-semibold text-slate-800">{label}</span>
      <div
        className="relative h-1.5 rounded-full bg-slate-100 mt-1.5"
        title={`Position within the shown lines (${Math.round(pos)}% — left = lowest shown, right = highest)`}
      >
        <div
          className="absolute top-1/2 -translate-y-1/2 -translate-x-1/2 w-2.5 h-2.5 rounded-full bg-indigo-500 border-2 border-white shadow-sm"
          style={{ left: `${pos}%` }}
        />
      </div>
    </div>
  );
};

const Metric: React.FC<{ label: string; value: string }> = ({ label, value }) => (
  <div>
    <span className="text-[11px] text-slate-500 block">{label}</span>
    <span className="text-base font-bold text-slate-800">{value}</span>
  </div>
);
