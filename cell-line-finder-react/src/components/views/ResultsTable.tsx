import React, { useState, useEffect } from 'react';
import { QueryParams, RankedCellLine, ResultMeta } from '../../types';
import { deriveVerdict } from '../../lib/evidence';
import { TierPill } from '../common/TierPill';
import { EvidenceSpine } from '../common/EvidenceSpine';
import { VerdictBanner } from '../common/VerdictBanner';
import { Trophy, ChevronRight, AlertCircle, Loader2, Info, Zap } from 'lucide-react';

interface ResultsTableProps {
  queryParams: QueryParams;
  results: RankedCellLine[];
  meta: ResultMeta | null;
  loading: boolean;
  error: string | null;
  onRetry: () => void;
  onSelectCellLine: (cellLine: string) => void;
}

export const ResultsTable: React.FC<ResultsTableProps> = ({
  queryParams,
  results,
  meta,
  loading,
  error,
  onRetry,
  onSelectCellLine,
}) => {
  const [selectedInDropdown, setSelectedInDropdown] = useState<string>('');

  useEffect(() => {
    if (results.length > 0) setSelectedInDropdown(results[0].cellLine);
  }, [results]);

  if (loading) {
    return (
      <div className="max-w-4xl mx-auto text-center py-24 space-y-4">
        <Loader2 className="w-8 h-8 mx-auto animate-spin text-indigo-600" />
        <p className="text-sm text-slate-500">
          Querying live ranking for <strong>{queryParams.targets[0]}</strong>
          {queryParams.exclusions[0] && <> vs <strong>{queryParams.exclusions[0]}</strong></>}…
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

  const verdict = deriveVerdict(results);

  if (results.length === 0) {
    return (
      <div className="max-w-3xl mx-auto py-16">
        <VerdictBanner verdict={verdict} />
      </div>
    );
  }

  const mode = results[0].mode;
  const topPick = results[0];

  const ignoredBits: string[] = [];
  if (queryParams.lineage !== 'Any') ignoredBits.push('lineage filtering');
  if (queryParams.targets.length > 1) ignoredBits.push('additional targets');
  if (queryParams.exclusions.length > 1) ignoredBits.push('additional exclusions');

  return (
    <div className="space-y-6 max-w-6xl mx-auto">
      {/* Header */}
      <div>
        <h2 className="text-2xl font-bold text-slate-900">Ranked Recommendations</h2>
        <p className="text-sm text-slate-500 mt-1">
          {meta?.mode === 'selectivity' ? (
            <>
              <strong>Selective for {meta.primaryGene}</strong> against{' '}
              <strong>{meta.excludedGene}</strong>
            </>
          ) : (
            <>
              <strong>Target:</strong> {meta?.primaryGene ?? queryParams.targets[0]}
              {meta?.ensg && <span className="text-slate-400"> ({meta.ensg})</span>}
            </>
          )}
          {meta && (
            <> {' '}| Top {results.length} of {meta.total.toLocaleString()} ranked</>
          )}
        </p>
      </div>

      {/* Abstention / low-separation verdict, surfaced above the table */}
      <VerdictBanner verdict={verdict} />

      {/* Method banner (selectivity) */}
      {meta?.mode === 'selectivity' && (
        <div className="flex items-start gap-2 p-3 bg-indigo-50 border border-indigo-200 text-indigo-800 rounded-lg text-xs">
          <Info className="w-4 h-4 mt-0.5 shrink-0" />
          <p>
            Ranked by selectivity: <code className="font-mono">{meta.formula}</code> — rewards high{' '}
            {meta.primaryGene} and low {meta.excludedGene}.
          </p>
        </div>
      )}

      {/* What's still not applied */}
      {ignoredBits.length > 0 && (
        <div className="flex items-start gap-2 p-3 bg-amber-50 border border-amber-200 text-amber-800 rounded-lg text-xs">
          <Info className="w-4 h-4 mt-0.5 shrink-0" />
          <p>Not applied by the live endpoint yet: {ignoredBits.join(', ')}.</p>
        </div>
      )}

      <div className="grid grid-cols-1 lg:grid-cols-3 gap-6">
        {/* Table */}
        <div className="lg:col-span-2 bg-white rounded-xl shadow-sm border border-slate-200 overflow-hidden">
          <table className="w-full text-left text-sm text-slate-700">
            <thead className="bg-slate-50 border-b border-slate-200 text-xs font-semibold text-slate-600 uppercase">
              <tr>
                <th className="px-4 py-3">Rank</th>
                <th className="px-4 py-3">Cell Line</th>
                {mode === 'single' ? (
                  <>
                    <th className="px-4 py-3 text-center">Evidence</th>
                    <th className="px-4 py-3">Tier</th>
                    <th className="px-4 py-3">Score</th>
                  </>
                ) : (
                  <>
                    <th className="px-4 py-3">Component scores</th>
                    <th className="px-4 py-3">Selectivity</th>
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
                      {r.mode === 'single' && r.driverAlteration && (
                        <span
                          title="Driver alteration present"
                          className="inline-flex items-center gap-0.5 text-[10px] font-bold text-purple-700 bg-purple-50 border border-purple-200 px-1.5 py-0.5 rounded"
                        >
                          <Zap className="w-3 h-3" /> DRIVER
                        </span>
                      )}
                    </div>
                    <span className="text-[11px] text-slate-400">{r.modelId}</span>
                  </td>

                  {r.mode === 'single' ? (
                    <>
                      <td className="px-4 py-3">
                        <EvidenceSpine
                          tracks={r.tracks}
                          nLayers={r.nLayers}
                          driverAlteration={r.driverAlteration}
                        />
                      </td>
                      <td className="px-4 py-3">
                        <TierPill tier={r.tier} />
                      </td>
                      <td className="px-4 py-3">
                        <ScoreBar width={r.relativeScore} color="bg-indigo-500" label={r.score.toFixed(4)} />
                      </td>
                    </>
                  ) : (
                    <>
                      <td className="px-4 py-3">
                        <div className="flex items-center gap-3 text-xs font-mono">
                          <span className="text-emerald-600" title={`${r.geneHigh} score`}>
                            {r.geneHigh} {r.scoreHigh.toFixed(3)}
                          </span>
                          <span className="text-rose-500" title={`${r.geneLow} score`}>
                            {r.geneLow} {r.scoreLow.toFixed(3)}
                          </span>
                        </div>
                      </td>
                      <td className="px-4 py-3">
                        <ScoreBar width={r.relativeScore} color="bg-indigo-500" label={r.selectivity.toFixed(4)} />
                      </td>
                    </>
                  )}

                  <td className="px-4 py-3 text-right">
                    <button
                      onClick={() => onSelectCellLine(r.cellLine)}
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

          {/* Track legend — only meaningful in single mode */}
          {mode === 'single' && <SpineLegend />}
        </div>

        {/* Top pick + quick dive */}
        <div className="space-y-4">
          <div className="bg-gradient-to-br from-indigo-50 to-indigo-100/50 border border-indigo-100 rounded-xl p-5 space-y-4">
            <div className="flex items-center space-x-2 text-indigo-900 font-bold">
              <Trophy className="w-5 h-5 text-indigo-600" />
              <span>Top Recommendation</span>
            </div>
            <div>
              <h3 className="text-2xl font-black text-slate-900">{topPick.cellLine}</h3>
              <p className="text-xs text-slate-600">{topPick.modelId}</p>
            </div>
            <div className="grid grid-cols-2 gap-2 pt-2 border-t border-indigo-200/60">
              {topPick.mode === 'single' ? (
                <>
                  <div>
                    <span className="text-[11px] text-slate-500 block mb-1">Confidence</span>
                    <TierPill tier={topPick.tier} />
                  </div>
                  <Metric label="Evidence Layers" value={String(topPick.nLayers)} />
                </>
              ) : (
                <>
                  <Metric label="Selectivity" value={topPick.selectivity.toFixed(4)} />
                  <Metric
                    label={`${topPick.geneHigh}↑ / ${topPick.geneLow}↓`}
                    value={`${topPick.scoreHigh.toFixed(2)} / ${topPick.scoreLow.toFixed(2)}`}
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
                <option key={r.modelId} value={r.cellLine}>
                  {r.cellLine}
                </option>
              ))}
            </select>
            <button
              onClick={() => onSelectCellLine(selectedInDropdown)}
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

const ScoreBar: React.FC<{ width: number; color: string; label: string }> = ({ width, color, label }) => (
  <div className="flex items-center space-x-2">
    <div className="w-20 bg-slate-200 rounded-full h-2 overflow-hidden">
      <div className={`h-full rounded-full ${color}`} style={{ width: `${width}%` }} />
    </div>
    <span className="text-xs font-mono text-slate-500">{label}</span>
  </div>
);

const Metric: React.FC<{ label: string; value: string }> = ({ label, value }) => (
  <div>
    <span className="text-[11px] text-slate-500 block">{label}</span>
    <span className="text-base font-bold text-slate-800">{value}</span>
  </div>
);

const SpineLegend: React.FC = () => (
  <div className="px-4 py-3 border-t border-slate-100 text-[11px] text-slate-400">
    Each filled bar is one evidence layer backing the model. Which specific layers
    fired appears here once cell-line detail is available.
  </div>
);