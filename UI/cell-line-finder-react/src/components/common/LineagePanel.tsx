import React, { useEffect, useState } from 'react';
import { JointSelectivityRanked, MultiRanked, RankedCellLine, ResultMode, SelectivityRanked, SingleRanked } from '../../types';
import {
  contextFields,
  fetchJointSelectivityLineageRanking,
  fetchLineageRanking,
  fetchMultiLineageRanking,
  fetchSelectivityLineageRanking,
  toLineageGrouped,
} from '../../lib/api';
import { ChevronRight, Layers, Loader2 } from 'lucide-react';

interface LineagePanelProps {
  gene: string;              // primary target: the gene for single mode, gene_a for selectivity
  mode: ResultMode;
  genes?: string[];          // full target list, multi mode only
  excludedGenes?: string[];  // exclusion list, selectivity mode only
  onSelect: (modelId: string) => void;
}

const metricOf = (r: RankedCellLine) =>
  r.mode === 'selectivity'
    ? r.selectivity
    : r.mode === 'jointSelectivity'
    ? r.combinedScore
    : r.mode === 'multi'
    ? r.jointScore
    : r.score;

type GroupLine = { ranked: RankedCellLine; rankGlobal: number };
type Group = { lineage: string; lines: GroupLine[] };
type RawRow = { lineage: string; ranked: RankedCellLine; rankGlobal: number };

const groupByLineage = (rows: RawRow[]): Group[] => {
  const byLineage = new Map<string, GroupLine[]>();
  const order: string[] = [];
  for (const { lineage, ranked, rankGlobal } of rows) {
    if (!byLineage.has(lineage)) {
      byLineage.set(lineage, []);
      order.push(lineage);
    }
    byLineage.get(lineage)!.push({ ranked, rankGlobal });
  }
  return order.map((lineage) => ({ lineage, lines: byLineage.get(lineage)! }));
};

/** Server-side, full-panel lineage ranking for every query mode: top 10
 * lineages by their best candidate, top 5 lines within each — not just
 * whatever's on the current results page. The number on the left of each
 * row is rank-within-lineage (1-5); "global #" is this line's rank across
 * every lineage pooled together, for the same gene(s)/query.
 * `unassigned` = additional scoreable lines whose cell line has no lineage
 * on record, so they can't appear grouped here (flagged, not silently
 * dropped). */
export const LineagePanel: React.FC<LineagePanelProps> = ({ gene, mode, genes, excludedGenes, onSelect }) => {
  // Display label reflecting the whole query, not just the primary gene --
  // "BRAF + TP53" for multi, "BRAF vs KRAS, TP53" for selectivity.
  const queryLabel =
    mode === 'jointSelectivity' && genes && genes.length > 1 && excludedGenes && excludedGenes.length > 0
      ? `${genes.join(' + ')} vs ${excludedGenes.join(', ')}`
      : mode === 'multi' && genes && genes.length > 1
      ? genes.join(' + ')
      : mode === 'selectivity' && excludedGenes && excludedGenes.length > 0
      ? `${gene} vs ${excludedGenes.join(', ')}`
      : gene;

  const [groups, setGroups] = useState<Group[] | null>(null);
  const [unassigned, setUnassigned] = useState<number>(0);
  const [error, setError] = useState<string | null>(null);

  useEffect(() => {
    let cancelled = false;
    const ctrl = new AbortController();
    setGroups(null);
    setError(null);

    const onError = (err: unknown) => {
      if (!cancelled) setError(err instanceof Error ? err.message : 'Could not load lineage ranking.');
    };

    if (mode === 'single') {
      fetchLineageRanking(gene, ctrl.signal)
        .then((resp) => {
          if (cancelled) return;
          const data = toLineageGrouped(resp);
          setGroups(
            groupByLineage(
              data.lines.map((l) => ({
                lineage: l.lineage,
                rankGlobal: l.rankGlobal,
                ranked: {
                  mode: 'single',
                  rank: l.rankWithinLineage,
                  modelId: l.modelId,
                  cellLine: l.cellLine,
                  lineage: l.lineage,
                  score: l.coreScore,
                  relativeScore: 0,
                  lineageScore: 0,
                  primaryDisease: l.primaryDisease,
                  subtype: l.subtype,
                  growthPattern: l.growthPattern,
                  sex: l.sex,
                  rrid: l.rrid,
                } satisfies SingleRanked,
              }))
            )
          );
          setUnassigned(data.lineageUnassignedCount);
        })
        .catch(onError);
    } else if (mode === 'multi') {
      if (!genes || genes.length < 2) return;
      fetchMultiLineageRanking(genes, ctrl.signal)
        .then((resp) => {
          if (cancelled) return;
          setGroups(
            groupByLineage(
              resp.lines.map((l) => ({
                lineage: l.lineage,
                rankGlobal: l.rank_global,
                ranked: {
                  mode: 'multi',
                  rank: l.rank_within_lineage,
                  modelId: l.model_id.toUpperCase(),
                  cellLine: l.name && l.name !== 'None' ? l.name.toUpperCase() : l.model_id.toUpperCase(),
                  lineage: l.lineage,
                  jointScore: l.joint_score,
                  genes: resp.genes,
                  geneScores: l.scores,
                  limitingGene: l.limiting_gene,
                  relativeScore: 0,
                  lineageScore: 0,
                  ...contextFields(l.metadata),
                } satisfies MultiRanked,
              }))
            )
          );
          setUnassigned(resp.lineage_unassigned_count);
        })
        .catch(onError);
    } else if (mode === 'selectivity') {
      if (!excludedGenes || excludedGenes.length === 0) return;
      fetchSelectivityLineageRanking(gene, excludedGenes, ctrl.signal)
        .then((resp) => {
          if (cancelled) return;
          setGroups(
            groupByLineage(
              resp.lines.map((l) => ({
                lineage: l.lineage,
                rankGlobal: l.rank_global,
                ranked: {
                  mode: 'selectivity',
                  rank: l.rank_within_lineage,
                  modelId: l.model_id.toUpperCase(),
                  cellLine: l.name && l.name !== 'None' ? l.name.toUpperCase() : l.model_id.toUpperCase(),
                  lineage: l.lineage,
                  selectivity: l.selectivity,
                  geneHigh: resp.gene_a,
                  scoreHigh: l.score_a,
                  excludedGenes: resp.excluded_genes,
                  exclusionScores: l.exclusion_scores,
                  relativeScore: 0,
                  lineageScore: 0,
                  ...contextFields(l.metadata),
                } satisfies SelectivityRanked,
              }))
            )
          );
          setUnassigned(resp.lineage_unassigned_count);
        })
        .catch(onError);
    } else if (mode === 'jointSelectivity') {
      if (!genes || genes.length < 2 || !excludedGenes || excludedGenes.length === 0) return;
      fetchJointSelectivityLineageRanking(genes, excludedGenes, ctrl.signal)
        .then((resp) => {
          if (cancelled) return;
          setGroups(
            groupByLineage(
              resp.lines.map((l) => ({
                lineage: l.lineage,
                rankGlobal: l.rank_global,
                ranked: {
                  mode: 'jointSelectivity',
                  rank: l.rank_within_lineage,
                  modelId: l.model_id.toUpperCase(),
                  cellLine: l.name && l.name !== 'None' ? l.name.toUpperCase() : l.model_id.toUpperCase(),
                  lineage: l.lineage,
                  jointScore: l.joint_score,
                  combinedScore: l.combined_score,
                  genes: resp.genes,
                  geneScores: l.scores,
                  limitingGene: l.limiting_gene,
                  excludedGenes: resp.excluded_genes,
                  exclusionScores: l.exclusion_scores,
                  relativeScore: 0,
                  lineageScore: 0,
                  ...contextFields(l.metadata),
                } satisfies JointSelectivityRanked,
              }))
            )
          );
          setUnassigned(resp.lineage_unassigned_count);
        })
        .catch(onError);
    }

    return () => {
      cancelled = true;
      ctrl.abort();
    };
  }, [gene, mode, genes?.join(','), excludedGenes?.join(',')]);

  const [filter, setFilter] = useState<string>('__all__');

  if (error) {
    return (
      <div className="bg-white rounded-xl shadow-sm border border-slate-200 p-8 text-center text-sm text-slate-500">
        {error}
      </div>
    );
  }

  if (!groups) {
    return (
      <div className="bg-white rounded-xl shadow-sm border border-slate-200 p-12 text-center">
        <Loader2 className="w-6 h-6 mx-auto animate-spin text-indigo-600" />
        <p className="text-xs text-slate-500 mt-3">Ranking every lineage for {queryLabel}…</p>
      </div>
    );
  }

  const totalShown = groups.reduce((n, g) => n + g.lines.length, 0);

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
        <h3 className="text-lg font-bold text-slate-900">Lineage breakdown for {queryLabel}</h3>
        <p className="text-xs text-slate-500 mt-1">
          Top {groups.length} lineages by their single best candidate, top {Math.max(...groups.map((g) => g.lines.length))}{' '}
          lines within each — ranked across the whole panel, not just this page's results. Use the filters below to
          focus on a single lineage.
        </p>
        {unassigned > 0 && (
          <p className="text-xs text-amber-700 mt-1.5">
            {unassigned.toLocaleString()} additional scoreable line{unassigned === 1 ? '' : 's'} have no lineage on
            record, so they can't be grouped here — not silently dropped, just outside this view.
          </p>
        )}
      </div>

      {/* Lineage filter */}
      <div className="flex flex-wrap gap-2">
        <FilterPill label={`All (${totalShown})`} active={filter === '__all__'} onClick={() => setFilter('__all__')} />
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
              {g.lines.map(({ ranked: l, rankGlobal }) => (
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
                      <span className="text-[11px] text-slate-400" title="Rank across every lineage, pooled">
                        global #{rankGlobal}
                      </span>
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
