import React, { useState } from 'react';
import { ALL_GENES } from '../../lib/mockData';
import { useGeneUniverse } from '../../lib/useGeneUniverse';
import { QueryParams } from '../../types';
import { GeneSearch } from '../common/GeneSearch';
import { Search, X, Eraser, BookMarked, Target, Ban, SlidersHorizontal } from 'lucide-react';

interface QueryBuilderProps {
  initialParams: QueryParams;
  onSubmit: (params: QueryParams) => void;
  onOpenReference?: () => void;
}

export const QueryBuilder: React.FC<QueryBuilderProps> = ({ initialParams, onSubmit, onOpenReference }) => {
  const { symbols: geneUniverse, isKnown } = useGeneUniverse();
  const [targets, setTargets] = useState<string[]>(initialParams.targets);
  const [exclusions, setExclusions] = useState<string[]>(initialParams.exclusions);
  const [topK, setTopK] = useState<number>(initialParams.topK);
  const [error, setError] = useState<string | null>(null);

  const removeGene = (gene: string, list: string[], setList: (val: string[]) => void) => {
    setList(list.filter((g) => g !== gene));
  };

  // A gene can't be both a target and an exclusion — adding to one drops it from the other.
  const addTarget = (gene: string) => {
    setExclusions((prev) => prev.filter((g) => g !== gene));
    setTargets((prev) => (prev.includes(gene) ? prev : [...prev, gene]));
  };
  const addExclusion = (gene: string) => {
    setTargets((prev) => prev.filter((g) => g !== gene));
    setExclusions((prev) => (prev.includes(gene) ? prev : [...prev, gene]));
  };

  const handleSubmit = (e: React.FormEvent) => {
    e.preventDefault();
    if (targets.length === 0) {
      setError('Please select at least one target gene.');
      return;
    }
    setError(null);
    onSubmit({ targets, exclusions, lineage: 'Any', topK });
  };

  const nothingSelected = targets.length === 0 && exclusions.length === 0 && topK === 5;
  const handleClear = () => {
    setTargets([]);
    setExclusions([]);
    setTopK(5);
    setError(null);
  };

  // Query mode derived from the current selection (mirrors the routing in fetchRanking).
  const mode =
    targets.length >= 2 && exclusions.length > 0
      ? 'JOINT_SELECTIVITY'
      : targets.length >= 2
      ? 'MULTI'
      : targets.length === 1 && exclusions.length > 0
      ? 'SELECTIVITY'
      : targets.length === 1
      ? 'SINGLE'
      : '—';
  const modeHint =
    mode === 'JOINT_SELECTIVITY'
      ? 'Joint across all targets, high; low in the excluded gene(s)'
      : mode === 'MULTI'
      ? 'Joint ranking across all targets'
      : mode === 'SELECTIVITY'
      ? 'High in the target, low in the excluded gene(s)'
      : mode === 'SINGLE'
      ? 'Single-gene ranking'
      : 'Add a target to begin';

  const Chips: React.FC<{ list: string[]; setList: (v: string[]) => void; tone: 'target' | 'exclude' }> = ({ list, setList, tone }) =>
    list.length === 0 ? null : (
      <div className="flex flex-wrap gap-2 mt-2">
        {list.map((g) => (
          <span
            key={g}
            className={`inline-flex items-center gap-1.5 text-white text-xs font-semibold font-mono pl-2.5 pr-1.5 py-1 rounded-md ${
              tone === 'target' ? 'bg-mulberry-600' : 'bg-rose-600'
            }`}
          >
            {g}
            <button type="button" onClick={() => removeGene(g, list, setList)} aria-label={`Remove ${g}`} className="opacity-80 hover:opacity-100">
              <X className="w-3 h-3" />
            </button>
          </span>
        ))}
      </div>
    );

  const QuickPicks: React.FC<{ list: string[]; setList: (v: string[]) => void; add: (g: string) => void; hide?: string[]; tone: 'target' | 'exclude' }> = ({
    list,
    setList,
    add,
    hide = [],
    tone,
  }) => (
    <div className="flex flex-wrap items-center gap-1.5 mt-3">
      <span className="text-[11px] uppercase tracking-wide text-slate-400 mr-0.5">Common</span>
      {ALL_GENES.filter((g) => !hide.includes(g)).map((gene) => {
        const selected = list.includes(gene);
        const on = tone === 'target' ? 'bg-mulberry-600 text-white border-mulberry-600' : 'bg-rose-600 text-white border-rose-600';
        return (
          <button
            type="button"
            key={gene}
            onClick={() => (selected ? removeGene(gene, list, setList) : add(gene))}
            className={`px-2.5 py-1 rounded-md text-xs font-semibold font-mono border transition ${
              selected ? on : 'bg-white text-slate-600 border-slate-200 hover:border-mulberry-300 hover:text-mulberry-700'
            }`}
          >
            {gene}
          </button>
        );
      })}
    </div>
  );

  return (
    <div className="max-w-5xl mx-auto space-y-6">
      {/* Header */}
      <div>
        <h2 className="text-2xl font-bold text-slate-900 tracking-tight">Find the best cell lines</h2>
        <p className="text-slate-500 mt-1 max-w-2xl leading-relaxed">
          Rank cancer cell-line models for a target gene using integrated multi-omics evidence. Scores are
          percentile-style model outputs — they indicate relative, not absolute, confidence.
        </p>
      </div>

      {error && <div className="p-4 bg-rose-50 border border-rose-200 text-rose-700 rounded-lg text-sm">{error}</div>}

      <form onSubmit={handleSubmit} className="bg-white rounded-2xl shadow-sm border border-slate-200 overflow-hidden">
        {/* thin brand accent */}
        <div className="h-1 bg-gradient-to-r from-mulberry-600 to-mulberry-400" />

        <div className="p-6 grid grid-cols-1 lg:grid-cols-3 gap-x-8 gap-y-6">
          {/* Targets (wide) */}
          <div className="lg:col-span-2">
            <label className="flex items-center gap-2 text-sm font-semibold text-slate-800 mb-2">
              <Target className="w-4 h-4 text-mulberry-600" />
              Target gene(s) <span className="text-rose-500">*</span>
            </label>
            <GeneSearch
              suggestions={geneUniverse}
              isKnown={isKnown}
              exclude={[...targets, ...exclusions]}
              placeholder="Search genes — type to filter, or enter any symbol…"
              onPick={addTarget}
            />
            <p className="text-[11px] text-slate-400 mt-1.5">
              Not sure of a gene symbol?{' '}
              {onOpenReference ? (
                <button type="button" onClick={onOpenReference} className="inline-flex items-center gap-1 text-mulberry-600 hover:text-mulberry-800 font-medium">
                  <BookMarked className="w-3 h-3" /> Look it up in the Reference
                </button>
              ) : (
                <span>Look it up in the Reference tab.</span>
              )}
            </p>
            <Chips list={targets} setList={setTargets} tone="target" />
            <QuickPicks list={targets} setList={setTargets} add={addTarget} hide={exclusions} tone="target" />
          </div>

          {/* Right column: exclusions + count */}
          <div className="space-y-6">
            <div>
              <label className="flex items-center gap-2 text-sm font-semibold text-slate-800 mb-2">
                <Ban className="w-4 h-4 text-slate-400" />
                Exclusion gene(s)
                <span className="text-[11px] font-normal text-slate-400">optional</span>
              </label>
              <GeneSearch
                suggestions={geneUniverse}
                isKnown={isKnown}
                exclude={[...exclusions, ...targets]}
                placeholder="Gene to rank low…"
                onPick={addExclusion}
              />
              <p className="text-[11px] text-slate-400 mt-1.5">Ranks lines high in the target and low here.</p>
              <Chips list={exclusions} setList={setExclusions} tone="exclude" />
              <QuickPicks list={exclusions} setList={setExclusions} add={addExclusion} hide={targets} tone="exclude" />
            </div>

            <div className="pt-5 border-t border-slate-100">
              <div className="flex justify-between items-center mb-1">
                <label className="flex items-center gap-2 text-sm font-semibold text-slate-800">
                  <SlidersHorizontal className="w-4 h-4 text-slate-400" />
                  Recommendations
                </label>
                <span className="text-xs font-bold text-mulberry-700 bg-mulberry-50 border border-mulberry-100 px-2 py-0.5 rounded font-mono">
                  Top {topK}
                </span>
              </div>
              <input
                type="range"
                min={3}
                max={30}
                value={topK}
                onChange={(e) => setTopK(Number(e.target.value))}
                className="w-full h-2 bg-slate-200 rounded-lg appearance-none cursor-pointer accent-mulberry-600 mt-2"
              />
              <div className="flex justify-between text-[11px] text-slate-400 mt-1 font-mono">
                <span>3</span>
                <span>15</span>
                <span>30</span>
              </div>
            </div>
          </div>
        </div>

        {/* Footer bar */}
        <div className="flex flex-col sm:flex-row sm:items-center justify-between gap-3 px-6 py-4 bg-slate-50 border-t border-slate-100">
          <div className="text-xs text-slate-500">
            <span className="text-slate-400">Query mode: </span>
            <span className="font-mono font-bold text-slate-800">{mode}</span>
            <span className="text-slate-400"> · {modeHint}</span>
          </div>
          <div className="flex gap-3">
            <button
              type="button"
              onClick={handleClear}
              disabled={nothingSelected}
              className="flex items-center justify-center gap-2 border border-slate-300 text-slate-600 font-medium px-4 py-2.5 rounded-lg transition hover:bg-white hover:text-slate-800 disabled:opacity-40 disabled:cursor-not-allowed"
            >
              <Eraser className="w-4 h-4" />
              <span>Clear</span>
            </button>
            <button
              type="submit"
              className="flex items-center justify-center gap-2 bg-mulberry-600 hover:bg-mulberry-700 text-white font-semibold px-6 py-2.5 rounded-lg shadow-sm transition"
            >
              <Search className="w-4 h-4" />
              <span>Find Cell Lines</span>
            </button>
          </div>
        </div>
      </form>
    </div>
  );
};
