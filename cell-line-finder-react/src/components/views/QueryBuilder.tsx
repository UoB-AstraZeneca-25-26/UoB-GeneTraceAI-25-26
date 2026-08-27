import React, { useState } from 'react';
import { ALL_GENES } from '../../lib/mockData';
import { QueryParams } from '../../types';
import { GeneSearch } from '../common/GeneSearch';
import { Search, X, Eraser } from 'lucide-react';

interface QueryBuilderProps {
  initialParams: QueryParams;
  onSubmit: (params: QueryParams) => void;
}

export const QueryBuilder: React.FC<QueryBuilderProps> = ({ initialParams, onSubmit }) => {
  const [targets, setTargets] = useState<string[]>(initialParams.targets);
  const [exclusions, setExclusions] = useState<string[]>(initialParams.exclusions);
  const [topK, setTopK] = useState<number>(initialParams.topK);
  const [error, setError] = useState<string | null>(null);

  const toggleGene = (gene: string, list: string[], setList: (val: string[]) => void) => {
    if (list.includes(gene)) setList(list.filter((g) => g !== gene));
    else setList([...list, gene]);
  };

  const addGene = (gene: string, list: string[], setList: (val: string[]) => void) => {
    if (!list.includes(gene)) setList([...list, gene]);
  };

  const removeGene = (gene: string, list: string[], setList: (val: string[]) => void) => {
    setList(list.filter((g) => g !== gene));
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

  // Genes not in the quick-pick grid are shown as removable chips.
  const customTargets = targets.filter((t) => !ALL_GENES.includes(t));
  const customExclusions = exclusions.filter((t) => !ALL_GENES.includes(t));

  return (
    <div className="max-w-4xl mx-auto space-y-8">
      <div>
        <h2 className="text-2xl font-bold text-slate-900 font-display tracking-tight">Find the Best Cell Lines</h2>
        <p className="text-slate-600 mt-1">
          Identify cell lines that best match your biological targets and experimental requirements using
          integrated multi-omics evidence.
        </p>
      </div>

      {error && (
        <div className="p-4 bg-rose-50 border border-rose-200 text-rose-700 rounded-lg text-sm">{error}</div>
      )}

      <form onSubmit={handleSubmit} className="bg-white rounded-xl shadow-sm border border-slate-200 p-6 space-y-6">
        {/* Target Genes */}
        <div>
          <label className="block text-sm font-semibold text-slate-800 mb-2">
            1. Define Biological Target(s) <span className="text-rose-500">*</span>
          </label>
          <div className="mb-2">
            <GeneSearch
              suggestions={ALL_GENES}
              exclude={targets}
              placeholder="Type any gene symbol — e.g. NRAS, CDK4…"
              onPick={(g) => addGene(g, targets, setTargets)}
            />
          </div>
          {customTargets.length > 0 && (
            <div className="flex flex-wrap gap-2 mb-2">
              {customTargets.map((g) => (
                <span key={g} className="inline-flex items-center gap-1 bg-indigo-600 text-white text-xs font-semibold font-mono px-2.5 py-1 rounded-lg">
                  {g}
                  <button type="button" onClick={() => removeGene(g, targets, setTargets)} aria-label={`Remove ${g}`}>
                    <X className="w-3 h-3" />
                  </button>
                </span>
              ))}
            </div>
          )}
          <div className="flex flex-wrap gap-2">
            {ALL_GENES.map((gene) => {
              const selected = targets.includes(gene);
              return (
                <button
                  type="button"
                  key={gene}
                  onClick={() => toggleGene(gene, targets, setTargets)}
                  className={`px-3 py-1.5 rounded-lg text-xs font-semibold font-mono border transition ${
                    selected
                      ? 'bg-indigo-600 text-white border-indigo-600 shadow-sm'
                      : 'bg-slate-50 text-slate-700 border-slate-300 hover:bg-slate-100'
                  }`}
                >
                  {selected ? `✓ ${gene}` : gene}
                </button>
              );
            })}
          </div>
        </div>

        {/* Exclusion Criteria */}
        <div>
          <label className="block text-sm font-semibold text-slate-800 mb-2">2. Define Exclusion Criteria</label>
          <div className="mb-2">
            <GeneSearch
              suggestions={ALL_GENES}
              exclude={[...exclusions, ...targets]}
              placeholder="Type any gene symbol to exclude…"
              onPick={(g) => addGene(g, exclusions, setExclusions)}
            />
          </div>
          {customExclusions.length > 0 && (
            <div className="flex flex-wrap gap-2 mb-2">
              {customExclusions.map((g) => (
                <span key={g} className="inline-flex items-center gap-1 bg-rose-600 text-white text-xs font-semibold font-mono px-2.5 py-1 rounded-lg">
                  {g}
                  <button type="button" onClick={() => removeGene(g, exclusions, setExclusions)} aria-label={`Remove ${g}`}>
                    <X className="w-3 h-3" />
                  </button>
                </span>
              ))}
            </div>
          )}
          <div className="flex flex-wrap gap-2">
            {ALL_GENES.filter((g) => !targets.includes(g)).map((gene) => {
              const selected = exclusions.includes(gene);
              return (
                <button
                  type="button"
                  key={gene}
                  onClick={() => toggleGene(gene, exclusions, setExclusions)}
                  className={`px-3 py-1.5 rounded-lg text-xs font-semibold font-mono border transition ${
                    selected
                      ? 'bg-rose-600 text-white border-rose-600 shadow-sm'
                      : 'bg-slate-50 text-slate-700 border-slate-300 hover:bg-slate-100'
                  }`}
                >
                  {selected ? `✕ ${gene}` : gene}
                </button>
              );
            })}
          </div>
        </div>

        {/* Number of recommendations */}
        <div className="pt-4 border-t border-slate-100">
          <div className="max-w-md">
            <div className="flex justify-between items-center mb-2">
              <label className="text-sm font-semibold text-slate-800">Number of Recommendations</label>
              <span className="text-xs font-bold text-indigo-600 bg-indigo-50 px-2 py-0.5 rounded font-mono">Top {topK}</span>
            </div>
            <input
              type="range"
              min={3}
              max={30}
              value={topK}
              onChange={(e) => setTopK(Number(e.target.value))}
              className="w-full h-2 bg-slate-200 rounded-lg appearance-none cursor-pointer accent-indigo-600 mt-2"
            />
            <div className="flex justify-between text-[11px] text-slate-600 mt-1 font-mono">
              <span>3</span>
              <span>15</span>
              <span>30</span>
            </div>
          </div>
        </div>

        <div className="flex flex-col sm:flex-row gap-3">
          <button
            type="button"
            onClick={handleClear}
            disabled={nothingSelected}
            className="sm:w-44 flex items-center justify-center gap-2 border border-slate-300 text-slate-600 font-medium py-3 rounded-lg transition hover:bg-slate-50 hover:text-slate-800 disabled:opacity-40 disabled:cursor-not-allowed disabled:hover:bg-transparent"
          >
            <Eraser className="w-4 h-4" />
            <span>Clear filters</span>
          </button>
          <button
            type="submit"
            className="flex-1 flex items-center justify-center space-x-2 bg-indigo-600 hover:bg-indigo-700 text-white font-medium py-3 rounded-lg shadow-sm transition"
          >
            <Search className="w-4 h-4" />
            <span>Find Cell Lines</span>
          </button>
        </div>
      </form>
    </div>
  );
};
