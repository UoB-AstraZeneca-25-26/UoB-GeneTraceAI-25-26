import React, { useState } from 'react';
import { ALL_GENES, ALL_LINEAGES } from '../../lib/mockData';
import { QueryParams } from '../../types';
import { Search, Sparkles } from 'lucide-react';

interface QueryBuilderProps {
  initialParams: QueryParams;
  onSubmit: (params: QueryParams) => void;
}

export const QueryBuilder: React.FC<QueryBuilderProps> = ({ initialParams, onSubmit }) => {
  const [targets, setTargets] = useState<string[]>(initialParams.targets);
  const [exclusions, setExclusions] = useState<string[]>(initialParams.exclusions);
  const [lineage, setLineage] = useState<string>(initialParams.lineage);
  const [topK, setTopK] = useState<number>(initialParams.topK);
  const [error, setError] = useState<string | null>(null);

  const toggleGene = (gene: string, list: string[], setList: (val: string[]) => void) => {
    if (list.includes(gene)) {
      setList(list.filter(g => g !== gene));
    } else {
      setList([...list, gene]);
    }
  };

  const handleSubmit = (e: React.FormEvent) => {
    e.preventDefault();
    if (targets.length === 0) {
      setError("Please select at least one target gene.");
      return;
    }
    setError(null);
    onSubmit({ targets, exclusions, lineage, topK });
  };

  return (
    <div className="max-w-4xl mx-auto space-y-8">
      <div>
        <h2 className="text-2xl font-bold text-slate-900">Find the Best Cell Lines</h2>
        <p className="text-slate-600 mt-1">
          Identify cell lines that best match your biological targets and experimental requirements using integrated multi-omics evidence.
        </p>
      </div>

      {error && (
        <div className="p-4 bg-rose-50 border border-rose-200 text-rose-700 rounded-lg text-sm">
          {error}
        </div>
      )}

      <form onSubmit={handleSubmit} className="bg-white rounded-xl shadow-sm border border-slate-200 p-6 space-y-6">
        {/* Target Genes */}
        <div>
          <label className="block text-sm font-semibold text-slate-800 mb-2">
            1. Define Biological Target(s) <span className="text-rose-500">*</span>
          </label>
          <div className="flex flex-wrap gap-2">
            {ALL_GENES.map(gene => {
              const selected = targets.includes(gene);
              return (
                <button
                  type="button"
                  key={gene}
                  onClick={() => toggleGene(gene, targets, setTargets)}
                  className={`px-3 py-1.5 rounded-lg text-xs font-semibold border transition ${
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
          <label className="block text-sm font-semibold text-slate-800 mb-2">
            2. Define Exclusion Criteria
          </label>
          <div className="flex flex-wrap gap-2">
            {ALL_GENES.filter(g => !targets.includes(g)).map(gene => {
              const selected = exclusions.includes(gene);
              return (
                <button
                  type="button"
                  key={gene}
                  onClick={() => toggleGene(gene, exclusions, setExclusions)}
                  className={`px-3 py-1.5 rounded-lg text-xs font-semibold border transition ${
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

        {/* Biological Context (Lineage) & Top K */}
        <div className="grid grid-cols-1 md:grid-cols-2 gap-6 pt-4 border-t border-slate-100">
          <div>
            <label className="block text-sm font-semibold text-slate-800 mb-2">
              3. Tissue / Lineage Context
            </label>
            <select
              value={lineage}
              onChange={(e) => setLineage(e.target.value)}
              className="w-full bg-slate-50 border border-slate-300 text-slate-800 rounded-lg p-2.5 text-sm focus:ring-2 focus:ring-indigo-500 focus:outline-none"
            >
              {ALL_LINEAGES.map(l => (
                <option key={l} value={l}>{l}</option>
              ))}
            </select>
          </div>

          <div>
            <div className="flex justify-between items-center mb-2">
              <label className="text-sm font-semibold text-slate-800">
                Number of Recommendations
              </label>
              <span className="text-xs font-bold text-indigo-600 bg-indigo-50 px-2 py-0.5 rounded">
                Top {topK}
              </span>
            </div>
            <input
              type="range"
              min={3}
              max={10}
              value={topK}
              onChange={(e) => setTopK(Number(e.target.value))}
              className="w-full h-2 bg-slate-200 rounded-lg appearance-none cursor-pointer accent-indigo-600 mt-2"
            />
            <div className="flex justify-between text-[11px] text-slate-600 mt-1">
              <span>3</span>
              <span>5</span>
              <span>10</span>
            </div>
          </div>
        </div>

        {/* Submit */}
        <button
          type="submit"
          className="w-full flex items-center justify-center space-x-2 bg-indigo-600 hover:bg-indigo-700 text-white font-medium py-3 rounded-lg shadow-sm transition"
        >
          <Search className="w-4 h-4" />
          <span>Find Cell Lines</span>
        </button>
      </form>
    </div>
  );
};