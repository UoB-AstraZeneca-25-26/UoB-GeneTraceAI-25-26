import React, { useMemo, useState } from 'react';
import { CELL_LINE_DIRECTORY } from '../../lib/mockApi';
import { GENE_REF } from '../../lib/mockData';
import { Search, Dna, FlaskConical, Info } from 'lucide-react';

type Tab = 'lines' | 'genes';

/**
 * A lookup table so users can find the ACH id for a cell-line name (and vice
 * versa), and confirm gene symbols. This is a representative preview; the full
 * DepMap / gene tables will be wired in behind the same search UI later.
 */
export const ReferencePage: React.FC = () => {
  const [tab, setTab] = useState<Tab>('lines');
  const [q, setQ] = useState('');

  const query = q.trim().toLowerCase();

  const lines = useMemo(() => {
    if (!query) return CELL_LINE_DIRECTORY;
    return CELL_LINE_DIRECTORY.filter(
      (l) =>
        l.ach.toLowerCase().includes(query) ||
        l.name.toLowerCase().includes(query) ||
        l.lineage.toLowerCase().includes(query) ||
        l.disease.toLowerCase().includes(query)
    );
  }, [query]);

  const genes = useMemo(() => {
    if (!query) return GENE_REF;
    return GENE_REF.filter(
      (g) =>
        g.symbol.toLowerCase().includes(query) ||
        g.name.toLowerCase().includes(query) ||
        g.ensg.toLowerCase().includes(query)
    );
  }, [query]);

  return (
    <div className="max-w-4xl mx-auto space-y-6">
      <div>
        <h2 className="text-2xl font-bold text-slate-900 font-display tracking-tight">Reference lookup</h2>
        <p className="text-slate-600 mt-1">
          Find a cell line’s <span className="font-mono text-slate-700">ACH</span> id from its name (or the other
          way round), and confirm gene symbols before you search.
        </p>
      </div>

      <div className="flex items-start gap-2 text-xs text-slate-500 bg-mulberry-50/60 border border-mulberry-100 rounded-lg p-3">
        <Info className="w-4 h-4 text-mulberry-500 shrink-0 mt-0.5" />
        <span>
          This is a representative preview. The full cell-line and gene tables will be connected here later — the
          search box will work the same way against the complete list.
        </span>
      </div>

      {/* Tabs */}
      <div className="inline-flex bg-slate-100 rounded-lg p-0.5">
        <button
          onClick={() => setTab('lines')}
          className={`inline-flex items-center gap-1.5 px-3 py-1.5 rounded-md text-sm font-medium transition ${
            tab === 'lines' ? 'bg-white text-slate-900 shadow-sm' : 'text-slate-500 hover:text-slate-700'
          }`}
        >
          <FlaskConical className="w-4 h-4" /> Cell lines
        </button>
        <button
          onClick={() => setTab('genes')}
          className={`inline-flex items-center gap-1.5 px-3 py-1.5 rounded-md text-sm font-medium transition ${
            tab === 'genes' ? 'bg-white text-slate-900 shadow-sm' : 'text-slate-500 hover:text-slate-700'
          }`}
        >
          <Dna className="w-4 h-4" /> Genes
        </button>
      </div>

      {/* Search */}
      <div className="relative">
        <Search className="w-4 h-4 absolute left-3 top-1/2 -translate-y-1/2 text-slate-400 pointer-events-none" />
        <input
          value={q}
          onChange={(e) => setQ(e.target.value)}
          placeholder={tab === 'lines' ? 'Search by ACH id, name, lineage or disease…' : 'Search by gene symbol or name…'}
          className="w-full bg-white border border-slate-300 rounded-lg pl-9 pr-3 py-2.5 text-sm focus:ring-2 focus:ring-mulberry-500 focus:outline-none"
        />
      </div>

      {/* Table */}
      <div className="bg-white rounded-xl shadow-sm border border-slate-200 overflow-hidden">
        {tab === 'lines' ? (
          <table className="w-full text-left text-sm">
            <thead className="bg-slate-50 border-b border-slate-200 text-xs font-semibold text-slate-600 uppercase">
              <tr>
                <th className="px-4 py-3">ACH id</th>
                <th className="px-4 py-3">Name</th>
                <th className="px-4 py-3">Lineage</th>
                <th className="px-4 py-3">Primary disease</th>
              </tr>
            </thead>
            <tbody className="divide-y divide-slate-100">
              {lines.map((l) => (
                <tr key={l.ach} className="hover:bg-slate-50/70">
                  <td className="px-4 py-2.5 font-mono font-semibold text-slate-800">{l.ach}</td>
                  <td className="px-4 py-2.5 text-slate-700">{l.name}</td>
                  <td className="px-4 py-2.5 capitalize text-slate-500">{l.lineage}</td>
                  <td className="px-4 py-2.5 text-slate-500">{l.disease}</td>
                </tr>
              ))}
              {lines.length === 0 && <EmptyRow cols={4} />}
            </tbody>
          </table>
        ) : (
          <table className="w-full text-left text-sm">
            <thead className="bg-slate-50 border-b border-slate-200 text-xs font-semibold text-slate-600 uppercase">
              <tr>
                <th className="px-4 py-3">Symbol</th>
                <th className="px-4 py-3">Ensembl id</th>
                <th className="px-4 py-3">Name</th>
              </tr>
            </thead>
            <tbody className="divide-y divide-slate-100">
              {genes.map((g) => (
                <tr key={g.symbol} className="hover:bg-slate-50/70">
                  <td className="px-4 py-2.5 font-mono font-semibold text-slate-800">{g.symbol}</td>
                  <td className="px-4 py-2.5 font-mono text-slate-500">{g.ensg}</td>
                  <td className="px-4 py-2.5 text-slate-700">{g.name}</td>
                </tr>
              ))}
              {genes.length === 0 && <EmptyRow cols={3} />}
            </tbody>
          </table>
        )}
      </div>

      <p className="text-[11px] text-slate-400">
        {tab === 'lines' ? `${lines.length} cell line${lines.length === 1 ? '' : 's'}` : `${genes.length} gene${genes.length === 1 ? '' : 's'}`} shown.
      </p>
    </div>
  );
};

const EmptyRow: React.FC<{ cols: number }> = ({ cols }) => (
  <tr>
    <td colSpan={cols} className="px-4 py-8 text-center text-sm text-slate-400">
      No matches — try a different search.
    </td>
  </tr>
);
