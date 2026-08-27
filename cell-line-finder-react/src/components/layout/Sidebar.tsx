import React from 'react';
import { Search, BarChart3, Info, RotateCcw, Dna, BookOpen } from 'lucide-react';
import { ViewType } from '../../types';

interface SidebarProps {
  currentView: ViewType;
  setView: (view: ViewType) => void;
  onReset: () => void;
  hasResults: boolean;
}

export const Sidebar: React.FC<SidebarProps> = ({ currentView, setView, onReset, hasResults }) => {
  return (
    <aside className="w-64 bg-slate-900 text-slate-100 flex flex-col justify-between p-4 shrink-0 border-r border-slate-800">
      <div className="space-y-6">
        {/* Header */}
        <div className="flex items-center space-x-3 px-2">
          <div className="p-2 bg-indigo-600 rounded-lg text-white">
            <Dna className="w-6 h-6" />
          </div>
          <div>
            <h1 className="font-bold text-base leading-tight">Cell Line Finder</h1>
            <p className="text-xs text-slate-400">Multi-Omics Selection</p>
          </div>
        </div>

        {/* Navigation items */}
        <nav className="space-y-1">
          <button
            onClick={() => setView('query')}
            className={`w-full flex items-center space-x-3 px-3 py-2.5 rounded-lg text-sm font-medium transition ${
              currentView === 'query' ? 'bg-indigo-600 text-white' : 'text-slate-400 hover:bg-slate-800 hover:text-slate-200'
            }`}
          >
            <Search className="w-4 h-4" />
            <span>Find Cell Lines</span>
          </button>

          <button
            onClick={() => setView('results')}
            disabled={!hasResults}
            className={`w-full flex items-center space-x-3 px-3 py-2.5 rounded-lg text-sm font-medium transition ${
              !hasResults ? 'opacity-40 cursor-not-allowed text-slate-500' :
              currentView === 'results' || currentView === 'profile' ? 'bg-indigo-600 text-white' : 'text-slate-400 hover:bg-slate-800 hover:text-slate-200'
            }`}
          >
            <BarChart3 className="w-4 h-4" />
            <span>Results</span>
          </button>

          <button
            onClick={() => setView('guide')}
            className={`w-full flex items-center space-x-3 px-3 py-2.5 rounded-lg text-sm font-medium transition ${
              currentView === 'guide' ? 'bg-indigo-600 text-white' : 'text-slate-400 hover:bg-slate-800 hover:text-slate-200'
            }`}
          >
            <BookOpen className="w-4 h-4" />
            <span>How to use</span>
          </button>

          <button
            onClick={() => setView('about')}
            className={`w-full flex items-center space-x-3 px-3 py-2.5 rounded-lg text-sm font-medium transition ${
              currentView === 'about' ? 'bg-indigo-600 text-white' : 'text-slate-400 hover:bg-slate-800 hover:text-slate-200'
            }`}
          >
            <Info className="w-4 h-4" />
            <span>About</span>
          </button>
        </nav>
      </div>

      {/* Bottom Controls */}
      <div className="border-t border-slate-800 pt-4 space-y-3">
        <button
          onClick={onReset}
          className="w-full flex items-center justify-center space-x-2 px-3 py-2 rounded-lg text-sm font-medium text-slate-400 bg-slate-800 hover:bg-slate-700 hover:text-white transition"
        >
          <RotateCcw className="w-4 h-4" />
          <span>New Query</span>
        </button>
      </div>
    </aside>
  );
};