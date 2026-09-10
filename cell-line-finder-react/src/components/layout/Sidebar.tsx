import React from 'react';
import { Search, BarChart3, Info, RotateCcw, Dna, BookOpen, BookMarked } from 'lucide-react';
import { ViewType } from '../../types';

interface SidebarProps {
  currentView: ViewType;
  setView: (view: ViewType) => void;
  onReset: () => void;
  hasResults: boolean;
}

export const Sidebar: React.FC<SidebarProps> = ({ currentView, setView, onReset, hasResults }) => {
  const navClass = (active: boolean, disabled = false) =>
    `w-full flex items-center space-x-3 px-3 py-2.5 rounded-lg text-sm font-medium transition ${
      disabled
        ? 'opacity-40 cursor-not-allowed text-plum-muted'
        : active
        ? 'bg-mulberry-600 text-white shadow-sm'
        : 'text-plum-muted hover:bg-plum-accent hover:text-white'
    }`;

  return (
    <aside className="w-64 bg-plum text-white flex flex-col justify-between p-4 shrink-0 border-r border-plum-border">
      <div className="space-y-7">
        {/* Header */}
        <div className="flex items-center space-x-3 px-2 pt-1">
          <div className="p-2 bg-mulberry-600 rounded-lg text-white">
            <Dna className="w-6 h-6" />
          </div>
          <div>
            {/* <p className="text-[10px] font-semibold tracking-[0.14em] text-mulberry-300 uppercase">AstraZeneca R&amp;D</p> */}
            <h1 className="font-bold text-base leading-tight -mt-0.5">Cell Line Finder</h1>
          </div>
        </div>

        {/* Navigation items */}
        <nav className="space-y-1">
          <button onClick={() => setView('query')} className={navClass(currentView === 'query')}>
            <Search className="w-4 h-4" />
            <span>Find Cell Lines</span>
          </button>

          <button
            onClick={() => setView('results')}
            disabled={!hasResults}
            className={navClass(currentView === 'results' || currentView === 'profile', !hasResults)}
          >
            <BarChart3 className="w-4 h-4" />
            <span>Results</span>
          </button>

          <button onClick={() => setView('reference')} className={navClass(currentView === 'reference')}>
            <BookMarked className="w-4 h-4" />
            <span>Reference</span>
          </button>

          <button onClick={() => setView('guide')} className={navClass(currentView === 'guide')}>
            <BookOpen className="w-4 h-4" />
            <span>How to use</span>
          </button>

          <button onClick={() => setView('about')} className={navClass(currentView === 'about')}>
            <Info className="w-4 h-4" />
            <span>About</span>
          </button>
        </nav>
      </div>

      {/* Bottom Controls */}
      <div className="border-t border-plum-border pt-4 space-y-3">
        <button
          onClick={onReset}
          className="w-full flex items-center justify-center space-x-2 px-3 py-2 rounded-lg text-sm font-medium text-white/80 bg-plum-accent hover:bg-plum-border hover:text-white transition"
        >
          <RotateCcw className="w-4 h-4" />
          <span>New Query</span>
        </button>
        <p className="text-[11px] leading-relaxed text-plum-muted px-1">
          Internal research tool. Rankings are model predictions from multi-omics evidence — confirm
          experimentally before use.
        </p>
      </div>
    </aside>
  );
};
