import React, { useState, useRef, useEffect, useMemo } from 'react';
import { Search } from 'lucide-react';

interface GeneSearchProps {
  suggestions: string[];        // known symbols for type-ahead
  exclude?: string[];           // symbols to hide (already chosen)
  placeholder?: string;
  onPick: (symbol: string) => void;
}

/* Type-ahead gene picker with keyboard nav. Unlike the fixed button grid, this
   accepts ANY symbol the user types — the live API isn't limited to the demo
   gene list, so neither should the input be. */
export const GeneSearch: React.FC<GeneSearchProps> = ({ suggestions, exclude = [], placeholder, onPick }) => {
  const [value, setValue] = useState('');
  const [open, setOpen] = useState(false);
  const [hi, setHi] = useState(0);
  const boxRef = useRef<HTMLDivElement>(null);

  const matches = useMemo(() => {
    const q = value.trim().toUpperCase();
    return suggestions
      .filter((s) => !exclude.includes(s))
      .filter((s) => (q ? s.toUpperCase().includes(q) : true))
      .slice(0, 8);
  }, [value, suggestions, exclude]);

  // When the typed text isn't an exact known symbol, surface an explicit
  // "Add <value>" option so users can see that a custom symbol is committable
  // (not just discover it by blindly pressing Enter).
  const q = value.trim().toUpperCase();
  const showCustomAdd = q.length > 0 && !exclude.includes(q) && !matches.includes(q);
  const options = showCustomAdd ? [...matches, q] : matches;

  useEffect(() => {
    const h = (e: MouseEvent) => {
      if (boxRef.current && !boxRef.current.contains(e.target as Node)) setOpen(false);
    };
    document.addEventListener('mousedown', h);
    return () => document.removeEventListener('mousedown', h);
  }, []);

  const commit = (raw: string) => {
    const sym = raw.trim().toUpperCase();
    if (!sym) return;
    onPick(sym);
    setValue('');
    setOpen(false);
    setHi(0);
  };

  return (
    <div className="relative" ref={boxRef}>
      <Search className="w-4 h-4 absolute left-3 top-1/2 -translate-y-1/2 text-slate-400 pointer-events-none" />
      <input
        value={value}
        placeholder={placeholder || 'Type any gene symbol…'}
        spellCheck={false}
        onFocus={() => setOpen(true)}
        onChange={(e) => {
          setValue(e.target.value);
          setOpen(true);
          setHi(0);
        }}
        onKeyDown={(e) => {
          if (e.key === 'ArrowDown') {
            e.preventDefault();
            setHi((h) => Math.min(h + 1, options.length - 1));
          } else if (e.key === 'ArrowUp') {
            e.preventDefault();
            setHi((h) => Math.max(h - 1, 0));
          } else if (e.key === 'Enter') {
            e.preventDefault();
            commit(options[hi] ?? value);
          } else if (e.key === 'Escape') {
            setOpen(false);
          }
        }}
        className="w-full bg-slate-50 border border-slate-300 rounded-lg pl-9 pr-3 py-2 text-sm font-mono focus:ring-2 focus:ring-indigo-500 focus:outline-none"
      />
      {open && options.length > 0 && (
        <div className="absolute z-30 top-[calc(100%+4px)] left-0 right-0 bg-white border border-slate-200 rounded-lg shadow-lg p-1 max-h-64 overflow-auto">
          {options.map((m, i) => {
            const isCustom = showCustomAdd && i === options.length - 1;
            return (
              <button
                type="button"
                key={m}
                onMouseEnter={() => setHi(i)}
                onClick={() => commit(m)}
                className={`w-full text-left px-3 py-2 rounded-md text-sm font-mono ${
                  i === hi ? 'bg-indigo-50 text-indigo-700' : 'text-slate-700 hover:bg-slate-50'
                }`}
              >
                {isCustom ? (
                  <span>
                    Add <span className="font-semibold">"{m}"</span>
                  </span>
                ) : (
                  m
                )}
              </button>
            );
          })}
        </div>
      )}
    </div>
  );
};
