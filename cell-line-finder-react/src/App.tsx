import { useRef, useState } from 'react';
import { Sidebar } from './components/layout/Sidebar';
import { QueryBuilder } from './components/views/QueryBuilder';
import { ResultsTable } from './components/views/ResultsTable';
import { ProfileView } from './components/views/ProfileView';
import { AboutPage } from './components/views/AboutPage';
import { GuidePage } from './components/views/GuidePage';
import { ReferencePage } from './components/views/ReferencePage';
// Assistant chatbot page removed from routing/navigation — kept for possible reuse.
// import { AssistantPanel } from './components/views/AssistantPanel';
import { fetchRanking } from './lib/api';
import { lookupGeneSuggestion } from './lib/assistant';
import { InspectTarget, QueryParams, RankedCellLine, ResultMeta, ViewType } from './types';

// Matches the scoring API's 422 body, e.g. {"detail":"Gene not found: 'p53'"}
const GENE_NOT_FOUND = /Gene not found: '([^']+)'/;

export default function App() {
  const [currentView, setCurrentView] = useState<ViewType>('query');
  const [queryParams, setQueryParams] = useState<QueryParams>({
    targets: ['BRAF'],
    exclusions: ['KRAS'],
    lineage: 'Any',
    topK: 5
  });
  const [hasSearched, setHasSearched] = useState<boolean>(false);
  const [inspectTarget, setInspectTarget] = useState<InspectTarget | null>(null);

  // Live query state
  const [results, setResults] = useState<RankedCellLine[]>([]);
  const [meta, setMeta] = useState<ResultMeta | null>(null);
  const [loading, setLoading] = useState<boolean>(false);
  const [error, setError] = useState<string | null>(null);

  const abortRef = useRef<AbortController | null>(null);

  const runQuery = async (params: QueryParams) => {
    const target = params.targets[0];
    if (!target) {
      setError('Select at least one target gene.');
      return;
    }

    abortRef.current?.abort();
    const controller = new AbortController();
    abortRef.current = controller;

    setLoading(true);
    setError(null);
    try {
      // 2+ targets -> joint ranking; 1 target + exclusion(s) -> selectivity; else single.
      const { results: ranked, meta: m } = await fetchRanking(params.targets, params.exclusions, controller.signal);
      setResults(ranked);
      setMeta(m);
    } catch (e) {
      if ((e as Error).name === 'AbortError') return;
      const message = e instanceof Error ? e.message : 'Failed to fetch ranking.';
      const notFound = GENE_NOT_FOUND.exec(message);
      if (notFound) {
        // Ask the agent whether the typed name is a known alias/synonym before
        // giving up — "Gene not found" alone isn't actionable for the user.
        const suggestion = await lookupGeneSuggestion(notFound[1]);
        setError(
          suggestion.found && suggestion.symbol
            ? `Is this the gene you're trying to look for — ${suggestion.symbol}${
                suggestion.fullName ? ` (${suggestion.fullName})` : ''
              }?`
            : 'This gene is not present in the gene pool.'
        );
      } else {
        setError(message);
      }
      setResults([]);
      setMeta(null);
    } finally {
      if (abortRef.current === controller) setLoading(false);
    }
  };

  const handleSearchSubmit = (newParams: QueryParams) => {
    setQueryParams(newParams);
    setHasSearched(true);
    setCurrentView('results');
    void runQuery(newParams);
  };

  const handleInspect = (t: InspectTarget) => {
    setInspectTarget(t);
    setCurrentView('profile');
  };

  const handleReset = () => {
    abortRef.current?.abort();
    setQueryParams({ targets: ['BRAF'], exclusions: [], lineage: 'Any', topK: 5 });
    setResults([]);
    setMeta(null);
    setError(null);
    setLoading(false);
    setHasSearched(false);
    setInspectTarget(null);
    setCurrentView('query');
  };

  // topK isn't a server param; slice for display.
  const visibleResults = results.slice(0, queryParams.topK);

  return (
    <div className="flex h-screen w-screen overflow-hidden bg-[#fdfbfc] text-slate-900">
      <Sidebar
        currentView={currentView}
        setView={setCurrentView}
        onReset={handleReset}
        hasResults={hasSearched}
      />

      <main className="flex-1 overflow-y-auto p-8">
        {currentView === 'query' && (
          <QueryBuilder
            initialParams={queryParams}
            onSubmit={handleSearchSubmit}
            onOpenReference={() => setCurrentView('reference')}
          />
        )}

        {currentView === 'results' && (
          <ResultsTable
            queryParams={queryParams}
            results={visibleResults}
            meta={meta}
            loading={loading}
            error={error}
            onRetry={() => void runQuery(queryParams)}
            onInspect={handleInspect}
          />
        )}

        {currentView === 'profile' && inspectTarget && (
          <ProfileView
            target={inspectTarget}
            onBack={() => setCurrentView('results')}
            onInspect={handleInspect}
          />
        )}

        {currentView === 'about' && <AboutPage />}

        {currentView === 'guide' && (
          <GuidePage onStart={() => setCurrentView('query')} />
        )}

        {currentView === 'reference' && <ReferencePage />}
      </main>
    </div>
  );
}
