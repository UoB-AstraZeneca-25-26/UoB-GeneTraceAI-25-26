import { useRef, useState } from 'react';
import { Sidebar } from './components/layout/Sidebar';
import { QueryBuilder } from './components/views/QueryBuilder';
import { ResultsTable } from './components/views/ResultsTable';
import { ProfileView } from './components/views/ProfileView';
import { AboutPage } from './components/views/AboutPage';
import { MOCK_DATA } from './lib/mockData';
import { fetchRanking } from './lib/api';
import { QueryParams, RankedCellLine, ResultMeta, ViewType } from './types';

export default function App() {
  const [currentView, setCurrentView] = useState<ViewType>('query');
  const [queryParams, setQueryParams] = useState<QueryParams>({
    targets: ['BRAF'],
    exclusions: ['KRAS'],
    lineage: 'Any',
    topK: 5
  });
  const [hasSearched, setHasSearched] = useState<boolean>(false);
  const [selectedCellLine, setSelectedCellLine] = useState<string | null>(null);

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
      // One exclusion routes to the selectivity endpoint; none routes to /gene.
      const exclusion = params.exclusions[0];
      const { results: ranked, meta: m } = await fetchRanking(target, exclusion, controller.signal);
      setResults(ranked);
      setMeta(m);
    } catch (e) {
      if ((e as Error).name === 'AbortError') return;
      setError(e instanceof Error ? e.message : 'Failed to fetch ranking.');
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

  const handleReset = () => {
    abortRef.current?.abort();
    setQueryParams({ targets: ['BRAF'], exclusions: [], lineage: 'Any', topK: 5 });
    setResults([]);
    setMeta(null);
    setError(null);
    setLoading(false);
    setHasSearched(false);
    setSelectedCellLine(null);
    setCurrentView('query');
  };

  // topK isn't a server param; slice for display.
  const visibleResults = results.slice(0, queryParams.topK);

  return (
    <div className="flex h-screen w-screen overflow-hidden bg-slate-50 text-slate-900">
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
            onSelectCellLine={(cl) => {
              setSelectedCellLine(cl);
              setCurrentView('profile');
            }}
          />
        )}

        {currentView === 'profile' && selectedCellLine && (
          <ProfileView
            cellLine={selectedCellLine}
            dataset={MOCK_DATA}
            onBack={() => setCurrentView('results')}
          />
        )}

        {currentView === 'about' && <AboutPage />}
      </main>
    </div>
  );
}