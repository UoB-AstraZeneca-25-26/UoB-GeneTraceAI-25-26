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
import { CellLineDetailApiResponse, fetchRanking } from './lib/api';
import { lookupGeneSuggestion } from './lib/assistant';
import { GeneMatchSuggestion, InspectTarget, QueryParams, RankedCellLine, ResultMeta, ViewType } from './types';

// Matches the scoring API's 422 body, e.g. {"detail":"Gene not found: 'p53'"}
const GENE_NOT_FOUND = /Gene not found: '([^']+)'/;

export default function App() {
  const [currentView, setCurrentView] = useState<ViewType>('query');
  const [queryParams, setQueryParams] = useState<QueryParams>({
    targets: [],
    exclusions: [],
    lineage: 'Any',
    topK: 5
  });
  const [hasSearched, setHasSearched] = useState<boolean>(false);
  const [inspectTarget, setInspectTarget] = useState<InspectTarget | null>(null);

  // Back-navigation history for the profile view, populated only when a
  // profile is reached by clicking an RNA-similar alternative from another
  // profile (target.viaAlternative). A fresh inspect from ResultsTable clears
  // it -- 'back' from that normal path stays a hardcoded jump to 'results',
  // same as before this history existed.
  const [profileHistory, setProfileHistory] = useState<InspectTarget[]>([]);

  // Live query state
  const [results, setResults] = useState<RankedCellLine[]>([]);
  const [meta, setMeta] = useState<ResultMeta | null>(null);
  const [loading, setLoading] = useState<boolean>(false);
  const [error, setError] = useState<string | null>(null);
  const [suggestion, setSuggestion] = useState<GeneMatchSuggestion | null>(null);

  const abortRef = useRef<AbortController | null>(null);

  // Cache of /gene/detail responses keyed by "gene::modelId", shared between
  // ResultsTable's EvidenceMatrix and ProfileView so revisiting a line already
  // inspected in this session (or switching between the two views) doesn't
  // re-hit the API -- avoids redundant Lambda cold-starts and gives page
  // continuity when the user drills into the same line twice.
  const detailCacheRef = useRef<Map<string, CellLineDetailApiResponse>>(new Map());

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
    setSuggestion(null);
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
        const found = await lookupGeneSuggestion(notFound[1]);
        if (found.found && found.symbol) {
          // Kept as a structured suggestion (not baked into the error string)
          // so the UI can offer a one-click "search this instead" that swaps
          // the exact typed text for the resolved symbol and re-runs —
          // see handleAcceptSuggestion.
          setSuggestion({ originalQuery: notFound[1], symbol: found.symbol, fullName: found.fullName });
          setError(`"${notFound[1]}" isn't a recognized gene symbol.`);
        } else {
          setError('This gene is not present in the gene pool.');
        }
      } else {
        setError(message);
      }
      setResults([]);
      setMeta(null);
    } finally {
      if (abortRef.current === controller) setLoading(false);
    }
  };

  // Swaps the exact text the user typed (suggestion.originalQuery) for the
  // agent-resolved symbol wherever it appears — target or exclusion — then
  // re-runs immediately. Case-insensitive match since GeneSearch/manual entry
  // both uppercase on commit, but the failed API call's error text may not.
  const handleAcceptSuggestion = () => {
    if (!suggestion) return;
    const isMatch = (g: string) => g.trim().toUpperCase() === suggestion.originalQuery.trim().toUpperCase();
    const nextParams: QueryParams = {
      ...queryParams,
      targets: queryParams.targets.map((g) => (isMatch(g) ? suggestion.symbol : g)),
      exclusions: queryParams.exclusions.map((g) => (isMatch(g) ? suggestion.symbol : g)),
    };
    setQueryParams(nextParams);
    void runQuery(nextParams);
  };

  const handleSearchSubmit = (newParams: QueryParams) => {
    setQueryParams(newParams);
    setHasSearched(true);
    setCurrentView('results');
    setProfileHistory([]);
    void runQuery(newParams);
  };

  const handleInspect = (t: InspectTarget) => {
    // Reached via an RNA-alternative click from within a profile -- push the
    // profile we're leaving so "back" can return to it, instead of always
    // hard-jumping to 'results'. A normal inspect from ResultsTable (no
    // viaAlternative flag) starts a fresh chain: clear any stale history so
    // a later alternative-click in this new chain doesn't pop into a profile
    // left over from a previous line.
    if (t.viaAlternative && inspectTarget) {
      setProfileHistory((prev) => [...prev, inspectTarget]);
    } else {
      setProfileHistory([]);
    }
    setInspectTarget(t);
    setCurrentView('profile');
  };

  // 'Back' from a profile: pop the alternative-chain history if there is one,
  // otherwise fall back to 'results' -- preserves the original Path A
  // behavior exactly when no alternative was ever involved.
  const handleProfileBack = () => {
    if (profileHistory.length === 0) {
      setCurrentView('results');
      return;
    }
    const previous = profileHistory[profileHistory.length - 1];
    setProfileHistory(profileHistory.slice(0, -1));
    setInspectTarget(previous);
    setCurrentView('profile');
  };

  const handleReset = () => {
    abortRef.current?.abort();
    setQueryParams({ targets: [], exclusions: [], lineage: 'Any', topK: 5 });
    setResults([]);
    setMeta(null);
    setError(null);
    setSuggestion(null);
    setLoading(false);
    setHasSearched(false);
    setInspectTarget(null);
    setProfileHistory([]);
    setCurrentView('query');
    detailCacheRef.current.clear();
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
            value={queryParams}
            onChange={setQueryParams}
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
            suggestion={suggestion}
            onAcceptSuggestion={handleAcceptSuggestion}
            onRetry={() => void runQuery(queryParams)}
            onInspect={handleInspect}
            detailCache={detailCacheRef.current}
          />
        )}

        {currentView === 'profile' && inspectTarget && (
          <ProfileView
            target={inspectTarget}
            onBack={handleProfileBack}
            onInspect={handleInspect}
            detailCache={detailCacheRef.current}
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
