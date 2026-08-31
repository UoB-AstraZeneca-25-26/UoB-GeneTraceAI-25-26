import { useEffect, useState } from 'react';
import { fetchGeneList } from './api';
import { GENE_UNIVERSE } from './mockData';

/** The real, full gene universe (~19K symbols + ENSG ids), fetched once and
 * cached for the session -- lets the query builder validate a typed symbol
 * against the actual backend list instead of only finding out on submit.
 * Falls back to the small GENE_UNIVERSE demo list while loading or if the
 * fetch fails, so the search box is never empty. */

type GeneRow = { symbol: string; ensg: string; name: string | null };
type GeneCache = { symbols: string[]; symbolSet: Set<string>; ensgSet: Set<string>; rows: GeneRow[] };

let cache: GeneCache | null = null;
let inflight: Promise<GeneCache | null> | null = null;

function load() {
  if (cache) return Promise.resolve(cache);
  if (!inflight) {
    inflight = fetchGeneList()
      .then((resp) => {
        const rows = [...resp.genes].sort((a, b) => a.symbol.localeCompare(b.symbol));
        const symbols = rows.map((g) => g.symbol);
        cache = {
          symbols,
          symbolSet: new Set(symbols.map((s) => s.toUpperCase())),
          ensgSet: new Set(rows.map((g) => g.ensg.toUpperCase())),
          rows,
        };
        return cache;
      })
      .catch(() => null); // enrichment only -- GENE_UNIVERSE fallback covers this
  }
  return inflight;
}

export function useGeneUniverse() {
  const [ready, setReady] = useState(cache !== null);

  useEffect(() => {
    if (cache) return;
    let cancelled = false;
    load().then(() => {
      if (!cancelled) setReady(true);
    });
    return () => {
      cancelled = true;
    };
  }, []);

  const symbols = cache?.symbols ?? GENE_UNIVERSE;

  // Undefined while the real list hasn't loaded yet -- callers should treat
  // "unknown" as "don't judge it" rather than flagging a false positive.
  const isKnown = cache
    ? (query: string) => {
        const q = query.trim().toUpperCase();
        return cache!.symbolSet.has(q) || cache!.ensgSet.has(q);
      }
    : undefined;

  return { symbols, isKnown, ready, rows: cache?.rows ?? null };
}
