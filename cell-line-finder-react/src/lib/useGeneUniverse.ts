import { useEffect, useState } from 'react';
import { fetchGeneList } from './api';
import { GENE_UNIVERSE } from './mockData';

/** The real, full gene universe (~19K symbols + ENSG ids), fetched once and
 * cached for the session -- lets the query builder validate a typed symbol
 * against the actual backend list instead of only finding out on submit.
 * Falls back to the small GENE_UNIVERSE demo list while loading or if the
 * fetch fails, so the search box is never empty. */

let cache: { symbols: string[]; symbolSet: Set<string>; ensgSet: Set<string> } | null = null;
let inflight: Promise<typeof cache> | null = null;

function load() {
  if (cache) return Promise.resolve(cache);
  if (!inflight) {
    inflight = fetchGeneList()
      .then((resp) => {
        const symbols = resp.genes.map((g) => g.symbol).sort();
        cache = {
          symbols,
          symbolSet: new Set(symbols.map((s) => s.toUpperCase())),
          ensgSet: new Set(resp.genes.map((g) => g.ensg.toUpperCase())),
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

  return { symbols, isKnown, ready };
}
