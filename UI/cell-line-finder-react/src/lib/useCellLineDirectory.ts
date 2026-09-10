import { useEffect, useState } from 'react';
import { fetchCellLineList } from './api';
import { CELL_LINE_DIRECTORY } from './mockApi';

/** The real, full cell-line directory (~1.8K lines), fetched once and cached
 * for the session -- powers the Reference lookup page's "Cell lines" tab.
 * Falls back to the small CELL_LINE_DIRECTORY demo list while loading or if
 * the fetch fails, so the table is never empty. */

type CellLineRow = { ach: string; name: string; lineage: string; disease: string };

let cache: CellLineRow[] | null = null;
let inflight: Promise<CellLineRow[] | null> | null = null;

function load() {
  if (cache) return Promise.resolve(cache);
  if (!inflight) {
    inflight = fetchCellLineList()
      .then((resp) => {
        const rows = resp.cell_lines.map((c) => ({
          ach: c.model_id,
          name: c.name ? c.name.toUpperCase() : c.model_id,
          lineage: (c.lineage ?? 'unknown').replace(/_/g, ' '),
          disease: c.primary_disease ?? 'unknown',
        }));
        cache = rows;
        return rows;
      })
      .catch(() => null); // enrichment only -- CELL_LINE_DIRECTORY fallback covers this
  }
  return inflight;
}

export function useCellLineDirectory() {
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

  return { lines: cache ?? CELL_LINE_DIRECTORY, ready };
}
