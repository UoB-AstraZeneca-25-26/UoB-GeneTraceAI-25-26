import { RankedCellLine, ResultMeta } from '../types';

export interface ReportData {
  title: string;
  queryLine: string;
  totalsLine: string;
  generatedAt: string;
  columns: string[];
  rows: (string | number)[][];
}

const fmt4 = (n: number) => n.toFixed(4);
const fmt3 = (n: number) => n.toFixed(3);
const geneCol = (v: number | null) => (v == null ? '—' : fmt3(v));

// Cell-line context columns, shared across every mode -- inserted right after
// Lineage so the researcher-facing "is this a sensible model" info sits next
// to the tissue it's already grouped by, before the mode-specific scores.
const CONTEXT_COLUMNS = ['Primary Disease', 'Subtype', 'Growth Pattern', 'Sex', 'RRID'];
const contextRow = (r: RankedCellLine): (string | number)[] => [
  r.primaryDisease ?? '—',
  r.subtype ?? '—',
  r.growthPattern ?? '—',
  r.sex ?? '—',
  r.rrid ? r.rrid.toUpperCase() : '—',
];

/** Normalizes any of the four result modes into one flat table shape --
 *  same column/row logic ResultsTable.tsx already renders on screen, just
 *  serialized instead of laid out as JSX. CSV and PDF export both build off
 *  this single source of truth so the two formats never drift apart. */
export function buildReportData(results: RankedCellLine[], meta: ResultMeta | null): ReportData {
  const mode = results[0]?.mode ?? meta?.mode ?? 'single';

  const queryLine =
    mode === 'jointSelectivity'
      ? `Joint ranking: ${meta?.genes?.join(' + ') ?? ''} — selective against ${meta?.excludedGenes?.join(', ') ?? ''}`
      : mode === 'multi'
      ? `Joint ranking: ${meta?.genes?.join(' + ') ?? ''}`
      : mode === 'selectivity'
      ? `Selective for ${meta?.primaryGene ?? ''} against ${meta?.excludedGenes?.join(', ') ?? ''}`
      : `Target: ${meta?.primaryGene ?? ''}`;

  const totalsLine = meta
    ? `Showing ${results.length} of ${meta.total.toLocaleString()} ${mode === 'multi' || mode === 'jointSelectivity' ? 'passing' : 'ranked'}`
    : `${results.length} lines`;

  let columns: string[];
  let rows: (string | number)[][];

  if (mode === 'multi') {
    const genes = (results as Extract<RankedCellLine, { mode: 'multi' }>[])[0]?.genes ?? meta?.genes ?? [];
    columns = ['Rank', 'Model ID', 'Cell Line', 'Lineage', ...CONTEXT_COLUMNS, ...genes, 'Limiting Gene', 'Joint Score', 'In-Lineage Score'];
    rows = (results as Extract<RankedCellLine, { mode: 'multi' }>[]).map((r) => [
      r.rank,
      r.modelId,
      r.cellLine,
      r.lineage,
      ...contextRow(r),
      ...genes.map((g) => geneCol(r.geneScores[g])),
      r.limitingGene ?? '—',
      fmt4(r.jointScore),
      fmt3(r.lineageScore),
    ]);
  } else if (mode === 'selectivity') {
    const rowsTyped = results as Extract<RankedCellLine, { mode: 'selectivity' }>[];
    const geneHigh = rowsTyped[0]?.geneHigh ?? meta?.primaryGene ?? 'Target';
    const excluded = rowsTyped[0]?.excludedGenes ?? meta?.excludedGenes ?? [];
    columns = ['Rank', 'Model ID', 'Cell Line', 'Lineage', ...CONTEXT_COLUMNS, geneHigh, ...excluded, 'Selectivity', 'In-Lineage Score'];
    rows = rowsTyped.map((r) => [
      r.rank,
      r.modelId,
      r.cellLine,
      r.lineage,
      ...contextRow(r),
      fmt3(r.scoreHigh),
      ...excluded.map((g) => fmt3(r.exclusionScores[g] ?? 0)),
      fmt4(r.selectivity),
      fmt3(r.lineageScore),
    ]);
  } else if (mode === 'jointSelectivity') {
    const rowsTyped = results as Extract<RankedCellLine, { mode: 'jointSelectivity' }>[];
    const genes = rowsTyped[0]?.genes ?? meta?.genes ?? [];
    const excluded = rowsTyped[0]?.excludedGenes ?? meta?.excludedGenes ?? [];
    columns = [
      'Rank', 'Model ID', 'Cell Line', 'Lineage', ...CONTEXT_COLUMNS,
      ...genes, 'Limiting Gene', ...excluded, 'Combined Score', 'In-Lineage Score',
    ];
    rows = rowsTyped.map((r) => [
      r.rank,
      r.modelId,
      r.cellLine,
      r.lineage,
      ...contextRow(r),
      ...genes.map((g) => geneCol(r.geneScores[g])),
      r.limitingGene ?? '—',
      ...excluded.map((g) => fmt3(r.exclusionScores[g] ?? 0)),
      fmt4(r.combinedScore),
      fmt3(r.lineageScore),
    ]);
  } else {
    const rowsTyped = results as Extract<RankedCellLine, { mode: 'single' }>[];
    columns = ['Rank', 'Model ID', 'Cell Line', 'Lineage', ...CONTEXT_COLUMNS, 'Score', 'In-Lineage Score'];
    rows = rowsTyped.map((r) => [r.rank, r.modelId, r.cellLine, r.lineage, ...contextRow(r), fmt4(r.score), fmt3(r.lineageScore)]);
  }

  return {
    title: 'Cell Line Finder — Ranked Results',
    queryLine,
    totalsLine,
    generatedAt: new Date().toLocaleString(),
    columns,
    rows,
  };
}

function csvEscape(v: string | number): string {
  const s = String(v);
  return /[",\n]/.test(s) ? `"${s.replace(/"/g, '""')}"` : s;
}

function downloadBlob(blob: Blob, filename: string) {
  const url = URL.createObjectURL(blob);
  const a = document.createElement('a');
  a.href = url;
  a.download = filename;
  document.body.appendChild(a);
  a.click();
  document.body.removeChild(a);
  URL.revokeObjectURL(url);
}

const reportFilename = (ext: string) => `cell-line-finder-results-${new Date().toISOString().slice(0, 10)}.${ext}`;

export function exportReportCsv(results: RankedCellLine[], meta: ResultMeta | null) {
  const data = buildReportData(results, meta);
  const lines = [
    `# ${data.title}`,
    `# ${data.queryLine}`,
    `# ${data.totalsLine}`,
    `# Generated: ${data.generatedAt}`,
    data.columns.map(csvEscape).join(','),
    ...data.rows.map((row) => row.map(csvEscape).join(',')),
  ];
  downloadBlob(new Blob([lines.join('\n')], { type: 'text/csv;charset=utf-8' }), reportFilename('csv'));
}

export async function exportReportPdf(results: RankedCellLine[], meta: ResultMeta | null) {
  const [{ default: jsPDF }, autoTable] = await Promise.all([
    import('jspdf'),
    import('jspdf-autotable').then((m) => m.default),
  ]);
  const data = buildReportData(results, meta);

  const doc = new jsPDF({ orientation: data.columns.length > 7 ? 'landscape' : 'portrait' });
  doc.setFontSize(14);
  doc.text(data.title, 14, 16);
  doc.setFontSize(9);
  doc.setTextColor(90);
  doc.text(data.queryLine, 14, 23);
  doc.text(`${data.totalsLine}  ·  Generated ${data.generatedAt}`, 14, 28);

  autoTable(doc, {
    startY: 34,
    head: [data.columns],
    body: data.rows,
    styles: { fontSize: 8, cellPadding: 2 },
    headStyles: { fillColor: [130, 20, 75] }, // matches the app's mulberry accent
    margin: { left: 14, right: 14 },
  });

  doc.save(reportFilename('pdf'));
}
