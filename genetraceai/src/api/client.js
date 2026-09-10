/* ------------------------------------------------------------------ */
/*  API client — all external calls live here                          */
/*                                                                     */
/*  Scoring API: AWS Lambda (deployed)                                 */
/*  Agent API:   FastAPI on localhost:8000 (optional, may be offline)   */
/* ------------------------------------------------------------------ */

const SCORING_API =
  "https://9368clqa34.execute-api.eu-north-1.amazonaws.com";
const AGENT_API = "http://localhost:8000";

/* ---- Scoring & Ranking ---- */

/**
 * Search for a gene and get the top-30 ranked cell lines.
 * GET /gene?query=BRAF
 */
export async function searchGene(query) {
  const res = await fetch(
    `${SCORING_API}/gene?query=${encodeURIComponent(query)}`
  );
  if (!res.ok) throw new Error(`Gene search failed (${res.status})`);
  return res.json();
}

/**
 * Get the rank and score for one specific gene × cell-line pair,
 * plus RNA-similarity alternatives for that cell line.
 * GET /gene?query=BRAF&cell_line=ACH-000219
 */
export async function getGeneCellLine(query, cellLine) {
  const res = await fetch(
    `${SCORING_API}/gene?query=${encodeURIComponent(query)}&cell_line=${encodeURIComponent(cellLine)}`
  );
  if (!res.ok) throw new Error(`Cell line lookup failed (${res.status})`);
  return res.json();
}

/**
 * Get cell lines ranked by selectivity: high expression of gene_a
 * AND low expression of gene_b.
 * GET /exclude?gene_a=BRAF&gene_b=EGFR
 */
export async function getSelectivity(geneA, geneB) {
  const res = await fetch(
    `${SCORING_API}/exclude?gene_a=${encodeURIComponent(geneA)}&gene_b=${encodeURIComponent(geneB)}`
  );
  if (!res.ok) throw new Error(`Selectivity query failed (${res.status})`);
  return res.json();
}

/* ---- AI Agent ---- */

/**
 * Query the AI agent for alias lookups, score explanations,
 * or dataset descriptions.
 * POST http://localhost:8000/v1/agent/query
 */
export async function queryAgent(queryString) {
  const res = await fetch(`${AGENT_API}/v1/agent/query`, {
    method: "POST",
    headers: { "Content-Type": "application/json" },
    body: JSON.stringify({ query: queryString }),
  });
  if (!res.ok) throw new Error(`Agent query failed (${res.status})`);
  return res.json();
}
