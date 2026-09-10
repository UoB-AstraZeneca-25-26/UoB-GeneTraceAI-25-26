import { AssistantContext, ChatMessage, QueryParams, RankedCellLine, ResultMeta } from '../types';
import { AGENT_API_URL } from '../config';

/** Build the context snapshot handed to the agent. Summarises the current query
 *  and the top handful of results so the agent can explain the ranking. */
export function buildAssistantContext(
  queryParams: QueryParams,
  meta: ResultMeta | null,
  results: RankedCellLine[]
): AssistantContext {
  return {
    mode: meta?.mode ?? null,
    genes: meta?.genes ?? queryParams.targets,
    exclusions: queryParams.exclusions,
    lineage: queryParams.lineage,
    primaryGene: meta?.primaryGene ?? queryParams.targets[0] ?? null,
    total: meta?.total ?? null,
    formula: meta?.formula ?? null,
    topLines: results.slice(0, 8).map((r) => ({
      rank: r.rank,
      cellLine: r.cellLine,
      modelId: r.modelId,
      score:
        r.mode === 'selectivity'
          ? r.selectivity
          : r.mode === 'jointSelectivity'
          ? r.combinedScore
          : r.mode === 'multi'
          ? r.jointScore
          : r.score,
      tier: null,
    })),
  };
}

export type SendMessage = (
  messages: ChatMessage[],
  context: AssistantContext | null
) => Promise<string>;

/* ---------------------------------------------------------------------------
 * Wired to the GeneTraceAI agent (Claude Haiku on Bedrock) behind a Lambda
 * Function URL, via the shared AGENT_API_URL (see ../config.ts) — the dev
 * build routes through Vite's same-origin proxy since the Lambda's CORS
 * allowlist doesn't cover localhost; a production build calls it directly.
 * We POST a single { query } string and read back the { answer } field. The
 * current on-screen ranking is folded into the query so questions like "why
 * is the top hit #1?" resolve to concrete cell lines and genes.
 * ------------------------------------------------------------------------- */

const AGENT_URL = `${AGENT_API_URL}/v1/agent/query`;

export type GeneSuggestion = { found: boolean; symbol: string | null; fullName: string | null };

/** Resolve a gene name the scoring API rejected as unknown, via the agent's
 *  gene_alias_lookup (local 19,213-gene panel + Ensembl/HGNC fallback). Used
 *  to offer "did you mean <SYMBOL>?" instead of a bare 404. */
export async function lookupGeneSuggestion(query: string): Promise<GeneSuggestion> {
  try {
    const res = await fetch(AGENT_URL, {
      method: 'POST',
      headers: { 'Content-Type': 'application/json' },
      body: JSON.stringify({ query: `What are the aliases for ${query}` }),
    });
    if (!res.ok) return { found: false, symbol: null, fullName: null };
    const data = await res.json();
    return {
      found: Boolean(data?.data?.found),
      symbol: data?.data?.symbol ?? null,
      fullName: data?.data?.full_name ?? null,
    };
  } catch {
    return { found: false, symbol: null, fullName: null };
  }
}

export const sendMessage: SendMessage = async (messages, context) => {
  const lastUser = [...messages].reverse().find((m) => m.role === 'user');
  const question = (lastUser?.content ?? '').trim();
  if (!question) return 'Ask me anything about the ranking or the underlying data.';

  let query = question;
  if (context && context.topLines.length > 0) {
    const genes = context.genes.join(', ');
    const top = context.topLines
      .slice(0, 5)
      .map((l) => `#${l.rank} ${l.cellLine} (${l.modelId}) score=${l.score.toFixed(3)}`)
      .join('; ');
    query = `Current ${context.mode ?? 'ranking'} for ${genes}. Top lines: ${top}. Question: ${question}`;
  }

  const res = await fetch(AGENT_URL, {
    method: 'POST',
    headers: { 'Content-Type': 'application/json' },
    body: JSON.stringify({ query }),
  });
  if (!res.ok) {
    const detail = await res.text().catch(() => '');
    throw new Error(`Agent error: ${res.status}${detail ? ` — ${detail.slice(0, 200)}` : ''}`);
  }
  const data = await res.json();
  return (data.answer as string) ?? 'The agent returned no answer.';
};

export const SUGGESTED_PROMPTS = [
  'How was this ranking done?',
  'Why is the top hit ranked #1?',
  'What does the tier mean?',
  'Explain the score for me',
];
