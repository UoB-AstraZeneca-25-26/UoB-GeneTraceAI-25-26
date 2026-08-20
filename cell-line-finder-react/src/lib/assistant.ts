import { AssistantContext, ChatMessage, QueryParams, RankedCellLine, ResultMeta } from '../types';

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
      score: r.mode === 'selectivity' ? r.selectivity : r.mode === 'multi' ? r.jointScore : r.score,
      tier: r.mode === 'single' ? r.tier : null,
    })),
  };
}

export type SendMessage = (
  messages: ChatMessage[],
  context: AssistantContext | null
) => Promise<string>;

/* ---------------------------------------------------------------------------
 * INTEGRATION SEAM — point this at your own agent.
 *
 * Never call an LLM provider directly from the browser (it would expose your API
 * key). POST to your own backend / agent endpoint, which holds the key and can
 * run retrieval, tools, etc. Example:
 *
 *   export const sendMessage: SendMessage = async (messages, context) => {
 *     const res = await fetch('/api/agent', {
 *       method: 'POST',
 *       headers: { 'Content-Type': 'application/json' },
 *       body: JSON.stringify({ messages, context }),
 *     });
 *     if (!res.ok) throw new Error(`Agent error: ${res.status}`);
 *     const data = await res.json();
 *     return data.reply as string;
 *   };
 *
 * The default below is a local stub so the UI is fully usable with no backend.
 * It answers a few ranking questions from the context and otherwise explains how
 * to wire the real agent. Replace it when your agent is ready.
 * ------------------------------------------------------------------------- */

const MODE_BLURB: Record<string, (c: AssistantContext) => string> = {
  single: (c) =>
    `This is a single-gene ranking for ${c.primaryGene}. Each of the ${c.total ?? 'many'} cell lines was scored on how strongly it supports ${c.primaryGene} across the available evidence layers, then ranked by that score. Tier (HIGH / MEDIUM / LOW / CONTEXT) is a separate confidence band — it is not derived from the score.`,
  selectivity: (c) =>
    `This is a selectivity ranking. Lines are scored by ${c.formula ?? 'score_high × (1 − score_low)'} — rewarding a high ${c.primaryGene} score together with a low score for the excluded gene. So the top hits are lines where your target is high and the thing you want to avoid is low.`,
  multi: (c) =>
    `This is a joint multi-gene ranking for ${c.genes.join(', ')}. Each line's joint score combines its per-gene scores. Note: a line missing a score for one of the genes is ranked on the genes it does have — so a high rank doesn't always mean the line is high across all of them. Watch for the "partial coverage" flag.`,
};

export const sendMessage: SendMessage = async (messages, context) => {
  await new Promise((r) => setTimeout(r, 450)); // mimic network latency

  const last = (messages[messages.length - 1]?.content ?? '').toLowerCase();
  const asksRanking = /(rank|score|why|how|selectiv|tier|joint|work)/.test(last);

  if (context && context.topLines.length > 0 && asksRanking) {
    const top = context.topLines[0];
    const blurb = MODE_BLURB[context.mode ?? 'single']?.(context) ?? '';
    const topLine = `Right now the top hit is ${top.cellLine} (${top.modelId})${
      top.tier ? `, tier ${top.tier}` : ''
    }, with a score of ${top.score.toFixed(4)}.`;
    return `${blurb}\n\n${topLine}\n\n(This is the built-in demo assistant. Wire your own agent in src/lib/assistant.ts to get real answers and follow-up reasoning.)`;
  }

  return "I'm the built-in demo assistant — I can sketch how the current ranking works, but I'm not wired to a real model yet. Point `sendMessage` in src/lib/assistant.ts at your agent endpoint to enable full answers. Meanwhile, try asking \u201chow was this ranking done?\u201d";
};

export const SUGGESTED_PROMPTS = [
  'How was this ranking done?',
  'Why is the top hit ranked #1?',
  'What does the tier mean?',
  'Explain the score for me',
];