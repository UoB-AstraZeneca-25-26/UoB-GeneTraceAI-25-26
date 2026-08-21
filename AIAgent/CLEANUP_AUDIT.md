# CellLineFinder Agent — Cleanup Audit

Generated as Step 1 of the cleanup pass. **Nothing has been removed yet** — this
is a findings report only. File:line references are against the current
working tree (several files already have uncommitted changes from an
in-progress Groq→Bedrock / local-first→live-first migration — `git status`
shows `AgentDevelopment/llm.py`, `api/settings.py`, `api/routes.py`,
`FunctionCalling.py`, `gene_index.py`, several `Tests/*.py`, `agent_details.md`,
`AIAgent-testResults.md` as modified, and `plan.md` as already deleted).

---

## A. Unused LLM provider code (`llm.py`)

| # | Finding | Location |
|---|---------|----------|
| A1 | `groq` branch in `get_llm()` imports `langchain_groq`, which **is no longer declared in either dependency file** — it was removed from `requirement_agents.txt` and `requirements.lock.txt` in the uncommitted diff (`langchain-groq==1.1.3` deleted). A fresh `pip install -r requirements.lock.txt` (i.e. the Docker image) will **not** have `langchain_groq` installed. `LLM_PROVIDER=groq` will now raise `ModuleNotFoundError` in any freshly built container. It only "works" on this dev machine because `langchain_groq` is a leftover global install. **This is a live bug in the rollback path, not just dead code.** | `AgentDevelopment/llm.py:37-48` |
| A2 | `get_llm()`'s default provider is still `"groq"` (`os.getenv("LLM_PROVIDER", "groq")`), and `Settings.llm_provider` defaults to `"groq"` too (`api/settings.py:36`). Both defaults are masked today because `.env` explicitly sets `LLM_PROVIDER=bedrock`, but an env without that var falls through to the now-broken Groq path (see A1). | `AgentDevelopment/llm.py:19`, `api/settings.py:36` |
| A3 | Default model fallback `"qwen-qwq-32b"` (`llm.py:20`) doesn't match the comment two lines below it, which talks about `"openai/gpt-oss-120b"` as "the default model" (lines 25-30). These two have been inconsistent with each other since before the Bedrock migration — not itself a Groq/Bedrock issue, but it's misleading either way and should be reconciled or dropped since `.env` always sets `LLM_MODEL` explicitly. | `AgentDevelopment/llm.py:20-30` |
| A4 | The raw `groq` SDK (not `langchain-groq`) is still declared in `requirements.lock.txt:34` and imported directly in `api/routes.py:21` for `GroqRateLimitError`. This import still resolves (package present), but the `except GroqRateLimitError` handlers (`routes.py:213`, `routes.py:242`) are only reachable if a `ChatGroq` instance is actually running — which A1 shows is no longer possible in a fresh install. So this is defensive code around a currently-unreachable code path. |  `api/routes.py:21,213,242` |
| A5 | `api/errors.py:42` docstring: `"""Raised when the LLM provider (Groq) 429s..."""` — `UpstreamRateLimitedError` is now raised for **both** Groq 429s and Bedrock `ThrottlingException` (`routes.py:216-217,245-246`), so the docstring undersells what the class now covers. |
| A6 | `huggingface`/`openai` stub branches (`llm.py:69-73`) raise `NotImplementedError` — harmless placeholders, not migration debris. No action needed. |

**Decision needed before touching any of A1-A5:** `Endpoint/test_results_report.md` §5 ("Provider Comparison") explicitly keeps Groq as a documented, tested alternative to Bedrock — this looks like a deliberate architectural choice (env-var-swappable LLM backend), not leftover cruft. See "Recommendation" at the bottom of this report.

---

## B. Stale environment variables

| # | Finding | Location |
|---|---------|----------|
| B1 | `.env` lines 1-16 are a **commented-out duplicate** of the same AWS/Bedrock config that's active in lines 18-23 below it — dead weight, and notably the commented block has a plaintext `AWS_ACCESS_KEY_ID` / `AWS_SECRET_ACCESS_KEY` pair sitting in it. Even commented out and gitignored, this is worth deleting rather than leaving around. | `.env:1-16` |
| B2 | `.env.example` was deleted (uncommitted, `D .env.example`) but it only ever contained a single blank line (`git diff` shows `-` with no content) — deleting it was correct, nothing to restore. No action needed, just confirming it's not an accidental loss. |
| B3 | `groq_api_key: SecretStr = SecretStr("")` in `Settings` (`api/settings.py:35`) — already made optional during the Bedrock migration. Whether to keep depends on the same Groq-rollback decision as category A. |
| B4 | `AWS_BEARER_TOKEN_BEDROCK` only appears in a documentation comment (`api/settings.py:44`) listing boto3's standard credential chain, and in `Tests/test_agent.py:17` as one of several credential checks for skipping live tests. Both are accurate, harmless mentions of a real boto3 credential mechanism — **not stale**, no action needed. |

---

## C. Dead code paths in `FunctionCalling.py`

| # | Finding | Location |
|---|---------|----------|
| C1 | **Contrary to the task brief's assumption**, `SCORING_API_URL` is **not** set anywhere in this repo (`.env` doesn't define it, nothing hardcodes it). `score_explainer` is still genuinely running in methodology-only mode as its *default*, not as a rare degradation path. The comment block at `FunctionCalling.py:220-225` ("Until that endpoint is wired up, the tool runs in methodology-only mode") is accurate as written and should **not** be changed to "fallback" framing yet — doing so would misdescribe current behaviour. |
| C2 | **Query parameter mismatch (flagged as a bug per the task brief, not verified against a live spec since no scoring API is reachable from this repo):** `score_explainer` builds `f"{scoring_url}?gene={ensg_id}&model={model_id}"` (`FunctionCalling.py:350`), and this exact contract is asserted by two tests (`Tests/test_score_explainer.py:147`, `:216-217`). If the real deployed scoring API actually expects `?query={gene}&cell_line={model_id}` (per the task brief), every API-mode call will silently 404/empty and fall back to methodology mode without any visible error. **This needs confirmation from whoever owns the scoring API spec before changing** — flagging only, not fixing, since I have no way to verify the real contract from this codebase. |
| C3 | No `_MOCK_SCORES` or hardcoded numeric scores anywhere in `AgentDevelopment/` — confirmed clean, and actively guarded by `Tests/test_score_explainer.py::test_no_hardcoded_scores_in_module`. No action needed. |

---

## D. Stale comments referencing outdated architecture

| # | Finding | Location |
|---|---------|----------|
| D1 | `gene_index.py` and `FunctionCalling.py` docstrings/comments already say "live-first, local fallback" — this part of the migration is **done** (uncommitted diff already rewrote them). No action needed. |
| D2 | Comments/docs reference `plan.md`'s gap numbers (`G1`, `G4`, `R7`) as if the file still exists, but `plan.md` was deleted (uncommitted `D plan.md`). Dangling references: `AgentDevelopment/FunctionCalling.py:20` ("gap G1"), `Tests/test_gene_alias.py:9` ("plan.md R7", "gap G4"), `AIAgent-testResults.md:13` ("`Endpoint/plan.md`'s gap analysis (G4/R7)"), `AIAgent-testResults.md:33` ("`Knowledge/` casing (G1) is fixed"), `AIAgent-testResults.md:35` ("Full gap list... live in `Endpoint/plan.md`" — that path never existed; `plan.md` lived at `AIAgent/plan.md`). These should keep the *why* (case-sensitivity bug, live-suite gap) but drop the pointer to a file that no longer exists. |
| D3 | `api/app.py:52` — "Prewarm: pay the **Groq** TLS handshake now" — stale, should be provider-neutral (the prewarm call goes through whichever `get_llm()` returns). |
| D4 | `api/security.py:6` — "so a runaway client can't silently burn the whole **Groq** quota" — stale, should be provider-neutral (Bedrock is billed per-token too, but the specific vendor name is now wrong/incomplete). |
| D5 | `Tests/conftest.py:7` and `Tests/test_latency.py:4,46,64` — comments/strings mentioning "Groq" in illustrative contexts (console encoding, a fake slow-LLM test double). Low priority, cosmetic; not incorrect, just vendor-specific phrasing for what's now a generic concern. |
| D6 | No remaining "0% success" / "unreliable" language about Ensembl anywhere — already cleaned up in the uncommitted diff (old comment removed from `FunctionCalling.py`). No action needed. |

---

## E. Outdated documentation files

### `agent_details.md`
- Architecture diagram still labels the LLM box "gpt-oss-120b (Groq)" (line 23) — no Bedrock/Claude Haiku 4.5 mention anywhere in this file.
- Line 319: "9 require live Ensembl/HGNC/**Groq** connections" — should say Bedrock.
- Line 337 (tech stack table): "LLM | OpenAI gpt-oss-120b via Groq | 120B MoE model, strong tool-calling, free tier sufficient for demo" — entirely describes the retired setup.
- Line 343: "The production code supports Groq; stubs exist for OpenAI and HuggingFace. Switching to a self-hosted model (Ollama on AWS EC2) requires changing one line in `.env`." — Bedrock isn't mentioned at all despite being the actual production path; Ollama was never implemented.
- Line 350 (Known gaps table): "Groq free tier rate limit... Mitigated by max_tokens=512" — both the provider and the max_tokens value (now 2048, per `llm.py:31`) are stale.
- The gene-resolution section (lines 91, 339-340) **has already been updated** to live-first language in the uncommitted diff — good, no action needed there.

### `AIAgent-testResults.md`
- Lines 1-6 (the very top, before the "Production readiness update" section): "**Model:** `openai/gpt-oss-120b` via Groq Cloud (free tier)" / architecture diagram at lines 43-81 showing "ChatGroq (gpt-oss-120b) ↔ ChatOllama (swap)" — this is the **original** pre-migration snapshot, and the file already has an honest disclaimer immediately below it ("Production readiness update (2026-08-17)... What changed... Gene resolution is live-API-first...") explaining that Tests 1-4 are historical. This is a deliberate, already-annotated "before" snapshot, not an oversight — recommend leaving the Test 1-4 walkthrough as-is (it's clearly marked as legacy) but see below for two items that should still change.
- Line 87: "LLM backend is swappable via `.env` (Groq → Ollama → OpenAI)" — Ollama/OpenAI were never implemented (`llm.py` only has Groq/Bedrock working + two `NotImplementedError` stubs); should say Bedrock, not Ollama.
- Line 258: "Tool calling latency: ~10–20s per query on Groq free tier (production deployment on Ollama/SageMaker would reduce this)" — stale; Bedrock is the actual production deployment now and its latency is already measured (`Endpoint/test_results_report.md` §5: "~2.3s (from playground)").

### `Endpoint/test_results_report.md` (dated 2026-08-20, **not** git-tracked — added to `.gitignore` in the uncommitted diff)
- This is the accurate, current report: correctly documents Bedrock/Claude Haiku 4.5 as the live provider, live-first gene resolution as primary, and explicitly keeps a "Provider Comparison" table (§5) comparing Groq vs. Bedrock — this is the strongest evidence in the repo that keeping the Groq code path is a **deliberate** decision, not an oversight.
- Minor inconsistency: Test 1's breakdown (line 36) shows `sources=[local_index]` for a plain `TP53` lookup — that's the old local-first default; under the now-live-first resolution order this specific historical run would be expected to show `sources=[ensembl, hgnc]` (or similar) unless it happened to hit the fallback. Likely just a stale capture from before the live-first change landed. Low priority since this file is gitignored/local-only.
- §3's category description for gene-alias tests ("Local-index hits (symbol and Ensembl-ID lookup)...") also undersells the now-primary network-first test coverage added in the uncommitted `Tests/test_gene_alias.py` diff.

### `agent_api.md` (not in the task's named list, but same issue)
- §"LLM-side (Groq free tier)" (line 622) and surrounding text (lines 558, 627, 631, 633) present Groq as the sole/current provider with no mention of Bedrock at all. Flagging as a bonus finding — same treatment as `agent_details.md` would apply here.

---

## F. Test suite

Already brought up to date in the uncommitted diff — confirmed by reading the
diffs directly, not just the task's assumptions:

- `Tests/test_agent.py` — now detects `LLM_PROVIDER` and checks the right
  credential env vars per provider (Bedrock or Groq) before skipping live
  tests; added `_extract_output_text()` to handle Bedrock's list-of-blocks
  `.content` shape. **No action needed.**
- `Tests/test_api.py` — added `test_valid_query_handles_bedrock_list_content_chunks`
  and `test_bedrock_throttling_returns_503_with_retry_after` **alongside**
  (not replacing) the existing Groq-rate-limit test. Consistent with keeping
  Groq as a documented, tested path. **No action needed** unless the Groq
  decision (see Recommendation) changes.
- `Tests/test_gene_alias.py` — already rewritten so `sources: ["ensembl", "hgnc"]`
  is the default-path assertion, and `sources: ["local_index"]` is now only
  asserted in the dedicated
  `test_gene_alias_lookup_falls_back_to_local_when_network_down` test. This
  is exactly what the task brief asked for — **already done, no action
  needed.**
- `Tests/test_tool.py` — already migrated off a direct `ChatGroq(...)`
  instantiation to `get_llm()`. **No action needed.**
- No duplicate or superseded test files found.
- Full offline run confirms all of the above: `python -m pytest Tests/ -m "not live" -q` → **50 passed, 9 deselected** (up from the "44 tests total" figure quoted in the stale docs above — another reason those counts need updating).

---

## G. Unused imports and dead functions

Ran `ruff check --select F401,F811,F841` (unused imports, redefinitions,
unused locals) across `AgentDevelopment/`, `api/`, `Schema/`, `main.py`,
`lambda_handler.py`: **all checks passed, nothing flagged.** No dead
functions found by cross-referencing definitions against call sites either.
This category is clean — no action needed.

---

## Step 5 summary (preliminary — nothing has been removed yet)

**Kept, and recommended to keep:** the Groq provider branch in `llm.py`,
`GroqRateLimitError` handling, and `groq_api_key` in `Settings`. This is a
deliberate architectural choice (LLM backend swappable via `.env`), and
`Endpoint/test_results_report.md` §5 documents Groq vs. Bedrock as a tested
comparison, not an accident. **Recommendation: don't delete the Groq
branch — fix it instead.** A1 shows the rollback path is currently broken
(missing dependency), which defeats the purpose of keeping it. Suggest:
restore `langchain-groq` to `requirement_agents.txt` /
`requirements.lock.txt` (small pin, cheap to keep valid), and flip the
*defaults* in `llm.py`/`settings.py` from `"groq"` to `"bedrock"` (A2) so a
misconfigured environment fails toward the currently-supported provider
instead of the broken one.

**Bug found, not yet fixed (needs owner confirmation):** C2, the
`score_explainer` query parameter names (`gene`/`model` vs. a possibly
different real contract `query`/`cell_line`). Not independently verifiable
from this repo since `SCORING_API_URL` is never actually configured here —
needs sign-off from whoever owns the scoring API before changing, since
getting it wrong silently degrades every API-mode call to methodology mode
with no visible error.

**Correction to the task's premise:** C1 — the real scoring API is **not**
actually connected in this codebase (`SCORING_API_URL` is unset
everywhere); methodology-only mode is still the genuine default, not a
degradation path. Comments describing it that way are accurate and were
left unchanged.

**Safe to clean up now, low risk, no dependency on the Groq decision:**
- D2 (dangling `plan.md` gap-number references)
- D3, D4 (Groq-specific wording in `app.py`/`security.py` comments)
- B1 (duplicate commented-out block in `.env`, including a plaintext AWS key)
- E: `agent_details.md`, `AIAgent-testResults.md` (and optionally
  `agent_api.md`) — update stale Groq/local-first/test-count claims, keep
  file structure and the already-good "here's what changed" framing intact.

**Not touched, and shouldn't be:** `reference/gene_lookup.parquet`
dependency, the local parquet fallback in `gene_alias_lookup`, any tool
definitions/agent wiring/routes, the frontend.

Next step: confirm the Groq rollback-vs-removal call and the C2 query-param
question with the user, then proceed to Step 3 (apply the low-risk D/B/E
fixes first, regardless of the Groq answer).
