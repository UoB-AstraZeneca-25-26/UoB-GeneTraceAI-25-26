import React, { useState, useMemo, useRef, useEffect } from "react";

/* ------------------------------------------------------------------ */
/*  Evidence tracks — the six layers every cell line is scored on     */
/* ------------------------------------------------------------------ */
/* Six tracks. `kind` records what the number actually is, because the views
   must not imply a magnitude where the pipeline only has a fact:
     level  -> within-gene percentile of a measured quantity (real gradient)
     ratio  -> deviation from diploid copy number (real gradient)
     call   -> categorical evidence; shown at full or half strength, never
               interpolated into a fake magnitude                              */
const TRACKS = [
  { key: "expression", label: "Expression", color: "#3E63DD", kind: "level" },
  { key: "proteomics", label: "Proteomics", color: "#30A46C", kind: "level" },
  { key: "cna",        label: "Copy number", color: "#0B7C86", kind: "ratio" },
  { key: "mutation",   label: "Mutation",   color: "#E5484D", kind: "call" },
  { key: "fusion",     label: "Fusion",     color: "#8E4EC6", kind: "call" },
  { key: "signature",  label: "Signature",  color: "#F5A623", kind: "call" },
];
const TRACK_BY_KEY = Object.fromEntries(TRACKS.map((t) => [t.key, t]));

/* One line of plain english for a track on one cell line. */
function trackFinding(key, d) {
  if (!d || !d.present) return "not measured";
  if (key === "expression") {
    const x = d.detail || {};
    const bits = [];
    if (x.depmap_log2tpm != null) bits.push(`DepMap log2TPM ${x.depmap_log2tpm}`);
    if (x.hpa_ntpm != null) bits.push(`HPA nTPM ${x.hpa_ntpm}`);
    return bits.length ? bits.join(" · ") : (d.sources || []).join(" + ");
  }
  if (key === "proteomics") {
    const x = d.detail || {};
    const bits = [];
    if (x.ccle != null) bits.push(`CCLE ${x.ccle}`);
    if (x.procan != null) bits.push(`ProCan ${x.procan}`);
    return bits.length ? bits.join(" · ") : (d.sources || []).join(" + ");
  }
  if (key === "cna") return `copy number ${d.total_cn}${d.direction ? " — " + d.direction : ""}`;
  if (key === "mutation") return d.driver_alteration
    ? "driver alteration present in this line"
    : "mutations measured, no driver call";
  if (key === "fusion") return "fusion calls available";
  if (key === "signature") return "mutational signatures available";
  return "";
}

const LIN_COLOR = (l) => ({
  Lung: "#3E63DD", Skin: "#E5484D", Colon: "#F5A623", Pancreas: "#8E4EC6",
  Breast: "#D6409F", Stomach: "#30A46C", Lymphoid: "#0B7C86",
}[l] || "#5D7183");

/* ------------------------------------------------------------------ */
/*  Data — grounded in real cell-line biology                          */
/* ------------------------------------------------------------------ */
/* ------------------------------------------------------------------ */
/*  Data layer — real pipeline output, fetched as static JSON          */
/*                                                                     */
/*  Written by:  python src/pipeline/export_web.py --genes BRAF ...    */
/*  Served from: public/data/<SYMBOL>.json  (Vite serves public/ at /) */
/*                                                                     */
/*  No backend: core_score.parquet is 29.8M rows, but one gene's       */
/*  top-N slice is ~48 KB, so each gene is a static file.              */
/* ------------------------------------------------------------------ */
const DATA_BASE = "/data";

/* Per-track values are REAL: expression and proteomics are within-gene
   percentiles of the measured quantity, copy number is deviation from
   diploid, and the three categorical layers carry a call rather than an
   invented magnitude. See `kind` on TRACKS above. */
function adaptTracks(t) {
  if (!t) return {};
  const out = {};
  for (const tr of TRACKS) {
    const d = t[tr.key] || {};
    out[tr.key] = { ...d, score: d.present ? (d.score ?? 0) : 0,
                    finding: trackFinding(tr.key, d) };
  }
  return out;
}

/* core_score is a WITHIN-GENE percentile in [0,1]. It is a real measured
   quantity and is what drives every bar width and node radius below.
   It is NOT a calibrated confidence — the ordinal tier carries that, and
   the pipeline deliberately exports no 0-1 confidence float (design
   record C1/C5). */
function adaptPayload(j) {
  return {
    symbol: j.symbol,
    ensg: j.ensg,
    name: j.gene?.name || "",
    role: j.gene?.["COSMIC role"] || "unknown",
    chr: j.gene?.["locus type"] || "",
    geneMeta: j.gene || {},
    klass: j.class,
    nCandidates: j.n_candidates,
    verdict: j.verdict,
    lines: (j.lines || []).map((l) => ({
      name: l.name || l.model_id,
      model_id: (l.model_id || "").toUpperCase(),
      lineage: cap(l.lineage) || "Unknown",
      subtype: l.subtype || l.disease || "—",
      score: l.core_score,          // real: within-gene percentile
      tier: l.confidence_tier,      // ordinal band from Stage 5
      reason: l.confidence_reason,
      nLayers: l.n_layers,
      nModalities: l.n_modalities,
      metadata: l.metadata || {},
      tracks: adaptTracks(l.tracks),
    })),
  };
}

const cap = (x) => (x ? String(x).replace(/\b\w/g, (c) => c.toUpperCase()) : x);

async function fetchIndex() {
  const r = await fetch(`${DATA_BASE}/index.json`);
  if (!r.ok) throw new Error(`index.json ${r.status}`);
  return (await r.json()).genes || [];
}

async function fetchGene(symbol) {
  const r = await fetch(`${DATA_BASE}/${symbol}.json`);
  if (!r.ok) return null;
  return adaptPayload(await r.json());
}

/* ------------------------------------------------------------------ */
/*  Helpers                                                            */
/* ------------------------------------------------------------------ */
/* The tier is an ordinal band produced by Stage 5 (high / moderate / low /
   unknown / prior). It is NOT derived from the score here — deriving a band
   from a number would re-introduce exactly the calibrated-confidence claim
   the pipeline refuses to make. */
const confTier = (tier) => ({
  high:     { label: "High confidence", cls: "high" },
  moderate: { label: "Moderate",        cls: "strong" },
  low:      { label: "Low",             cls: "mod" },
  prior:    { label: "Prior only",      cls: "mod" },
  unknown:  { label: "Unvalidated",     cls: "weak" },
}[tier] || { label: tier || "Unvalidated", cls: "weak" });
const pct = (c) => Math.round((c || 0) * 100);
const segAlpha = (s) => 0.28 + 0.72 * s;

/* ------------------------------------------------------------------ */
/*  Small building blocks                                              */
/* ------------------------------------------------------------------ */
function Mark() {
  return (
    <svg width="26" height="26" viewBox="0 0 26 26" fill="none" aria-hidden="true">
      <rect x="1.5" y="1.5" width="23" height="23" rx="7" stroke="var(--accent)" strokeWidth="1.6" />
      <path d="M6 17 L10 9 L13 14 L16 6 L20 12" stroke="var(--accent)" strokeWidth="1.7" strokeLinecap="round" strokeLinejoin="round" fill="none" />
      <circle cx="20" cy="12" r="1.7" fill="var(--accent)" />
    </svg>
  );
}

function EvidenceSpine({ tracks }) {
  return (
    <div className="spine" role="img" aria-label="Evidence signature">
      {TRACKS.map((tr) => {
        const s = tracks[tr.key]?.score || 0;
        const on = s > 0;
        return (
          <span key={tr.key} className={`seg ${on ? "on" : "off"}`}
            title={`${tr.label}: ${on ? pct(s) + "%" : "no signal"}`}
            style={on ? { background: tr.color, opacity: segAlpha(s) } : undefined} />
        );
      })}
    </div>
  );
}

/* ------------------------------------------------------------------ */
/*  View 1 — Ranked list row                                           */
/* ------------------------------------------------------------------ */
function CellLineRow({ line, rank, active, onClick }) {
  const tier = confTier(line.tier);
  return (
    <button className={`row ${active ? "active" : ""}`} onClick={onClick}>
      <span className="row-main">
        <span className="row-top">
          <span className="cl-name">{line.name}</span>
          <span className="cl-id mono">{line.model_id}</span>
        </span>
        <span className="row-sub">
          <span className="lineage-dot" style={{ background: LIN_COLOR(line.lineage) }} />
          {line.lineage} · {line.subtype}
        </span>
      </span>
      <EvidenceSpine tracks={line.tracks} />
      <span className="conf-cell">
        <span className={`conf-num mono ${tier.cls}`}>{pct(line.score)}</span>
        <span className="conf-track">
          <span className={`conf-fill ${tier.cls}`} style={{ width: pct(line.score) + "%" }} />
        </span>
      </span>
    </button>
  );
}

/* ------------------------------------------------------------------ */
/*  View 2 — Correlation network (gene at centre, distance = strength) */
/* ------------------------------------------------------------------ */
function useNetworkLayout(lines, W, H) {
  return useMemo(() => {
    const cx = W / 2, cy = H / 2;
    const rMin = 78, rMax = Math.min(W, H) / 2 - 52;
    const lineages = [...new Set(lines.map((l) => l.lineage))];
    const nodes = lines.map((l, i) => {
      const sec = lineages.indexOf(l.lineage);
      const a = ((sec + 0.5) / lineages.length) * Math.PI * 2 + ((i % 3) - 1) * 0.22;
      const target = rMin + (1 - l.score) * (rMax - rMin);
      return { ...l, x: cx + Math.cos(a) * target, y: cy + Math.sin(a) * target, target };
    });
    for (let it = 0; it < 220; it++) {
      for (let i = 0; i < nodes.length; i++) {
        for (let j = i + 1; j < nodes.length; j++) {
          let dx = nodes[j].x - nodes[i].x, dy = nodes[j].y - nodes[i].y;
          let d = Math.hypot(dx, dy) || 0.01;
          const min = 62;
          if (d < min) {
            const push = ((min - d) / d) * 0.5;
            dx *= push; dy *= push;
            nodes[i].x -= dx; nodes[i].y -= dy;
            nodes[j].x += dx; nodes[j].y += dy;
          }
        }
      }
      for (const n of nodes) {
        let dx = n.x - cx, dy = n.y - cy; const d = Math.hypot(dx, dy) || 0.01;
        const k = ((n.target - d) / d) * 0.10;
        n.x += dx * k; n.y += dy * k;
        n.x = Math.max(34, Math.min(W - 34, n.x));
        n.y = Math.max(30, Math.min(H - 30, n.y));
      }
    }
    return nodes;
  }, [lines, W, H]);
}

function NetworkGraph({ gene, lines, selected, onSelect }) {
  const W = 620, H = 470;
  const [hover, setHover] = useState(null);
  const nodes = useNetworkLayout(lines, W, H);
  const cx = W / 2, cy = H / 2;
  const lineages = [...new Set(lines.map((l) => l.lineage))];

  return (
    <div className="figure">
      <div className="fig-caption">
        Each cell line sits at a distance set by its <b>core score</b> for {gene.symbol} — a
        within-gene percentile, so the nearest lines are the highest-ranked of the{" "}
        {gene.nCandidates?.toLocaleString()} scored. Node size follows the same score; colour marks
        tissue lineage. Distance is <b>not</b> a correlation and not a probability.
      </div>
      <svg viewBox={`0 0 ${W} ${H}`} className="net-svg" role="img" aria-label={`Correlation network for ${gene.symbol}`}>
        {nodes.map((n) => {
          const dim = hover && hover !== n.model_id && selected?.model_id !== n.model_id;
          return (
            <line key={"e" + n.model_id} x1={cx} y1={cy} x2={n.x} y2={n.y}
              stroke="var(--accent)" strokeWidth={1 + n.score * 3.6}
              strokeLinecap="round"
              opacity={dim ? 0.06 : 0.14 + n.score * 0.5} />
          );
        })}
        <g>
          <circle cx={cx} cy={cy} r={34} fill="var(--accent)" />
          <circle cx={cx} cy={cy} r={34} fill="none" stroke="#fff" strokeOpacity="0.25" strokeWidth="1" />
          <text x={cx} y={cy - 4} textAnchor="middle" className="net-gene">{gene.symbol}</text>
          <text x={cx} y={cy + 11} textAnchor="middle" className="net-gene-sub">gene</text>
        </g>
        {nodes.map((n) => {
          const active = hover === n.model_id || selected?.model_id === n.model_id;
          const dim = hover && hover !== n.model_id && selected?.model_id !== n.model_id;
          const r = 9 + n.score * 13;
          return (
            <g key={n.model_id} className="net-node" opacity={dim ? 0.4 : 1}
              onMouseEnter={() => setHover(n.model_id)} onMouseLeave={() => setHover(null)}
              onClick={() => onSelect(n)}>
              <circle cx={n.x} cy={n.y} r={active ? r + 3 : r}
                fill={LIN_COLOR(n.lineage)} stroke={active ? "var(--ink)" : "#fff"} strokeWidth={active ? 2 : 1.5} />
              {active && <text x={n.x} y={n.y + 1} textAnchor="middle" dominantBaseline="central" className="net-conf">{pct(n.score)}</text>}
              <text x={n.x} y={n.y + r + 13} textAnchor="middle" className="net-label">{n.name}</text>
            </g>
          );
        })}
      </svg>
      <div className="fig-legend">
        {lineages.map((l) => (
          <span className="lg" key={l}><span className="lg-dot" style={{ background: LIN_COLOR(l) }} />{l}</span>
        ))}
      </div>
    </div>
  );
}

/* ------------------------------------------------------------------ */
/*  View 3 — Omics map (Gene → omics layers → cell lines)              */
/* ------------------------------------------------------------------ */
function OmicsMap({ gene, lines, selected, onSelect }) {
  const W = 620;
  const rows = Math.max(lines.length, TRACKS.length);
  const H = 60 + rows * 60;
  const [hoverT, setHoverT] = useState(null);
  const [hoverL, setHoverL] = useState(null);

  const geneX = 66, geneR = 108, geneY = H / 2;
  const trackCX = 310, trackW = 112;
  const lineX = 486;

  const tY = (i) => (TRACKS.length === 1 ? H / 2 : 46 + i * ((H - 92) / (TRACKS.length - 1)));
  const lY = (i) => (lines.length === 1 ? H / 2 : 46 + i * ((H - 92) / (lines.length - 1)));

  const trackTotals = TRACKS.map((tr) =>
    lines.reduce((a, l) => a + (l.tracks[tr.key]?.score || 0), 0));

  const path = (x1, y1, x2, y2) => {
    const mx = (x1 + x2) / 2;
    return `M${x1},${y1} C${mx},${y1} ${mx},${y2} ${x2},${y2}`;
  };

  const linkState = (trKey, modelId) => {
    if (hoverT) return hoverT === trKey ? "hot" : "cold";
    if (hoverL) return hoverL === modelId ? "hot" : "cold";
    if (selected) return selected.model_id === modelId ? "hot" : "cold";
    return "base";
  };

  return (
    <div className="figure">
      <div className="fig-caption">
        Which layers actually carry {gene.symbol}'s link to each model. Hover a <b>layer</b> to see
        every model it supports, or a <b>model</b> to see which layers back it. Thickness encodes the
        measured value where one exists — <b>expression</b> and <b>proteomics</b> as a within-gene
        percentile, <b>copy number</b> as deviation from diploid. <b>Mutation</b>, <b>fusion</b> and
        <b>signature</b> are categorical: a full-strength line means a driver call, a half-strength
        one means the layer measured this model without one.
      </div>
      <div className="track-legend">
        {TRACKS.map((tr) => (
          <span key={tr.key} className="tl-item">
            <span className="tl-dot" style={{ background: tr.color }} />
            {tr.label}
            <span className="tl-kind mono">{tr.kind === "call" ? "call" : "measured"}</span>
          </span>
        ))}
      </div>
      <svg viewBox={`0 0 ${W} ${H}`} className="omics-svg" role="img" aria-label={`Omics interaction map for ${gene.symbol}`}>
        {/* gene -> track links */}
        {TRACKS.map((tr, i) => {
          const w = 1 + (trackTotals[i] / Math.max(1, lines.length)) * 5;
          const on = !hoverT || hoverT === tr.key;
          return (
            <path key={"gt" + tr.key} d={path(geneX + geneR / 2, geneY, trackCX - trackW / 2, tY(i))}
              fill="none" stroke={tr.color} strokeWidth={w}
              opacity={on ? 0.32 : 0.06} />
          );
        })}
        {/* track -> line links */}
        {TRACKS.map((tr, i) =>
          lines.map((l) => {
            const s = l.tracks[tr.key]?.score || 0;
            if (s <= 0) return null;
            const st = linkState(tr.key, l.model_id);
            const op = st === "hot" ? 0.92 : st === "cold" ? 0.05 : 0.26;
            return (
              <path key={"tl" + tr.key + l.model_id}
                d={path(trackCX + trackW / 2, tY(i), lineX, lY(lines.indexOf(l)))}
                fill="none" stroke={tr.color} strokeWidth={1 + s * 5} strokeLinecap="round" opacity={op} />
            );
          })
        )}
        {/* gene node */}
        <g>
          <rect x={geneX - geneR / 2} y={geneY - 24} width={geneR} height={48} rx={11} fill="var(--accent)" />
          <text x={geneX} y={geneY - 3} textAnchor="middle" className="om-gene">{gene.symbol}</text>
          <text x={geneX} y={geneY + 13} textAnchor="middle" className="om-gene-sub">query gene</text>
        </g>
        {/* track nodes */}
        {TRACKS.map((tr, i) => {
          const on = !hoverT || hoverT === tr.key;
          return (
            <g key={tr.key} className="om-track" onMouseEnter={() => setHoverT(tr.key)} onMouseLeave={() => setHoverT(null)} opacity={on ? 1 : 0.45}>
              <rect x={trackCX - trackW / 2} y={tY(i) - 15} width={trackW} height={30} rx={8} fill={tr.color} />
              <text x={trackCX} y={tY(i) + 1} textAnchor="middle" dominantBaseline="central" className="om-track-label">{tr.label}</text>
            </g>
          );
        })}
        {/* cell line nodes */}
        {lines.map((l, j) => {
          const active = hoverL === l.model_id || selected?.model_id === l.model_id;
          return (
            <g key={l.model_id} className="om-line" onMouseEnter={() => setHoverL(l.model_id)} onMouseLeave={() => setHoverL(null)} onClick={() => onSelect(l)}>
              <circle cx={lineX} cy={lY(j)} r={active ? 8 : 6} fill={LIN_COLOR(l.lineage)} stroke={active ? "var(--ink)" : "#fff"} strokeWidth={active ? 2 : 1.5} />
              <text x={lineX + 14} y={lY(j) + 1} dominantBaseline="central" className={`om-line-label ${active ? "on" : ""}`}>{l.name}</text>
              <text x={lineX + 14} y={lY(j) + 14} dominantBaseline="central" className="om-line-sub mono">{pct(l.score)}</text>
            </g>
          );
        })}
      </svg>
    </div>
  );
}

/* ------------------------------------------------------------------ */
/*  Detail panel                                                       */
/* ------------------------------------------------------------------ */
function DetailPanel({ gene, line, onClose }) {
  if (!line) {
    return (
      <aside className="detail empty">
        <div className="detail-empty-inner">
          <div className="pulse" />
          <p>Select a cell line to read its full evidence.</p>
        </div>
      </aside>
    );
  }
  const tier = confTier(line.tier);
  const fired = TRACKS.filter((tr) => (line.tracks[tr.key]?.score || 0) > 0);
  const silent = TRACKS.filter((tr) => !(line.tracks[tr.key]?.score || 0));

  return (
    <aside className="detail">
      <div className="detail-head">
        <button className="close" onClick={onClose} aria-label="Close">×</button>
        <div className="d-name">{line.name}</div>
        <div className="d-id mono">{line.model_id}</div>
        <div className="chips">
          <span className="chip" style={{ background: LIN_COLOR(line.lineage) + "22", color: LIN_COLOR(line.lineage) }}>{line.lineage}</span>
          <span className="chip subtle">{line.subtype}</span>
        </div>
      </div>

      <div className="conf-block">
        <div className="conf-block-top">
          <span className="conf-label">Confidence for {gene.symbol}</span>
          <span className={`tier-pill ${tier.cls}`}>{tier.label}</span>
        </div>
        <div className="conf-big">
          <span className={`mono ${tier.cls}`}>{pct(line.score)}</span>
          <span className="conf-big-unit">/ 100</span>
        </div>
        <div className="conf-track big">
          <span className={`conf-fill ${tier.cls}`} style={{ width: pct(line.score) + "%" }} />
        </div>
      </div>

      <div className="section-label">Evidence readout</div>
      <div className="readout">
        {fired.map((tr) => {
          const d = line.tracks[tr.key];
          return (
            <div className="ev" key={tr.key}>
              <span className="ev-dot" style={{ background: tr.color }} />
              <div className="ev-body">
                <div className="ev-top">
                  <span className="ev-name">{tr.label}</span>
                  <span className="ev-score mono">measured</span>
                </div>
                <div className="ev-finding">{d.finding || "layer covers this cell line"}</div>
              </div>
            </div>
          );
        })}
        {silent.length > 0 && (
          <div className="silent">
            <span className="silent-label">Not measured</span>
            {silent.map((tr) => <span key={tr.key} className="silent-chip">{tr.label}</span>)}
          </div>
        )}
      </div>

      {line.reason && (
        <p className="method-note"><b>Why this tier:</b> {line.reason}</p>
      )}

      {line.metadata && Object.keys(line.metadata).length > 0 && (
        <div className="line-meta">
          <div className="ev-name" style={{ marginBottom: 6 }}>Cell line metadata</div>
          {Object.entries(line.metadata).map(([k, v]) => (
            <div key={k}><span className="lm-k">{k}</span>{String(v)}</div>
          ))}
        </div>
      )}

      <p className="method-note">
        <b>core score</b> is a within-gene percentile across the {gene.nCandidates?.toLocaleString()} cell
        lines scored for {gene.symbol}. It is comparable between cell lines for this gene, not between
        genes, and it is not a calibrated probability — the tier above carries the confidence judgement.
      </p>
    </aside>
  );
}


/* ------------------------------------------------------------------ */
/*  Verdict banner — the abstention states, carried in the payload      */
/*                                                                     */
/*  The pipeline exports a `verdict` block per gene (evidence_state.py).*/
/*  A gene that cannot be ranked says so HERE, above the results, and   */
/*  suppresses or greys the table -- rather than drawing a confident    */
/*  looking chart with a footnote. GAPDH is the worked example: it      */
/*  scores 0.9989 at only 2.8x top-vs-median.                           */
/* ------------------------------------------------------------------ */
function VerdictBanner({ verdict }) {
  if (!verdict || verdict.state === "RANKED") {
    if (verdict && verdict.detail && verdict.detail.length) {
      return (
        <div className="verdict verdict-ok">
          <div className="verdict-head">{verdict.headline}</div>
          <ul>{verdict.detail.map((d, i) => <li key={i}>{d}</li>)}</ul>
        </div>
      );
    }
    return null;
  }
  const cls = verdict.state === "NO_EVIDENCE" ? "verdict-none" : "verdict-warn";
  return (
    <div className={`verdict ${cls}`}>
      <div className="verdict-tag mono">{verdict.state.replace("_", " ")}</div>
      <div className="verdict-head">{verdict.headline}</div>
      <ul>{verdict.detail.map((d, i) => <li key={i}>{d}</li>)}</ul>
      {verdict.qualify_ranking && (
        <div className="verdict-foot">
          The ranking below is shown for completeness only. It is <b>not</b> a recommendation.
        </div>
      )}
    </div>
  );
}

/* ------------------------------------------------------------------ */
/*  Search box                                                         */
/* ------------------------------------------------------------------ */
function SearchBox({ value, onChange, onSubmit, big, index = [] }) {
  const [open, setOpen] = useState(false);
  const [hi, setHi] = useState(0);
  const boxRef = useRef(null);
  const matches = useMemo(() => {
    const q = value.trim().toUpperCase();
    const keys = index.map((g) => g.symbol);
    return q ? keys.filter((k) => k.includes(q)) : keys;
  }, [value, index]);
  useEffect(() => {
    const h = (e) => { if (boxRef.current && !boxRef.current.contains(e.target)) setOpen(false); };
    document.addEventListener("mousedown", h);
    return () => document.removeEventListener("mousedown", h);
  }, []);
  const pick = (k) => { onChange(k); onSubmit(k); setOpen(false); };
  return (
    <div className={`search ${big ? "big" : ""}`} ref={boxRef}>
      <svg className="search-icon" width="18" height="18" viewBox="0 0 18 18" fill="none">
        <circle cx="8" cy="8" r="5.5" stroke="currentColor" strokeWidth="1.6" />
        <path d="M12.2 12.2L16 16" stroke="currentColor" strokeWidth="1.6" strokeLinecap="round" />
      </svg>
      <input value={value} placeholder={big ? "Enter a gene — e.g. EGFR, BRAF, KRAS" : "Search a gene"}
        spellCheck={false} onFocus={() => setOpen(true)}
        onChange={(e) => { onChange(e.target.value); setOpen(true); setHi(0); }}
        onKeyDown={(e) => {
          if (e.key === "ArrowDown") { e.preventDefault(); setHi((h) => Math.min(h + 1, matches.length - 1)); }
          else if (e.key === "ArrowUp") { e.preventDefault(); setHi((h) => Math.max(h - 1, 0)); }
          else if (e.key === "Enter") { matches[hi] ? pick(matches[hi]) : onSubmit(value); setOpen(false); }
          else if (e.key === "Escape") setOpen(false);
        }} />
      {open && matches.length > 0 && (
        <div className="suggest">
          {matches.map((k, i) => (
            <button key={k} className={`sug ${i === hi ? "hi" : ""}`} onMouseEnter={() => setHi(i)} onClick={() => pick(k)}>
              <span className="sug-sym">{k}</span>
              <span className="sug-role">exported</span>
            </button>
          ))}
        </div>
      )}
    </div>
  );
}

/* ------------------------------------------------------------------ */
/*  App                                                                */
/* ------------------------------------------------------------------ */
const VIEWS = [
  { key: "list", label: "Ranked list" },
  { key: "network", label: "Network" },
  { key: "omics", label: "Omics map" },
];

export default function App() {
  const [query, setQuery] = useState("");
  const [gene, setGene] = useState(null);
  const [selected, setSelected] = useState(null);
  const [lineageFilter, setLineageFilter] = useState("All");
  const [minConf, setMinConf] = useState(0);
  const [view, setView] = useState("list");
  const [showScoring, setShowScoring] = useState(false);

  const [index, setIndex] = useState([]);
  const [geneData, setGeneData] = useState(null);
  const [loading, setLoading] = useState(false);
  const [loadError, setLoadError] = useState(null);

  // which genes have been exported. `python src/pipeline/export_web.py --genes ...`
  useEffect(() => {
    fetchIndex().then(setIndex).catch((e) => setLoadError(String(e)));
  }, []);

  useEffect(() => {
    if (!gene) { setGeneData(null); return; }
    let cancelled = false;
    setLoading(true); setLoadError(null);
    fetchGene(gene)
      .then((d) => { if (!cancelled) { setGeneData(d); setLoading(false); } })
      .catch((e) => { if (!cancelled) { setLoadError(String(e)); setLoading(false); } });
    return () => { cancelled = true; };
  }, [gene]);

  const submit = (term) => {
    const key = String(term).trim().toUpperCase();
    setGene(key); setSelected(null); setLineageFilter("All"); setMinConf(0); setView("list");
  };

  const lineages = useMemo(() => {
    if (!geneData) return [];
    return ["All", ...Array.from(new Set(geneData.lines.map((l) => l.lineage)))];
  }, [geneData]);

  const visible = useMemo(() => {
    if (!geneData) return [];
    return geneData.lines
      .filter((l) => lineageFilter === "All" || l.lineage === lineageFilter)
      .filter((l) => l.score * 100 >= minConf)
      .sort((a, b) => b.score - a.score);
  }, [geneData, lineageFilter, minConf]);

  return (
    <div className="app">
      <style>{CSS}</style>

      <header className="topbar">
        <button className="brand" onClick={() => { setGene(null); setQuery(""); setSelected(null); }}>
          <Mark />
          <span className="wordmark">GeneTrace<span className="wm-accent">AI</span></span>
        </button>
        {gene && <div className="topbar-search"><SearchBox value={query} onChange={setQuery} onSubmit={submit} index={index} /></div>}
        <div className="meta mono">20 sources · 2,145 models · 29.8M scored pairs</div>
      </header>

      {/* landing */}
      {!gene && (
        <main className="landing">
          <div className="landing-inner">
            <div className="eyebrow">Target-to-model discovery</div>
            <h1>Trace any human gene<br />to the cancer cell lines<br />that depend on it.</h1>
            <p className="lede">
              Enter a gene symbol. GeneTraceAI ranks cancer cell-line models by how strongly the
              evidence connects them to your gene — across mutations, expression, fusions, mutational
              signatures, and proteomics.
            </p>
            <div className="landing-search"><SearchBox value={query} onChange={setQuery} onSubmit={submit} big index={index} /></div>
            <div className="examples">
              <span className="ex-label">Try</span>
              {index.map((g) => <button key={g.symbol} className="ex-chip" onClick={() => { setQuery(g.symbol); submit(g.symbol); }}>{g.symbol}</button>)}
            </div>
            <div className="stat-strip">
              <div className="stat"><span className="stat-n mono">19,176</span><span className="stat-l">scored genes</span></div>
              <div className="stat"><span className="stat-n mono">2,145</span><span className="stat-l">cancer cell-line models</span></div>
              <div className="stat"><span className="stat-n mono">20</span><span className="stat-l">harmonised sources</span></div>
            </div>
          </div>
        </main>
      )}

      {gene && loading && (
        <main className="noresult"><div className="noresult-inner">
          <h2>Loading {gene}…</h2></div></main>
      )}

      {/* no result */}
      {gene && !loading && !geneData && (
        <main className="noresult">
          <div className="noresult-inner">
            <div className="nr-glyph">{gene.slice(0, 2)}</div>
            <h2>No models on record for “{gene}”.</h2>
            <p>We don’t yet have cell-line evidence for this gene. Check the symbol, or try one of the examples below.</p>
            <div className="examples center">
              {index.map((g) => <button key={g.symbol} className="ex-chip" onClick={() => { setQuery(g.symbol); submit(g.symbol); }}>{g.symbol}</button>)}
            </div>
          </div>
        </main>
      )}

      {/* results */}
      {geneData && (
        <main className="results">
          <div className="results-col">
            <div className="gene-head">
              <div className="gene-head-main">
                <h2 className="gene-sym">{geneData.symbol}</h2>
                <div className="gene-facts">
                  <span className="mono gene-ensg">{geneData.ensg}</span>
                  <span className={`role-pill ${geneData.role === "oncogene" ? "onco" : "tsg"}`}>{geneData.role}</span>
                  <span className="gene-chr mono">{geneData.klass || "unclassified"}</span>
                  {geneData.geneMeta?.["signal spread"] && (
                    <span className="gene-chr mono">spread: {geneData.geneMeta["signal spread"]}</span>
                  )}
                </div>
                <div className="gene-name">{geneData.name}</div>
              </div>
              <div className="gene-head-right">
                <div className="gene-count">
                  <span className="gc-n mono">{visible.length}</span>
                  <span className="gc-l">of {geneData.nCandidates?.toLocaleString()} scored</span>
                </div>
                <button className="scoring-btn" onClick={() => setShowScoring((s) => !s)}>
                  {showScoring ? "Hide" : "How is this ranked?"}
                </button>
              </div>
            </div>

            <VerdictBanner verdict={geneData.verdict} />

            {showScoring && (
              <div className="scoring-panel">
                <div className="sp-title">Ranking basis</div>
                <div className="sp-grid">
                  <div><b>Alteration</b><span>Activating hotspot, amplification or fusion for oncogenes; loss-of-function for tumour suppressors.</span></div>
                  <div><b>Lineage context</b><span>Whether the gene is the established driver in that tissue.</span></div>
                  <div><b>Expression</b><span>mRNA percentile across all cell-line models.</span></div>
                  <div><b>Proteomics</b><span>Protein-level corroboration where measured.</span></div>
                  <div><b>Convergence</b><span>How many independent tracks agree — more agreement lifts confidence.</span></div>
                </div>
                <p className="sp-note">Demo scores are illustrative, assembled from established cell-line biology (DepMap / CCLE / COSMIC-type evidence) — not live pipeline output.</p>
              </div>
            )}

            {/* view switcher */}
            <div className="viewbar">
              <div className="tabs">
                {VIEWS.map((v) => (
                  <button key={v.key} className={`tab ${view === v.key ? "on" : ""}`} onClick={() => setView(v.key)}>{v.label}</button>
                ))}
              </div>
              <div className="filters">
                <label className="fld">
                  <span>Lineage</span>
                  <select value={lineageFilter} onChange={(e) => setLineageFilter(e.target.value)}>
                    {lineages.map((l) => <option key={l}>{l}</option>)}
                  </select>
                </label>
                <label className="fld slider-fld">
                  <span>Min core score <b className="mono">{minConf}</b></span>
                  <input type="range" min="0" max="100" value={minConf} onChange={(e) => setMinConf(Number(e.target.value))} />
                </label>
              </div>
            </div>

            {/* active view */}
            {view === "list" && (
              <>
                <div className="list-head">
                  <span>Model</span>
                  <span className="lh-spine">Evidence</span>
                  <span className="lh-conf">Core score</span>
                </div>
                <div className="list">
                  {visible.length === 0 && <div className="list-empty">No models match these filters. Lower the confidence threshold.</div>}
                  {visible.map((line, i) => (
                    <CellLineRow key={line.model_id} line={line} rank={i + 1}
                      active={selected?.model_id === line.model_id} onClick={() => setSelected(line)} />
                  ))}
                </div>
                <div className="track-legend">
                  {TRACKS.map((tr) => <span className="lg" key={tr.key}><span className="lg-dot" style={{ background: tr.color }} />{tr.label}</span>)}
                </div>
              </>
            )}
            {view === "network" && visible.length > 0 && (
              <NetworkGraph gene={geneData} lines={visible} selected={selected} onSelect={setSelected} />
            )}
            {view === "omics" && visible.length > 0 && (
              <OmicsMap gene={geneData} lines={visible} selected={selected} onSelect={setSelected} />
            )}
            {view !== "list" && visible.length === 0 && (
              <div className="list-empty">No models match these filters. Lower the confidence threshold.</div>
            )}
          </div>

          <DetailPanel gene={geneData} line={selected} onClose={() => setSelected(null)} />
        </main>
      )}
    </div>
  );
}

/* ------------------------------------------------------------------ */
/*  Styles                                                             */
/* ------------------------------------------------------------------ */
const CSS = `
.track-legend{display:flex;flex-wrap:wrap;gap:12px;margin:2px 0 10px;font-size:11.5px;color:var(--muted);}
.tl-item{display:inline-flex;align-items:center;gap:5px;}
.tl-dot{width:9px;height:9px;border-radius:3px;display:inline-block;}
.tl-kind{opacity:.55;font-size:10px;padding:1px 5px;border:1px solid currentColor;border-radius:3px;}

.verdict { border-radius: 10px; padding: 14px 16px; margin: 14px 0 18px; font-size: 13.5px; line-height: 1.5; }
.verdict-tag { font-size: 10.5px; letter-spacing: .09em; text-transform: uppercase; opacity: .85; margin-bottom: 6px; }
.verdict-head { font-weight: 650; font-size: 14.5px; margin-bottom: 6px; }
.verdict ul { margin: 0; padding-left: 18px; }
.verdict li { margin: 3px 0; opacity: .92; }
.verdict-foot { margin-top: 9px; padding-top: 8px; border-top: 1px solid currentColor; opacity: .8; }
.verdict-warn { background: rgba(245,166,35,.10); border: 1px solid rgba(245,166,35,.45); color: #8a5a00; }
.verdict-none { background: rgba(93,113,131,.10); border: 1px solid rgba(93,113,131,.40); color: #44515e; }
.verdict-ok   { background: rgba(48,163,108,.08); border: 1px solid rgba(48,163,108,.35); color: #1d6b48; }
.line-meta { margin-top: 10px; font-size: 12px; line-height: 1.7; }
.line-meta .lm-k { opacity: .6; display: inline-block; min-width: 132px; }

@import url('https://fonts.googleapis.com/css2?family=Space+Grotesk:wght@500;600;700&family=IBM+Plex+Sans:wght@400;500;600&family=IBM+Plex+Mono:wght@400;500&display=swap');

.app{
  --paper:#F4F6F8; --surface:#FFFFFF; --ink:#0F2130; --ink-2:#213445;
  --muted:#5D7183; --faint:#93A3B1; --line:#E4E9ED; --line-2:#EEF2F5;
  --accent:#0B7C86; --accent-d:#075A62; --accent-t:#E6F1F1;
  --high:#0B7C86; --strong:#2E8A8A; --mod:#C79A2E; --weak:#9AA9B4;
  --onco:#C0334D; --tsg:#3A5BB8;
  --disp:'Space Grotesk',sans-serif; --body:'IBM Plex Sans',sans-serif; --mono:'IBM Plex Mono',monospace;
  font-family:var(--body); color:var(--ink); background:var(--paper);
  min-height:100vh; -webkit-font-smoothing:antialiased;
}
.app *{box-sizing:border-box;}
.app .mono{font-family:var(--mono);}
.app button{font-family:inherit; cursor:pointer;}

/* top bar */
.topbar{display:flex; align-items:center; gap:20px; padding:14px 26px; background:var(--surface); border-bottom:1px solid var(--line); position:sticky; top:0; z-index:20;}
.brand{display:flex; align-items:center; gap:10px; background:none; border:none; padding:0;}
.wordmark{font-family:var(--disp); font-weight:700; font-size:19px; letter-spacing:-.02em; color:var(--ink);}
.wm-accent{color:var(--accent);}
.topbar-search{flex:0 1 340px; margin-left:8px;}
.meta{margin-left:auto; font-size:12px; color:var(--faint); letter-spacing:-.01em; display:flex; align-items:center; gap:8px;}
.demo-pill{font-family:var(--body); font-size:10px; letter-spacing:.04em; text-transform:uppercase; color:var(--mod); background:#FAF1DC; padding:2px 7px; border-radius:10px;}

/* search */
.search{position:relative; display:flex; align-items:center;}
.search-icon{position:absolute; left:12px; color:var(--faint); pointer-events:none;}
.search input{width:100%; font-family:var(--body); font-size:14px; color:var(--ink); padding:9px 12px 9px 36px; border:1px solid var(--line); border-radius:9px; background:var(--surface); outline:none; transition:border-color .15s, box-shadow .15s;}
.search input:focus{border-color:var(--accent); box-shadow:0 0 0 3px var(--accent-t);}
.search.big input{font-size:17px; padding:16px 16px 16px 46px; border-radius:13px;}
.search.big .search-icon{left:17px; width:20px; height:20px;}
.suggest{position:absolute; top:calc(100% + 6px); left:0; right:0; z-index:30; background:var(--surface); border:1px solid var(--line); border-radius:11px; box-shadow:0 12px 34px rgba(15,33,48,.12); overflow:hidden; padding:5px;}
.sug{display:flex; align-items:center; gap:10px; width:100%; text-align:left; padding:9px 11px; border:none; background:none; border-radius:8px; color:var(--ink);}
.sug.hi{background:var(--accent-t);}
.sug-sym{font-family:var(--disp); font-weight:600; font-size:15px;}
.sug-role{font-size:12px; color:var(--muted);}
.sug-count{margin-left:auto; font-size:12px; color:var(--faint);}

/* landing */
.landing{display:flex; justify-content:center; padding:8vh 26px 60px;}
.landing-inner{width:100%; max-width:780px;}
.eyebrow{font-size:12px; letter-spacing:.14em; text-transform:uppercase; color:var(--accent); font-weight:600; margin-bottom:18px;}
.landing h1{font-family:var(--disp); font-weight:700; font-size:clamp(32px,5vw,50px); line-height:1.05; letter-spacing:-.03em; margin:0 0 20px; color:var(--ink);}
.lede{font-size:17px; line-height:1.55; color:var(--muted); max-width:600px; margin:0 0 30px;}
.landing-search{margin-bottom:18px;}
.examples{display:flex; align-items:center; gap:9px; flex-wrap:wrap;}
.examples.center{justify-content:center; margin-top:8px;}
.ex-label{font-size:13px; color:var(--faint); margin-right:2px;}
.ex-chip{font-family:var(--mono); font-size:13px; color:var(--ink-2); padding:6px 13px; border:1px solid var(--line); border-radius:20px; background:var(--surface); transition:all .14s;}
.ex-chip:hover{border-color:var(--accent); color:var(--accent); background:var(--accent-t);}
.stat-strip{display:flex; gap:40px; margin-top:52px; padding-top:26px; border-top:1px solid var(--line);}
.stat{display:flex; flex-direction:column; gap:3px;}
.stat-n{font-size:26px; font-weight:500; color:var(--ink); letter-spacing:-.02em;}
.stat-l{font-size:13px; color:var(--muted);}

/* results layout */
.results{display:grid; grid-template-columns:1fr 380px; min-height:calc(100vh - 61px);}
.results-col{padding:26px 26px 40px; min-width:0; border-right:1px solid var(--line);}

/* gene header */
.gene-head{display:flex; align-items:flex-start; justify-content:space-between; gap:20px; margin-bottom:18px;}
.gene-sym{font-family:var(--disp); font-weight:700; font-size:40px; letter-spacing:-.03em; margin:0 0 8px; line-height:1;}
.gene-facts{display:flex; align-items:center; gap:10px; flex-wrap:wrap; margin-bottom:8px;}
.gene-ensg{font-size:13px; color:var(--muted);}
.role-pill{font-size:11px; font-weight:600; letter-spacing:.03em; text-transform:uppercase; padding:3px 9px; border-radius:5px;}
.role-pill.onco{color:var(--onco); background:#FBE9EC;}
.role-pill.tsg{color:var(--tsg); background:#E9EEFA;}
.gene-chr{font-size:12px; color:var(--faint);}
.gene-name{font-size:14px; color:var(--muted);}
.gene-head-right{display:flex; flex-direction:column; align-items:flex-end; gap:8px; flex-shrink:0;}
.gene-count{text-align:right;}
.gc-n{display:block; font-size:30px; font-weight:500; color:var(--accent); line-height:1;}
.gc-l{font-size:12px; color:var(--muted);}
.scoring-btn{font-size:12px; color:var(--muted); background:none; border:1px solid var(--line); padding:5px 11px; border-radius:16px; transition:all .14s;}
.scoring-btn:hover{border-color:var(--accent); color:var(--accent);}

/* scoring panel */
.scoring-panel{background:var(--surface); border:1px solid var(--line); border-radius:12px; padding:18px 20px; margin-bottom:20px;}
.sp-title{font-size:11px; letter-spacing:.08em; text-transform:uppercase; color:var(--faint); margin-bottom:12px;}
.sp-grid{display:grid; grid-template-columns:repeat(auto-fit,minmax(180px,1fr)); gap:14px 22px;}
.sp-grid div{display:flex; flex-direction:column; gap:3px;}
.sp-grid b{font-size:13px; color:var(--ink); font-weight:600;}
.sp-grid span{font-size:12.5px; color:var(--muted); line-height:1.45;}
.sp-note{font-size:12px; color:var(--faint); line-height:1.5; margin:16px 0 0; padding-top:12px; border-top:1px solid var(--line-2);}

/* view bar */
.viewbar{display:flex; align-items:flex-end; justify-content:space-between; gap:20px; flex-wrap:wrap; margin-bottom:18px;}
.tabs{display:inline-flex; background:var(--line-2); border-radius:10px; padding:3px;}
.tab{font-size:13px; font-weight:500; color:var(--muted); background:none; border:none; padding:7px 15px; border-radius:8px; transition:all .14s;}
.tab.on{background:var(--surface); color:var(--ink); box-shadow:0 1px 3px rgba(15,33,48,.10);}
.filters{display:flex; gap:22px; align-items:flex-end;}
.fld{display:flex; flex-direction:column; gap:6px; font-size:12px; color:var(--muted);}
.fld select{font-family:var(--body); font-size:13px; color:var(--ink); padding:7px 10px; border:1px solid var(--line); border-radius:8px; background:var(--surface); outline:none;}
.fld select:focus{border-color:var(--accent);}
.slider-fld{min-width:150px;}
.slider-fld b{color:var(--accent);}
.slider-fld input[type=range]{width:100%; accent-color:var(--accent);}

/* list */
.list-head{display:grid; grid-template-columns:1fr 110px 128px; gap:16px; padding:0 14px 9px; font-size:11px; letter-spacing:.06em; text-transform:uppercase; color:var(--faint); border-bottom:1px solid var(--line);}
.lh-spine{text-align:center;} .lh-conf{text-align:right;}
.list{display:flex; flex-direction:column;}
.list-empty{padding:32px 14px; color:var(--muted); font-size:14px;}
.row{display:grid; grid-template-columns:1fr 110px 128px; gap:16px; align-items:center; width:100%; text-align:left; background:none; border:none; padding:14px; border-bottom:1px solid var(--line-2); transition:background .12s;}
.row:hover{background:var(--surface);}
.row.active{background:var(--accent-t); box-shadow:inset 3px 0 0 var(--accent);}
.row-main{display:flex; flex-direction:column; gap:4px; min-width:0;}
.row-top{display:flex; align-items:baseline; gap:10px;}
.cl-name{font-family:var(--disp); font-weight:600; font-size:16px; letter-spacing:-.01em; color:var(--ink);}
.cl-id{font-size:11px; color:var(--faint); flex-shrink:0;}
.row-sub{display:flex; align-items:center; gap:7px; font-size:12.5px; color:var(--muted);}
.lineage-dot{width:7px; height:7px; border-radius:50%;}
.spine{display:flex; gap:3px; justify-content:center;}
.spine .seg{width:16px; height:9px;}
.seg{border-radius:2px; display:block; background:var(--line);}
.seg.off{background:transparent; box-shadow:inset 0 0 0 1px var(--line);}
.conf-cell{display:flex; flex-direction:column; align-items:flex-end; gap:5px;}
.conf-num{font-size:19px; font-weight:500; line-height:1;}
.conf-num.high{color:var(--high);} .conf-num.strong{color:var(--strong);} .conf-num.mod{color:var(--mod);} .conf-num.weak{color:var(--weak);}
.conf-track{width:100%; height:5px; border-radius:3px; background:var(--line-2); overflow:hidden;}
.conf-track.big{height:8px;}
.conf-fill{display:block; height:100%; border-radius:3px;}
.conf-fill.high{background:var(--high);} .conf-fill.strong{background:var(--strong);} .conf-fill.mod{background:var(--mod);} .conf-fill.weak{background:var(--weak);}
.track-legend{display:flex; gap:16px; flex-wrap:wrap; margin-top:16px; padding-top:14px; border-top:1px solid var(--line-2);}
.lg{display:flex; align-items:center; gap:6px; font-size:12px; color:var(--muted);}
.lg-dot{width:9px; height:9px; border-radius:3px;}

/* figures */
.figure{padding-top:4px;}
.fig-caption{font-size:13px; color:var(--muted); line-height:1.5; max-width:560px; margin-bottom:6px;}
.fig-caption b{color:var(--ink);}
.net-svg,.omics-svg{width:100%; height:auto; display:block; user-select:none;}
.fig-legend{display:flex; gap:14px; flex-wrap:wrap; margin-top:8px;}
.net-gene{font-family:var(--disp); font-weight:700; font-size:16px; fill:#fff;}
.net-gene-sub{font-family:var(--body); font-size:9px; fill:#fff; opacity:.7; letter-spacing:.05em; text-transform:uppercase;}
.net-node{cursor:pointer;}
.net-label{font-family:var(--body); font-size:11px; fill:var(--ink-2);}
.net-conf{font-family:var(--mono); font-size:10px; fill:#fff; font-weight:500;}

.om-gene{font-family:var(--disp); font-weight:700; font-size:16px; fill:#fff;}
.om-gene-sub{font-family:var(--body); font-size:9px; fill:#fff; opacity:.75; letter-spacing:.04em; text-transform:uppercase;}
.om-track{cursor:pointer;}
.om-track-label{font-family:var(--body); font-weight:600; font-size:12px; fill:#fff;}
.om-line{cursor:pointer;}
.om-line-label{font-family:var(--body); font-size:12px; fill:var(--ink-2);}
.om-line-label.on{fill:var(--ink); font-weight:600;}
.om-line-sub{font-size:9px; fill:var(--faint);}

/* detail panel */
.detail{background:var(--surface); padding:26px 24px 30px; position:sticky; top:61px; align-self:start; max-height:calc(100vh - 61px); overflow-y:auto;}
.detail.empty{display:flex; align-items:center; justify-content:center; color:var(--faint);}
.detail-empty-inner{text-align:center; font-size:14px; max-width:200px;}
.pulse{width:34px; height:34px; margin:0 auto 14px; border-radius:10px; border:2px solid var(--line); position:relative;}
.pulse::after{content:""; position:absolute; inset:9px; border-radius:50%; border:2px solid var(--accent); opacity:.5;}
.detail-head{position:relative; padding-bottom:18px; border-bottom:1px solid var(--line); margin-bottom:20px;}
.close{position:absolute; top:-6px; right:-4px; font-size:22px; line-height:1; color:var(--faint); background:none; border:none; padding:4px 8px;}
.close:hover{color:var(--ink);}
.d-name{font-family:var(--disp); font-weight:700; font-size:26px; letter-spacing:-.02em; line-height:1.05;}
.d-id{font-size:12px; color:var(--faint); margin-top:4px;}
.chips{display:flex; gap:7px; margin-top:12px; flex-wrap:wrap;}
.chip{font-size:12px; font-weight:500; padding:4px 10px; border-radius:6px;}
.chip.subtle{background:var(--line-2); color:var(--muted); font-weight:400;}
.conf-block{margin-bottom:24px;}
.conf-block-top{display:flex; align-items:center; justify-content:space-between; margin-bottom:8px;}
.conf-label{font-size:13px; color:var(--muted);}
.tier-pill{font-size:11px; font-weight:600; padding:3px 9px; border-radius:20px;}
.tier-pill.high{color:var(--high); background:var(--accent-t);}
.tier-pill.strong{color:var(--strong); background:#E4F0F0;}
.tier-pill.mod{color:var(--mod); background:#FAF1DC;}
.tier-pill.weak{color:var(--weak); background:var(--line-2);}
.conf-big{display:flex; align-items:baseline; gap:7px; margin-bottom:9px;}
.conf-big .mono{font-size:44px; font-weight:500; line-height:1; letter-spacing:-.02em;}
.conf-big .mono.high{color:var(--high);} .conf-big .mono.strong{color:var(--strong);} .conf-big .mono.mod{color:var(--mod);} .conf-big .mono.weak{color:var(--weak);}
.conf-big-unit{font-size:15px; color:var(--faint);}
.section-label{font-size:11px; letter-spacing:.08em; text-transform:uppercase; color:var(--faint); margin-bottom:12px;}
.readout{display:flex; flex-direction:column; gap:16px;}
.ev{display:flex; gap:11px;}
.ev-dot{width:10px; height:10px; border-radius:3px; margin-top:4px; flex-shrink:0;}
.ev-body{flex:1; min-width:0;}
.ev-top{display:flex; align-items:baseline; justify-content:space-between; margin-bottom:3px;}
.ev-name{font-size:13.5px; font-weight:600; color:var(--ink);}
.ev-score{font-size:13px; color:var(--muted);}
.ev-finding{font-size:13px; color:var(--muted); line-height:1.45; margin-bottom:7px;}
.ev-bar{height:4px; border-radius:3px; background:var(--line-2); overflow:hidden;}
.ev-bar span{display:block; height:100%; border-radius:3px;}
.silent{display:flex; align-items:center; gap:7px; flex-wrap:wrap; padding-top:6px; margin-top:2px; border-top:1px dashed var(--line);}
.silent-label{font-size:11px; text-transform:uppercase; letter-spacing:.05em; color:var(--faint);}
.silent-chip{font-size:11.5px; color:var(--faint); padding:2px 8px; border:1px solid var(--line); border-radius:12px;}
.detail-actions{display:flex; gap:9px; margin:24px 0 14px;}
.btn{font-size:13px; font-weight:500; padding:10px 15px; border-radius:9px; border:1px solid transparent; transition:all .14s;}
.btn.primary{background:var(--accent); color:#fff; flex:1;}
.btn.primary:hover{background:var(--accent-d);}
.btn.ghost{background:var(--surface); border-color:var(--line); color:var(--ink-2);}
.btn.ghost:hover{border-color:var(--accent); color:var(--accent);}
.method-note{font-size:12px; line-height:1.5; color:var(--faint); margin:0;}

/* no result */
.noresult{display:flex; justify-content:center; padding:12vh 26px;}
.noresult-inner{max-width:440px; text-align:center;}
.nr-glyph{width:64px; height:64px; margin:0 auto 20px; border-radius:16px; display:flex; align-items:center; justify-content:center; font-family:var(--disp); font-weight:700; font-size:24px; color:var(--faint); background:var(--surface); border:1px solid var(--line);}
.noresult h2{font-family:var(--disp); font-weight:600; font-size:24px; letter-spacing:-.02em; margin:0 0 10px;}
.noresult p{font-size:15px; color:var(--muted); line-height:1.55; margin:0 0 22px;}

/* responsive */
@media (max-width:940px){
  .results{grid-template-columns:1fr;}
  .results-col{border-right:none;}
  .detail{position:fixed; inset:auto 0 0 0; top:auto; max-height:78vh; border-top:1px solid var(--line); box-shadow:0 -12px 40px rgba(15,33,48,.14); z-index:40; border-radius:18px 18px 0 0;}
  .detail.empty{display:none;}
  .topbar-search{display:none;}
}
@media (max-width:560px){
  .topbar{padding:12px 16px;} .meta{display:none;}
  .results-col{padding:20px 16px 40px;}
  .list-head{grid-template-columns:1fr 84px;} .lh-conf{display:none;}
  .row{grid-template-columns:1fr 84px;}
  .conf-cell{display:none;}
  .gene-sym{font-size:32px;} .landing{padding:6vh 18px 40px;}
  .stat-strip{gap:22px; flex-wrap:wrap;}
  .viewbar{gap:14px;}
}
@media (prefers-reduced-motion:reduce){ .app *{transition:none !important;} }
`;
