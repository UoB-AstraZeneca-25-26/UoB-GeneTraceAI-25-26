/* ------------------------------------------------------------------ */
/*  MethodologyChart — the 7-step scoring pipeline as a vertical       */
/*  pipeline of clickable nodes.                                       */
/*                                                                     */
/*  Fed by QueryResponse.methodology from POST /v1/agent/query, which  */
/*  the agent returns as {steps:[{key,value,formula,status}, x7]}.     */
/*  Each node shows its label and formula; the 20-25 word explanation  */
/*  is revealed on click. A step that did not apply to this gene       */
/*  (TSG inversion on an oncogene) is drawn dashed and muted rather    */
/*  than dropped -- the pipeline stays auditable either way.           */
/* ------------------------------------------------------------------ */
import React, { useState } from "react";
import { AnimatePresence, motion, useReducedMotion } from "framer-motion";
import "./MethodologyChart.css";

/* ---- animation variants ---- */
/* Nodes and the connectors between them are siblings in one stagger, so a
   line always draws in the gap between the two nodes it joins. */
const container = {
  hidden: {},
  show: { transition: { staggerChildren: 0.06 } },
};

const item = {
  hidden: { opacity: 0, x: -30 },
  show: {
    opacity: 1,
    x: 0,
    transition: { type: "spring", stiffness: 300, damping: 24 },
  },
};

const linkWrap = { hidden: { opacity: 1 }, show: { opacity: 1 } };

const linkLine = {
  hidden: { pathLength: 0, opacity: 0 },
  show: {
    pathLength: 1,
    opacity: 1,
    transition: { duration: 0.24, ease: "easeOut" },
  },
};

/* ---- badge glyphs ---- */
function Glyph({ status, index, open }) {
  if (status === "skipped") {
    return (
      <svg viewBox="0 0 16 16" width="13" height="13" aria-hidden="true">
        <circle cx="8" cy="8" r="6" fill="none" stroke="currentColor" strokeWidth="1.6" />
        <line x1="4.2" y1="11.8" x2="11.8" y2="4.2" stroke="currentColor" strokeWidth="1.6" strokeLinecap="round" />
      </svg>
    );
  }
  if (status === "inverted") {
    return (
      <svg viewBox="0 0 16 16" width="13" height="13" aria-hidden="true">
        <path d="M3 5.5h8M8.5 3l2.5 2.5L8.5 8" fill="none" stroke="currentColor" strokeWidth="1.6" strokeLinecap="round" strokeLinejoin="round" />
        <path d="M13 10.5H5M7.5 8L5 10.5 7.5 13" fill="none" stroke="currentColor" strokeWidth="1.6" strokeLinecap="round" strokeLinejoin="round" />
      </svg>
    );
  }
  if (open) {
    return (
      <svg viewBox="0 0 16 16" width="13" height="13" aria-hidden="true">
        <path d="M3.5 8.5l3 3 6-7" fill="none" stroke="currentColor" strokeWidth="2" strokeLinecap="round" strokeLinejoin="round" />
      </svg>
    );
  }
  return <span className="mc-num">{index}</span>;
}

/* Why a step is greyed out. Only the TSG inversion can legitimately not
   apply, so name that reason when we can and stay generic otherwise. */
function skipReason(step, geneName) {
  const subject = geneName || "This gene";
  if (/inver|tsg|suppress/i.test(step.key)) {
    return `${subject} is not a tumour suppressor — inversion applies only to TSGs.`;
  }
  return `This step did not apply to ${geneName || "this gene"}.`;
}

const STATUS_LABEL = { active: "active", skipped: "skipped", inverted: "inverted" };

function Step({ step, index, open, onToggle, geneName, reduced }) {
  const status = STATUS_LABEL[step.status] ? step.status : "active";
  const hover = reduced
    ? {}
    : { scale: 1.02, boxShadow: "0 4px 12px rgba(13,148,136,0.15)" };

  return (
    <motion.div className="mc-step" variants={item}>
      <motion.button
        type="button"
        className={`mc-node mc-${status}${open ? " is-open" : ""}`}
        onClick={onToggle}
        whileHover={hover}
        whileTap={reduced ? {} : { scale: 0.995 }}
        aria-expanded={open}
      >
        <span className="mc-badge">
          <Glyph status={status} index={index} open={open} />
        </span>
        <span className="mc-body">
          <span className="mc-key">{step.key}</span>
          {step.formula && <span className="mc-formula">{step.formula}</span>}
        </span>
        <span
          className="mc-tag"
          title={status === "skipped" ? skipReason(step, geneName) : undefined}
        >
          {status}
        </span>
      </motion.button>

      <AnimatePresence initial={false}>
        {open && (
          <motion.div
            className="mc-detail"
            initial={{ height: 0, opacity: 0 }}
            animate={{ height: "auto", opacity: 1 }}
            exit={{ height: 0, opacity: 0 }}
            transition={{ duration: reduced ? 0 : 0.25 }}
          >
            <p className="mc-detail-text">{step.value}</p>
            {status === "skipped" && (
              <p className="mc-detail-note">{skipReason(step, geneName)}</p>
            )}
          </motion.div>
        )}
      </AnimatePresence>
    </motion.div>
  );
}

function Connector() {
  return (
    <motion.svg
      className="mc-link"
      variants={linkWrap}
      width="2"
      height="22"
      viewBox="0 0 2 22"
      preserveAspectRatio="none"
      aria-hidden="true"
    >
      <motion.line
        x1="1"
        y1="0"
        x2="1"
        y2="22"
        variants={linkLine}
        stroke="#cbd5e1"
        strokeWidth="2"
        strokeLinecap="round"
      />
    </motion.svg>
  );
}

export default function MethodologyChart({ steps, geneName }) {
  const [openKey, setOpenKey] = useState(null);
  const reduced = useReducedMotion();

  if (!steps || steps.length === 0) return null;

  return (
    <div className="mc">
      <div className="mc-head">
        Scoring methodology{geneName ? ` — ${geneName}` : ""}
      </div>
      <motion.div
        className="mc-pipe"
        variants={container}
        initial={reduced ? "show" : "hidden"}
        animate="show"
      >
        {steps.map((step, i) => (
          <React.Fragment key={`${i}-${step.key}`}>
            <Step
              step={step}
              index={i + 1}
              geneName={geneName}
              reduced={reduced}
              open={openKey === i}
              onToggle={() => setOpenKey(openKey === i ? null : i)}
            />
            {i < steps.length - 1 && <Connector />}
          </React.Fragment>
        ))}
      </motion.div>
      <p className="mc-foot">Click a step to read what it does.</p>
    </div>
  );
}

/* Shown while the agent is thinking (~8s). Same shape as the real chart, so
   the panel does not jump when the steps land. */
export function MethodologyChartSkeleton({ count = 7 }) {
  const reduced = useReducedMotion();
  return (
    <div className="mc">
      <div className="mc-head">Scoring methodology</div>
      <div className="mc-pipe">
        {Array.from({ length: count }, (_, i) => (
          <React.Fragment key={i}>
            <motion.div
              className="mc-skeleton"
              animate={reduced ? undefined : { backgroundPosition: ["200% 0", "-200% 0"] }}
              transition={{ duration: 1.5, repeat: Infinity, ease: "linear" }}
            />
            {i < count - 1 && <div className="mc-link mc-link-idle" />}
          </React.Fragment>
        ))}
      </div>
      <p className="mc-foot">Reading the pipeline…</p>
    </div>
  );
}
