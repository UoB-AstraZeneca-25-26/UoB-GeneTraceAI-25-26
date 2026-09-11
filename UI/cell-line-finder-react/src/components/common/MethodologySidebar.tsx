import React, { useEffect, useState } from 'react';
import { X } from 'lucide-react';
import { MethodologyStep } from '../../types';

interface MethodologySidebarProps {
  isOpen: boolean;
  onClose: () => void;
  geneName: string;
  ensgId: string;
  modelId: string;
  cellLineName: string;
}

// Static walkthrough of the general core_score method (architecture/Scoring/
// core_score.py, confidence_tiers.py, driver_routing.py) -- NOT a per-(gene,
// line) live computation. This used to call the agent for a bespoke
// explanation of this exact pair's score, but that endpoint wasn't
// returning real values, so this is the actual scoring pipeline's logic in
// plain terms instead: the same seven steps for every pair, not this pair's
// specific numbers.
const METHODOLOGY_STEPS: MethodologyStep[] = [
  {
    key: 'RNA standardization',
    formula: 'rna_std_z = (rna_z − median) / SD, per (gene, n_sources)',
    value:
      'RNA evidence is combined across up to 3 sources (DepMap, HPA, GEO), then standardized within each ' +
      '(gene, number-of-sources) group rather than globally per gene — the spread of RNA measurements genuinely ' +
      'differs depending on how many sources agree on a line, and treating that as one pool would over- or ' +
      'under-state confidence. Genes with too few lines in a group borrow strength from that gene\'s overall ' +
      'spread instead of using a noisy small-sample estimate. Capped at ±5 standard deviations either way.',
  },
  {
    key: 'Protein residualization',
    formula: 'prot_resid = prot_z − β·rna_z,  β = ρ × (SD_protein / SD_rna)',
    value:
      'Protein and RNA levels for the same gene are correlated, so protein evidence is first adjusted to remove ' +
      'the part that RNA already explains — what\'s left (the "residual") is the protein signal that\'s genuinely ' +
      'independent, not double-counting the same underlying biology. The correlation used per gene is blended ' +
      'toward a panel-wide typical value when a gene has few shared lines to estimate it from, so one noisy gene ' +
      'doesn\'t get an extreme, unreliable adjustment.',
  },
  {
    key: 'Weighted combination',
    formula: 'core_z = (W_RNA·rna_std_z + W_PROT·prot_resid_z) / √(W_RNA² + W_PROT²)',
    value:
      'The two standardized, independent signals are combined into one z-score. Each arm\'s weight reflects how ' +
      'much independent evidence it actually carries — RNA draws on up to 3 sources, protein on up to 2 — so an ' +
      'arm backed by more corroborating measurements counts for more, not an arbitrary 50/50 split. A line with ' +
      'only one arm available (no protein data, or a gene with no protein platform coverage) still gets scored, ' +
      'just from RNA alone, scaled consistently with how the two-arm case weights RNA.',
  },
  {
    key: 'Score mapping',
    formula: 'core_score = Φ(core_z)',
    value:
      'The combined z-score is mapped through the normal cumulative distribution function onto a 0-1 scale — ' +
      'this is what makes scores comparable across genes with very different raw z-score ranges. It\'s a ' +
      'relabelling, not a re-weighting: it doesn\'t change which lines rank higher than others, only the scale ' +
      'the ranking is displayed on.',
  },
  {
    key: 'Stratum ranking',
    formula: 'stratum_rank = rank(core_score) within (gene, lineage)',
    value:
      'Lines are ranked against other lines of the same tissue/lineage for this gene, not against the whole ' +
      'panel — a line that looks unremarkable pooled across all tissues can still be the clear best choice ' +
      'within its own lineage, which is usually the more useful comparison for picking a model.',
  },
  {
    key: 'Driver-alteration evidence',
    formula: 'has_driver_alteration = mutation_driver OR fusion_driver OR CNA_alteration',
    value:
      'Independently of the expression score above, each (gene, line) pair is checked for a known driver-class ' +
      'mutation, gene fusion, or copy-number alteration. This is genomic evidence, not expression evidence — a ' +
      'line can have a driver alteration with a modest expression score, or a high expression score with no ' +
      'detected alteration; the confidence tier below is what reconciles the two.',
  },
  {
    key: 'Confidence tier',
    formula: 'HIGH / MEDIUM / CONTEXT / LOW / NO_EVIDENCE',
    value:
      'HIGH = top-20% expression score AND a driver alteration is present. MEDIUM = top-20% score without an ' +
      'alteration. CONTEXT = an alteration is present but the expression score isn\'t in the top 20%. LOW = ' +
      'neither. For loss-of-function-relevant genes, "top" is direction-aware — a truncating mutation or ' +
      'deletion is expected to show LOW expression, not high, so the ranking direction flips for that gene\'s ' +
      'tier assignment. NO_EVIDENCE means there was no interpretable expression measurement at all for this ' +
      'pair, distinct from a real measurement that simply scored low.',
  },
];

export const MethodologySidebar: React.FC<MethodologySidebarProps> = ({
  isOpen,
  onClose,
  geneName,
  modelId,
  cellLineName,
}) => {
  const [expanded, setExpanded] = useState<number | null>(null);

  useEffect(() => {
    if (isOpen) document.body.style.overflow = 'hidden';
    else document.body.style.overflow = '';
    return () => {
      document.body.style.overflow = '';
    };
  }, [isOpen]);

  return (
    <>
      <div
        onClick={onClose}
        className={`fixed inset-0 bg-black/30 z-[999] transition-opacity duration-300 ${
          isOpen ? 'opacity-100 pointer-events-auto' : 'opacity-0 pointer-events-none'
        }`}
      />
      <aside
        className={`fixed top-0 right-0 h-screen w-screen sm:w-[420px] xl:w-[480px] max-w-[90vw] bg-white shadow-[-4px_0_24px_rgba(0,0,0,0.12)] z-[1000] overflow-y-auto p-6 transition-transform duration-300 ease-[cubic-bezier(0.16,1,0.3,1)] ${
          isOpen ? 'translate-x-0' : 'translate-x-full'
        }`}
      >
        <div className="flex items-start justify-between gap-4">
          <h3 className="text-lg font-bold text-slate-900 font-display">Scoring Methodology</h3>
          <button
            onClick={onClose}
            aria-label="Close"
            className="text-slate-400 hover:text-slate-700 transition shrink-0"
          >
            <X className="w-5 h-5" />
          </button>
        </div>
        <p className="text-xs text-slate-500 font-mono mb-1">
          {geneName} · {cellLineName} ({modelId})
        </p>
        <p className="text-xs text-slate-400 mb-6">
          General walkthrough of how every score is built — not this specific pair's own numbers.
        </p>

        <div>
          {METHODOLOGY_STEPS.map((step, i) => (
            <StepNode
              key={`${step.key}-${i}`}
              step={step}
              index={i}
              total={METHODOLOGY_STEPS.length}
              isExpanded={expanded === i}
              onToggle={() => setExpanded(expanded === i ? null : i)}
            />
          ))}
        </div>
      </aside>
    </>
  );
};

const StepNode: React.FC<{
  step: MethodologyStep;
  index: number;
  total: number;
  isExpanded: boolean;
  onToggle: () => void;
}> = ({ step, index, total, isExpanded, onToggle }) => {
  const isSkipped = step.status === 'skipped';
  const isInverted = step.status === 'inverted';

  const borderClass = isSkipped ? 'border-slate-300 border-dashed' : isInverted ? 'border-amber-300' : 'border-slate-200';
  const hoverBgClass = isSkipped ? 'hover:bg-slate-50' : isInverted ? 'hover:bg-amber-50' : 'hover:bg-violet-50';
  const bgClass = isExpanded ? (isSkipped ? 'bg-slate-50' : isInverted ? 'bg-amber-50' : 'bg-violet-50') : 'bg-white';
  const textClass = isSkipped ? 'text-slate-400' : 'text-slate-800';
  const circleActiveBg = isInverted ? 'bg-amber-500' : 'bg-mulberry-600';

  return (
    <div>
      <div
        onClick={onToggle}
        className={`flex items-start gap-3 px-4 py-3.5 rounded-lg border cursor-pointer transition-all hover:shadow-sm ${borderClass} ${bgClass} ${hoverBgClass}`}
      >
        <div
          className={`w-7 h-7 rounded-full flex items-center justify-center text-[13px] font-bold shrink-0 transition-colors ${
            isExpanded ? `${circleActiveBg} text-white` : `bg-slate-100 ${textClass}`
          }`}
        >
          {index + 1}
        </div>
        <div className="flex-1 min-w-0">
          <div className={`font-semibold text-sm flex items-center gap-2 flex-wrap ${textClass}`}>
            {step.key}
            {isSkipped && (
              <span className="text-[11px] font-normal text-slate-400 bg-slate-100 px-2 py-0.5 rounded">
                skipped
              </span>
            )}
            {isInverted && (
              <span className="text-[11px] font-normal text-amber-800 bg-amber-100 px-2 py-0.5 rounded">
                inverted
              </span>
            )}
          </div>
          {step.formula && <div className="font-mono text-xs text-slate-500 mt-1">{step.formula}</div>}
        </div>
        <div className={`text-xs text-slate-400 shrink-0 transition-transform ${isExpanded ? 'rotate-180' : ''}`}>
          ▼
        </div>
      </div>

      {isExpanded && step.value && (
        <div className="pl-[52px] pr-2 py-3 text-[13px] leading-relaxed text-slate-600 animate-fade-in">
          {step.value}
        </div>
      )}

      {index < total - 1 && (
        <div className={`w-0.5 h-4 ml-[42px] ${isSkipped ? 'bg-slate-200' : 'bg-mulberry-200'}`} />
      )}
    </div>
  );
};
