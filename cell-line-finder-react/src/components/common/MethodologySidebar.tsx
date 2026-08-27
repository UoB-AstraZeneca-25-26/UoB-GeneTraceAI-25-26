import React, { useEffect, useState } from 'react';
import { X, AlertTriangle } from 'lucide-react';
import { AGENT_API_URL } from '../../config';
import { MethodologyStep } from '../../types';

interface MethodologySidebarProps {
  isOpen: boolean;
  onClose: () => void;
  geneName: string;
  ensgId: string;
  modelId: string;
  cellLineName: string;
}

export const MethodologySidebar: React.FC<MethodologySidebarProps> = ({
  isOpen,
  onClose,
  geneName,
  ensgId,
  modelId,
  cellLineName,
}) => {
  const [steps, setSteps] = useState<MethodologyStep[] | null>(null);
  const [fallbackAnswer, setFallbackAnswer] = useState<string | null>(null);
  const [loading, setLoading] = useState(false);
  const [error, setError] = useState<string | null>(null);
  const [expanded, setExpanded] = useState<number | null>(null);
  const [attempt, setAttempt] = useState(0);

  useEffect(() => {
    if (!isOpen) return;
    let cancelled = false;

    setLoading(true);
    setError(null);
    setSteps(null);
    setFallbackAnswer(null);
    setExpanded(null);

    fetch(`${AGENT_API_URL}/v1/agent/query`, {
      method: 'POST',
      headers: { 'Content-Type': 'application/json' },
      body: JSON.stringify({ query: `Explain the score for ${ensgId} in ${modelId}` }),
    })
      .then((res) => res.json())
      .then((data) => {
        if (cancelled) return;
        if (data?.methodology?.steps) {
          setSteps(data.methodology.steps);
        } else if (data?.answer) {
          setFallbackAnswer(data.answer);
        } else {
          setError('Agent unavailable. Try again later.');
        }
      })
      .catch((err) => {
        console.error('Agent call failed:', err);
        if (!cancelled) setError('Agent unavailable. Try again later.');
      })
      .finally(() => {
        if (!cancelled) setLoading(false);
      });

    return () => {
      cancelled = true;
    };
  }, [isOpen, ensgId, modelId, attempt]);

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
        <p className="text-xs text-slate-500 font-mono mb-6">
          {geneName} · {cellLineName} ({modelId})
        </p>

        {loading && (
          <div className="space-y-2">
            {[0, 1, 2, 3, 4, 5, 6].map((i) => (
              <div
                key={i}
                className="h-16 rounded-lg bg-gradient-to-r from-slate-100 via-slate-200 to-slate-100 animate-shimmer"
                style={{ animationDelay: `${i * 0.1}s` }}
              />
            ))}
          </div>
        )}

        {!loading && error && (
          <div className="py-10 text-center text-slate-500">
            <AlertTriangle className="w-8 h-8 mx-auto mb-2 text-amber-500" />
            <p className="text-sm">{error}</p>
            <button
              onClick={() => setAttempt((a) => a + 1)}
              className="mt-3 px-4 py-1.5 border border-slate-300 rounded-lg text-xs bg-white hover:bg-slate-50 transition"
            >
              Retry
            </button>
          </div>
        )}

        {!loading && !error && steps && (
          <div>
            {steps.map((step, i) => (
              <StepNode
                key={`${step.key}-${i}`}
                step={step}
                index={i}
                total={steps.length}
                isExpanded={expanded === i}
                onToggle={() => setExpanded(expanded === i ? null : i)}
              />
            ))}
          </div>
        )}

        {!loading && !error && !steps && fallbackAnswer && (
          <div className="text-sm leading-relaxed text-slate-700 whitespace-pre-wrap">{fallbackAnswer}</div>
        )}
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
  const circleActiveBg = isInverted ? 'bg-amber-500' : 'bg-indigo-600';

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
        <div className={`w-0.5 h-4 ml-[42px] ${isSkipped ? 'bg-slate-200' : 'bg-indigo-200'}`} />
      )}
    </div>
  );
};
