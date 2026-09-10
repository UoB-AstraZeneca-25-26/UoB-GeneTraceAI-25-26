import React, { useEffect, useState } from 'react';
import { CellLineDetail, InspectTarget, RankedCellLine } from '../../types';
import { CellLineDetailApiResponse, fetchCellLineDetail, toCellLineDetail } from '../../lib/api';
import {
  confTier,
  pct,
  levelBand,
  TIER_PILL_CLASS,
  TIER_BAR_CLASS,
  EXPRESSION_SOURCES,
  PROTEOMICS_SOURCES,
} from '../../lib/evidence';
import { TierPill } from '../common/TierPill';
import { MethodologySidebar } from '../common/MethodologySidebar';
import { ArrowLeft, AlertCircle, Loader2, Sparkles, Zap, ChevronRight, GitBranch } from 'lucide-react';

interface ProfileViewProps {
  target: InspectTarget;
  onBack: () => void;
  onInspect: (t: InspectTarget) => void; // re-inspect an RNA alternative
  detailCache?: Map<string, CellLineDetailApiResponse>;
}

export const ProfileView: React.FC<ProfileViewProps> = ({ target, onBack, onInspect, detailCache }) => {
  const [detail, setDetail] = useState<CellLineDetail | null>(null);
  const [loading, setLoading] = useState<boolean>(true);
  const [error, setError] = useState<string | null>(null);
  const [showMethodology, setShowMethodology] = useState(false);

  // When the query had multiple TARGET genes (multi / jointSelectivity), evidence
  // (expression, proteomics, mutation, fusion, CNA) is only ever fetched for one
  // gene at a time -- /gene/detail takes a single gene param. Let the user switch
  // which target's evidence they're looking at, defaulting to the limiting gene
  // since that's usually the one that matters most for interpreting the score.
  const availableGenes = React.useMemo(() => {
    const r = target.ranked;
    if (r && (r.mode === 'multi' || r.mode === 'jointSelectivity') && r.genes.length > 1) return r.genes;
    return [target.gene];
  }, [target]);
  const limitingGene =
    target.ranked && (target.ranked.mode === 'multi' || target.ranked.mode === 'jointSelectivity')
      ? target.ranked.limitingGene
      : null;

  const [selectedGene, setSelectedGene] = useState<string>(limitingGene ?? availableGenes[0] ?? target.gene);

  // Reset which gene is selected whenever we navigate to a new inspect target
  // (different line, or a fresh query) -- otherwise a stale selectedGene from
  // the previous line could persist and silently fetch the wrong gene here.
  useEffect(() => {
    setSelectedGene(limitingGene ?? availableGenes[0] ?? target.gene);
    // The RNA-similar alternatives panel that leads here sits near the
    // bottom of a long profile page -- without this, clicking an alternative
    // swaps in the new line's data while the page stays scrolled down where
    // the click happened, so the header/tier/back-button above are off
    // screen. The actual scrolling element is <main> in App.tsx (it has
    // overflow-y-auto), not window, so window.scrollTo would silently do
    // nothing here.
    document.querySelector('main')?.scrollTo({ top: 0, behavior: 'auto' });
    // eslint-disable-next-line react-hooks/exhaustive-deps
  }, [target.modelId, target.gene]);

  // The rest of the original query's target/exclusion genes (if any), so
  // /gene/detail can scope rna_alternatives to lines that actually satisfy
  // the FULL query -- not just the one gene currently on screen.
  const otherGenes = React.useMemo(() => {
    const r = target.ranked;
    if (r && (r.mode === 'multi' || r.mode === 'jointSelectivity')) {
      return r.genes.filter((g) => g !== selectedGene);
    }
    return [];
  }, [target, selectedGene]);
  const excludeGenes = React.useMemo(() => {
    const r = target.ranked;
    if (r && (r.mode === 'selectivity' || r.mode === 'jointSelectivity')) return r.excludedGenes;
    return [];
  }, [target]);

  useEffect(() => {
    let cancelled = false;
    const ctrl = new AbortController();
    setLoading(true);
    setError(null);
    setDetail(null);
    fetchCellLineDetail(selectedGene, target.modelId, ctrl.signal, detailCache, otherGenes, excludeGenes)
      .then((resp) => {
        if (!cancelled) setDetail(toCellLineDetail(resp));
      })
      .catch((e: Error) => {
        if (!cancelled && e.name !== 'AbortError') setError(e.message);
      })
      .finally(() => {
        if (!cancelled) setLoading(false);
      });
    return () => {
      cancelled = true;
      ctrl.abort();
    };
    // eslint-disable-next-line react-hooks/exhaustive-deps
  }, [selectedGene, target.modelId, otherGenes.join(','), excludeGenes.join(',')]);

  const BackButton = (
    <button
      onClick={onBack}
      className="inline-flex items-center text-sm font-semibold text-slate-600 hover:text-slate-900 transition mb-3"
    >
      <ArrowLeft className="w-4 h-4 mr-1.5" />
      <span>Back to Results</span>
    </button>
  );

  if (loading) {
    return (
      <div className="max-w-5xl mx-auto">
        {BackButton}
        <div className="text-center py-24">
          <Loader2 className="w-8 h-8 mx-auto animate-spin text-mulberry-600" />
          <p className="text-sm text-slate-500 mt-3">
            Loading <strong className="font-mono">{target.cellLine}</strong> · {selectedGene}…
          </p>
        </div>
      </div>
    );
  }

  if (error || !detail) {
    return (
      <div className="max-w-5xl mx-auto">
        {BackButton}
        <div className="text-center py-16 space-y-3">
          <div className="inline-flex p-3 bg-rose-50 text-rose-600 rounded-full">
            <AlertCircle className="w-8 h-8" />
          </div>
          <h2 className="text-xl font-bold text-slate-800">Couldn't load this profile</h2>
          <p className="text-sm text-slate-500 max-w-md mx-auto">{error ?? 'No detail returned.'}</p>
        </div>
      </div>
    );
  }

  return (
    <>
    <div className="space-y-6 max-w-5xl mx-auto">
      <div>
        {BackButton}
        {target.viaAlternative && (
          <div
            className={`flex items-center gap-2 text-xs font-semibold px-3 py-2 rounded-lg border mb-3 ${TIER_PILL_CLASS.context}`}
          >
            <GitBranch className="w-3.5 h-3.5 shrink-0" />
            <span>Viewing an RNA-similar alternative — not part of the original ranked results.</span>
          </div>
        )}
        <div
          className={target.viaAlternative ? `border-l-4 pl-4 -ml-4 ${TIER_BAR_CLASS.context.replace('bg-', 'border-')}` : ''}
        >
          <div className="flex items-start justify-between gap-4">
            <div>
              <h2 className="text-3xl font-extrabold font-mono text-slate-900 tracking-tight">{detail.modelId}</h2>
              <p className="text-sm font-semibold text-slate-600">{detail.cellLine}</p>
              <GeneScores ranked={target.ranked} gene={detail.gene} score={detail.score} />
            </div>
            <TierPill tier={detail.tier} className="mt-2" />
          </div>
        </div>
      </div>

      {/* Stat cards — driven by the query the line was clicked from (falls back
          to single-gene detail when opened outside a ranking) */}
      <TopCards detail={detail} ranked={target.ranked} />

      <div className="flex items-center justify-between gap-4 flex-wrap">
        {availableGenes.length > 1 ? (
          <div className="flex items-center gap-2 flex-wrap">
            <span className="text-xs font-medium text-slate-500">Evidence below is for:</span>
            {availableGenes.map((g) => (
              <button
                key={g}
                onClick={() => setSelectedGene(g)}
                className={`inline-flex items-center gap-1 text-xs font-semibold px-3 py-1.5 rounded-full border transition ${
                  g === selectedGene
                    ? 'bg-mulberry-600 text-white border-mulberry-600'
                    : 'bg-white text-slate-600 border-slate-200 hover:border-mulberry-300 hover:text-mulberry-700'
                }`}
              >
                {g}
                {g === limitingGene && (
                  <span className={g === selectedGene ? 'text-[9px] text-white/80' : 'text-[9px] text-amber-600'}>▼ limiting</span>
                )}
              </button>
            ))}
          </div>
        ) : (
          <span />
        )}
        <button
          onClick={() => setShowMethodology(true)}
          className="inline-flex items-center gap-2 px-5 py-2.5 rounded-lg border border-mulberry-600 text-mulberry-600 text-sm font-semibold transition-colors hover:bg-mulberry-600 hover:text-white"
        >
          <Sparkles className="w-4 h-4" />
          Explain scoring methodology
        </button>
      </div>

      {/* Expression & proteomics — which datasets measured this line (lead section) */}
      <div className="bg-white rounded-xl shadow-sm border border-slate-200 p-6 space-y-4">
        <div className="flex items-center justify-between">
          <h3 className="text-sm font-bold text-slate-800">Expression &amp; proteomics</h3>
          <span className="text-[11px] text-slate-400">data sources</span>
        </div>

        <p className="text-[12px] text-slate-600 leading-relaxed bg-slate-50 border border-slate-100 rounded-lg p-3">
          <strong className="text-slate-700">How to read:</strong> each layer shows a level —{' '}
          <strong className="text-slate-700">Low, Moderate or High</strong> for {detail.gene} in this line{' '}
          relative to other cell lines — the same level is also shown per dataset below (Low/Moderate/High, or No
          if that dataset didn't measure it). Expression draws on DepMap, HPA and GEO; proteomics on ProCan and
          CCLE. Sources can disagree — that's real, not an error.
        </p>

        <div className="grid grid-cols-1 md:grid-cols-2 gap-4">
          <SourceRow
            label="Expression (RNA)"
            level={detail.expressionLevel}
            allSources={EXPRESSION_SOURCES}
            sourceLevels={detail.expressionBySource}
          />
          <SourceRow
            label="Proteomics (protein)"
            level={detail.proteomicsLevel}
            allSources={PROTEOMICS_SOURCES}
            sourceLevels={detail.proteomicsBySource}
          />
        </div>
      </div>

      {/* Alteration evidence — categorical layers */}
      <div className="bg-white rounded-xl shadow-sm border border-slate-200 p-6 space-y-4">
        <h3 className="text-sm font-bold text-slate-800">Alteration evidence</h3>
        <div className="grid grid-cols-1 sm:grid-cols-3 gap-3">
          <LayerCard
            label="Mutation"
            present={(detail.pMutation ?? 0) > 0}
            hint={detail.mutationDriver ? 'driver alteration' : (detail.pMutation ?? 0) > 0 ? 'variant detected' : 'no variant detected'}
            badge={detail.mutationDriver ? 'driver' : undefined}
          />
          <LayerCard
            label="Fusion"
            present={(detail.pFusion ?? 0) > 0}
            hint={detail.fusionDriver ? 'driver alteration' : (detail.pFusion ?? 0) > 0 ? 'fusion detected' : 'no fusion detected'}
            badge={detail.fusionDriver ? 'driver' : undefined}
          />
          <LayerCard
            label="Copy number"
            present={detail.hasCnaAlteration}
            hint={detail.hasCnaAlteration ? 'copy-number altered' : 'no alteration'}
            badge={detail.hasCnaAlteration ? 'driver' : undefined}
          />
        </div>
        <p className="text-[11px] text-slate-400">
          Alteration layers are categorical — each is simply present or absent for this line. A “driver”
          tag on a layer means that specific alteration (mutation, fusion, or copy-number change) is
          strong enough to be treated as a likely driver, not just any variant.
        </p>
      </div>

      {/* Cell-line metadata */}
      {detail.metadata.length > 0 && (
        <div className="bg-white rounded-xl shadow-sm border border-slate-200 p-6">
          <h3 className="text-sm font-bold text-slate-800 mb-4">Cell-line metadata</h3>
          <dl className="grid grid-cols-2 md:grid-cols-3 gap-x-6 gap-y-3">
            {detail.metadata.map((m) => (
              <div key={m.label}>
                <dt className="text-[11px] uppercase tracking-wide text-slate-400">{m.label}</dt>
                <dd className="text-sm text-slate-800 font-medium">{m.value}</dd>
              </div>
            ))}
          </dl>
        </div>
      )}

      {/* RNA-similar alternatives */}
      {detail.alternatives.length > 0 && (
        <div className="bg-white rounded-xl shadow-sm border border-slate-200 overflow-hidden">
          <div className="px-6 py-4 border-b border-slate-100">
            <h3 className="text-sm font-bold text-slate-800">RNA-similar alternatives</h3>
            <p className="text-xs text-slate-500 mt-0.5">
              Backup models with the closest transcriptomic profile. Click to inspect the same gene in that line.
            </p>
          </div>
          <ul className="divide-y divide-slate-100">
            {detail.alternatives.map((a) => {
              const p = pct(a.similarity);
              const simLabel = p >= 90 ? 'Very similar' : p >= 75 ? 'Similar' : 'Somewhat similar';
              const meta = [a.lineage, a.primaryDisease].filter(Boolean).join(' · ');
              return (
                <li key={a.modelId}>
                  <button
                    onClick={() =>
                      onInspect({
                        gene: detail.gene,
                        modelId: a.modelId,
                        cellLine: a.name,
                        ranked: target.ranked,
                        viaAlternative: true,
                      })
                    }
                    className="w-full flex items-center justify-between gap-4 px-6 py-3.5 hover:bg-slate-50 transition text-left"
                  >
                    <div className="min-w-0">
                      <div className="flex items-center gap-2">
                        <span className="font-bold font-mono text-slate-900">{a.modelId}</span>
                        <span className="text-[11px] text-slate-400">{a.name}</span>
                      </div>
                      {meta && <span className="block text-xs text-slate-500 capitalize mt-0.5 truncate">{meta}</span>}
                    </div>
                    <div className="flex items-center gap-3 shrink-0">
                      <div className="text-right">
                        <span className="block text-xs font-semibold text-mulberry-700">{simLabel}</span>
                        <span className="block text-[11px] text-slate-400 font-mono">{p}% RNA match</span>
                      </div>
                      <span className="inline-flex items-center text-xs font-semibold text-mulberry-600">
                        Inspect <ChevronRight className="w-4 h-4 ml-0.5" />
                      </span>
                    </div>
                  </button>
                </li>
              );
            })}
          </ul>
        </div>
      )}
    </div>

    <MethodologySidebar
      isOpen={showMethodology}
      onClose={() => setShowMethodology(false)}
      geneName={detail.gene}
      ensgId={detail.ensg}
      modelId={detail.modelId}
      cellLineName={detail.cellLine}
    />
    </>
  );
};

/** Context-aware gene score chips shown in the profile header — reflects the
 *  query the line was clicked from (single / joint / selectivity). */
const GeneScores: React.FC<{ ranked?: RankedCellLine; gene: string; score: number }> = ({ ranked, gene, score }) => {
  const chip = (label: string, tone: 'neutral' | 'high' | 'low' | 'limit', key: string) => {
    const cls =
      tone === 'high'
        ? 'bg-emerald-50 text-emerald-700 border-emerald-200'
        : tone === 'low'
        ? 'bg-rose-50 text-rose-700 border-rose-200'
        : tone === 'limit'
        ? 'bg-amber-50 text-amber-800 border-amber-200'
        : 'bg-slate-50 text-slate-700 border-slate-200';
    return (
      <span key={key} className={`inline-flex items-center text-[11px] font-mono px-2 py-0.5 rounded border ${cls}`}>
        {label}
      </span>
    );
  };

  let chips: React.ReactNode;
  if (!ranked || ranked.mode === 'single') {
    const v = ranked && ranked.mode === 'single' ? ranked.score : score;
    chips = chip(`${gene} ${v.toFixed(3)}`, 'neutral', gene);
  } else if (ranked.mode === 'multi') {
    chips = ranked.genes.map((g) => {
      const v = ranked.geneScores[g];
      const lim = ranked.limitingGene === g;
      return chip(`${g} ${v != null ? v.toFixed(3) : '—'}${lim ? ' ▼' : ''}`, lim ? 'limit' : 'neutral', g);
    });
  } else if (ranked.mode === 'jointSelectivity') {
    chips = [
      ...ranked.genes.map((g) => {
        const v = ranked.geneScores[g];
        const lim = ranked.limitingGene === g;
        return chip(`${g} ${v != null ? v.toFixed(3) : '—'}${lim ? ' ▼' : ''}`, lim ? 'limit' : 'high', g);
      }),
      ...ranked.excludedGenes.map((g) => chip(`${g} ${(ranked.exclusionScores[g] ?? 0).toFixed(3)} ↓`, 'low', g)),
    ];
  } else {
    chips = [
      chip(`${ranked.geneHigh} ${ranked.scoreHigh.toFixed(3)} ↑`, 'high', ranked.geneHigh),
      ...ranked.excludedGenes.map((g) => chip(`${g} ${(ranked.exclusionScores[g] ?? 0).toFixed(3)} ↓`, 'low', g)),
    ];
  }

  return <div className="flex flex-wrap gap-1.5 mt-2">{chips}</div>;
};

const TopCards: React.FC<{ detail: CellLineDetail; ranked?: RankedCellLine }> = ({ detail, ranked }) => {
  const lineageName = detail.metadata.find((m) => m.label === 'Lineage')?.value ?? ranked?.lineage ?? 'lineage';

  // Opened outside a ranking (RNA-similar alternative): single-gene detail is all we have.
  if (!ranked) {
    return (
      <div className="grid grid-cols-2 lg:grid-cols-4 gap-4">
        <StatCard label="Global rank" value={`#${detail.rank.toLocaleString()}`} sub={`of ${detail.total.toLocaleString()} for ${detail.gene}`} />
        <StatCard label="Global score" value={detail.score.toFixed(4)} sub={`${confTier(detail.tier).label} tier`} />
        <StatCard
          label="In-lineage"
          value={detail.lineageScore != null ? detail.lineageScore.toFixed(3) : '—'}
          sub={detail.lineageRank != null && detail.lineageTotal != null ? `#${detail.lineageRank} of ${detail.lineageTotal} in ${lineageName}` : 'within its lineage'}
        />
        <StatCard label="Evidence layers" value={String(detail.nLayers)} sub={detail.driverAlteration ? 'driver alteration present' : 'no driver alteration'} />
      </div>
    );
  }

  // From a ranking: show the score/rank the user actually clicked.
  const modeLabel =
    ranked.mode === 'selectivity'
      ? 'selectivity ranking'
      : ranked.mode === 'jointSelectivity'
      ? 'joint selectivity ranking'
      : ranked.mode === 'multi'
      ? 'joint ranking'
      : 'single-gene ranking';
  const metric =
    ranked.mode === 'selectivity'
      ? { label: 'Selectivity', value: ranked.selectivity }
      : ranked.mode === 'jointSelectivity'
      ? { label: 'Combined score', value: ranked.combinedScore }
      : ranked.mode === 'multi'
      ? { label: 'Joint score', value: ranked.jointScore }
      : { label: 'Score', value: ranked.score };
  const scoreSub =
    ranked.mode === 'selectivity'
      ? `${ranked.geneHigh} high vs ${ranked.excludedGenes.join('/')} low`
      : ranked.mode === 'jointSelectivity'
      ? `${ranked.genes.filter((g) => ranked.geneScores[g] != null).length}/${ranked.genes.length} targets scored, vs ${ranked.excludedGenes.join('/')} low`
      : ranked.mode === 'multi'
      ? `${ranked.genes.filter((g) => ranked.geneScores[g] != null).length}/${ranked.genes.length} genes scored`
      : `${confTier(detail.tier).label} tier`;

  return (
    <div className="grid grid-cols-2 lg:grid-cols-4 gap-4">
      <StatCard label="Rank" value={`#${ranked.rank}`} sub={`in this ${modeLabel}`} />
      <StatCard label={metric.label} value={metric.value.toFixed(4)} sub={scoreSub} />
      <StatCard label="In-lineage" value={ranked.lineageScore.toFixed(3)} sub={`in ${ranked.lineage}`} />
      <StatCard label="Evidence layers" value={String(detail.nLayers)} sub={detail.driverAlteration ? 'driver alteration present' : 'no driver alteration'} />
    </div>
  );
};

const StatCard: React.FC<{ label: string; value: string; sub?: string }> = ({ label, value, sub }) => (
  <div className="bg-white border border-slate-200 rounded-xl p-4">
    <span className="text-xs text-slate-500 block">{label}</span>
    <span className="text-2xl font-black text-slate-900 font-mono">{value}</span>
    {sub && <span className="text-[11px] text-slate-400 block mt-0.5">{sub}</span>}
  </div>
);

const LayerCard: React.FC<{ label: string; present: boolean; hint: string; badge?: string }> = ({
  label,
  present,
  hint,
  badge,
}) => (
  <div className={`rounded-lg border p-3 ${present ? 'border-slate-200 bg-white' : 'border-slate-100 bg-slate-50'}`}>
    <div className="flex items-center justify-between">
      <span className={`text-xs font-semibold ${present ? 'text-slate-700' : 'text-slate-400'}`}>{label}</span>
      {badge && (
        <span className="inline-flex items-center gap-0.5 text-[9px] font-bold text-purple-700 bg-purple-50 border border-purple-200 px-1.5 py-0.5 rounded">
          <Zap className="w-2.5 h-2.5" /> {badge}
        </span>
      )}
    </div>
    <div className="flex items-baseline gap-2 mt-1">
      <span className={`text-lg font-bold ${present ? 'text-emerald-600' : 'text-slate-400'}`}>
        {present ? 'Yes' : 'No'}
      </span>
      <span className="text-[11px] text-slate-400">{hint}</span>
    </div>
  </div>
);

const SourceRow: React.FC<{
  label: string;
  level: number | null;
  allSources: readonly string[];
  sourceLevels: Record<string, number | null>;
}> = ({ label, level, allSources, sourceLevels }) => {
  const measured = level != null;
  const band = measured ? levelBand(level) : null;
  return (
    <div className={`rounded-lg border p-4 space-y-3 ${measured ? 'border-slate-200 bg-white' : 'border-slate-100 bg-slate-50'}`}>
      <div className="flex items-center justify-between">
        <span className={`text-sm font-semibold ${measured ? 'text-slate-700' : 'text-slate-500'}`}>{label}</span>
        {band ? (
          <span className={`text-[11px] font-bold px-2.5 py-0.5 rounded-full border ${TIER_PILL_CLASS[band.tone]}`}>
            {band.label}
          </span>
        ) : (
          <span className="text-[11px] font-medium text-slate-400">Not measured</span>
        )}
      </div>
      <div>
        <span className="text-[10px] uppercase tracking-wide text-slate-400">Sources</span>
        <div className="flex flex-wrap gap-1.5 mt-1">
          {allSources.map((s) => {
            const srcLevel = sourceLevels[s];
            const srcBand = srcLevel != null ? levelBand(srcLevel) : null;
            return (
              <span
                key={s}
                title={srcBand ? `${s}: ${srcBand.label} for this line` : `${s}: not measured for this line`}
                className={`inline-flex items-center gap-1 text-[11px] font-medium px-2 py-0.5 rounded-full border ${
                  srcBand
                    ? TIER_PILL_CLASS[srcBand.tone]
                    : 'bg-slate-50 text-slate-400 border-dashed border-slate-200'
                }`}
              >
                {s}: {srcBand ? srcBand.label : 'No'}
              </span>
            );
          })}
        </div>
      </div>
    </div>
  );
};
