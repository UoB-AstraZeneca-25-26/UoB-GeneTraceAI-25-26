import React, { useEffect, useState } from 'react';
import { CellLineDetail, InspectTarget, SourceMeasurement } from '../../types';
import { fetchCellLineDetail, toCellLineDetail } from '../../lib/api';
import {
  confTier,
  pct,
  levelBand,
  TIER_PILL_CLASS,
  TRACKS,
  EXPRESSION_SOURCES,
  PROTEOMICS_SOURCES,
} from '../../lib/evidence';
import { TierPill } from '../common/TierPill';
import { EvidenceSpine } from '../common/EvidenceSpine';
import { MethodologySidebar } from '../common/MethodologySidebar';
import { ArrowLeft, AlertCircle, Loader2, Sparkles, Zap, ChevronRight } from 'lucide-react';

interface ProfileViewProps {
  target: InspectTarget;
  onBack: () => void;
  onInspect: (t: InspectTarget) => void; // re-inspect an RNA alternative
}

export const ProfileView: React.FC<ProfileViewProps> = ({ target, onBack, onInspect }) => {
  const [detail, setDetail] = useState<CellLineDetail | null>(null);
  const [loading, setLoading] = useState<boolean>(true);
  const [error, setError] = useState<string | null>(null);
  const [showMethodology, setShowMethodology] = useState(false);

  useEffect(() => {
    let cancelled = false;
    const ctrl = new AbortController();
    setLoading(true);
    setError(null);
    setDetail(null);
    fetchCellLineDetail(target.gene, target.modelId, ctrl.signal)
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
  }, [target.gene, target.modelId]);

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
          <Loader2 className="w-8 h-8 mx-auto animate-spin text-indigo-600" />
          <p className="text-sm text-slate-500 mt-3">
            Loading <strong className="font-mono">{target.cellLine}</strong> · {target.gene}…
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
        <div className="flex items-start justify-between gap-4">
          <div>
            <h2 className="text-3xl font-extrabold text-slate-900 font-display tracking-tight">{detail.cellLine}</h2>
            <p className="text-sm text-slate-500 font-mono">
              {detail.modelId} · {detail.gene}
              <span className="text-slate-400"> ({detail.ensg})</span>
            </p>
          </div>
          <TierPill tier={detail.tier} className="mt-2" />
        </div>
      </div>

      {/* Stat cards for THIS gene in THIS line */}
      <div className="grid grid-cols-1 md:grid-cols-3 gap-4">
        <StatCard label={`Rank for ${detail.gene}`} value={`#${detail.rank.toLocaleString()}`} sub={`of ${detail.total.toLocaleString()}`} />
        <StatCard label="Score" value={detail.score.toFixed(4)} sub={`${confTier(detail.tier).label} tier`} />
        <StatCard label="Evidence layers" value={String(detail.nLayers)} sub={detail.driverAlteration ? 'driver alteration present' : 'no driver alteration'} />
      </div>

      {/* Evidence fingerprint — labeled at-a-glance overview of all six layers */}
      <div className="bg-white rounded-xl shadow-sm border border-slate-200 p-5">
        <div className="flex flex-col sm:flex-row sm:items-start gap-3 sm:gap-6">
          <div className="shrink-0">
            <span className="text-[11px] uppercase tracking-wide font-semibold text-slate-500 block mb-1.5">
              Evidence fingerprint
            </span>
            <EvidenceSpine tracks={detail.tracks} />
          </div>
          <div className="flex flex-wrap gap-x-3 gap-y-1.5 text-[11px]">
            {TRACKS.map((t) => {
              const s = detail.tracks[t.key]?.score ?? 0;
              const on = s > 0;
              const readout = on ? (t.kind === 'call' ? 'present' : `${pct(s)}th pct`) : 'silent';
              return (
                <span key={t.key} className="inline-flex items-center gap-1.5">
                  <span className="w-2.5 h-2.5 rounded-sm" style={{ background: on ? t.color : '#e2e8f0' }} />
                  <span className={on ? 'text-slate-600' : 'text-slate-300'}>
                    {t.label} <span className="text-slate-400">· {readout}</span>
                  </span>
                </span>
              );
            })}
          </div>
        </div>
        <p className="text-[11px] text-slate-400 mt-3">
          One segment per evidence layer, in the order above: filled when the layer has signal for this line
          (brighter = stronger for expression &amp; proteomics; solid = present for mutation, fusion &amp; copy
          number), hollow when silent.
        </p>
      </div>

      <div className="flex justify-end">
        <button
          onClick={() => setShowMethodology(true)}
          className="inline-flex items-center gap-2 px-5 py-2.5 rounded-lg border border-indigo-600 text-indigo-600 text-sm font-semibold transition-colors hover:bg-indigo-600 hover:text-white"
        >
          <Sparkles className="w-4 h-4" />
          Explain scoring methodology
        </button>
      </div>

      {/* Expression & proteomics — measured levels, sources and raw values (lead section) */}
      <div className="bg-white rounded-xl shadow-sm border border-slate-200 p-6 space-y-4">
        <div className="flex items-center justify-between">
          <h3 className="text-sm font-bold text-slate-800">Expression &amp; proteomics</h3>
          <span className="text-[11px] text-slate-400">measured levels</span>
        </div>

        <p className="text-[12px] text-slate-600 leading-relaxed bg-slate-50 border border-slate-100 rounded-lg p-3">
          <strong className="text-slate-700">How to read:</strong> the band and meter are a{' '}
          <strong className="text-slate-700">within-gene percentile</strong> — where {detail.cellLine} sits
          among all cell lines for {detail.gene}. <strong className="text-slate-700">0 = lowest, 100 = highest.</strong>{' '}
          It’s a relative level, not an absolute amount or a probability: “High” means more of the transcript
          (expression) or protein (proteomics) than most other lines. Below each meter are the raw measurements
          from each dataset — bars are scaled per source, since their units differ.
        </p>

        <div className="grid grid-cols-1 md:grid-cols-2 gap-4">
          <LevelRow
            label="Expression (RNA)"
            level={detail.expressionLevel}
            allSources={EXPRESSION_SOURCES}
            presentSources={detail.expressionSources}
            measurements={detail.expressionMeasurements}
          />
          <LevelRow
            label="Proteomics (protein)"
            level={detail.proteomicsLevel}
            allSources={PROTEOMICS_SOURCES}
            presentSources={detail.proteomicsSources}
            measurements={detail.proteomicsMeasurements}
          />
        </div>
      </div>

      {/* Alteration evidence — categorical layers */}
      <div className="bg-white rounded-xl shadow-sm border border-slate-200 p-6 space-y-4">
        <h3 className="text-sm font-bold text-slate-800">Alteration evidence</h3>
        <div className="grid grid-cols-1 sm:grid-cols-3 gap-3">
          <LayerCard
            label="Mutation"
            present={detail.pMutation != null}
            hint={detail.driverAlteration ? 'driver alteration' : detail.pMutation != null ? 'variant detected' : 'no variant detected'}
            badge={detail.driverAlteration ? 'driver' : undefined}
          />
          <LayerCard
            label="Fusion"
            present={detail.pFusion != null}
            hint={detail.pFusion != null ? 'fusion detected' : 'no fusion detected'}
          />
          <LayerCard
            label="Copy number"
            present={detail.hasCnaAlteration}
            hint={detail.hasCnaAlteration ? 'copy-number altered' : 'no alteration'}
          />
        </div>
        <p className="text-[11px] text-slate-400">
          Alteration layers are categorical — each is simply present or absent for this line. A “driver”
          tag means the mutation is a known cancer driver, not just any variant.
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
            {detail.alternatives.map((a) => (
              <li key={a.modelId}>
                <button
                  onClick={() => onInspect({ gene: detail.gene, modelId: a.modelId, cellLine: a.name })}
                  className="w-full flex items-center justify-between px-6 py-3 hover:bg-slate-50 transition text-left"
                >
                  <div>
                    <span className="font-bold text-slate-900">{a.name}</span>
                    <span className="text-[11px] text-slate-400 font-mono ml-2">{a.modelId}</span>
                  </div>
                  <div className="flex items-center gap-3">
                    <div className="flex items-center gap-2">
                      <div className="w-24 bg-slate-200 rounded-full h-2 overflow-hidden">
                        <div className="h-full rounded-full bg-indigo-500" style={{ width: `${pct(a.similarity)}%` }} />
                      </div>
                      <span className="text-xs font-mono text-slate-500 w-10 text-right">{pct(a.similarity)}%</span>
                    </div>
                    <ChevronRight className="w-4 h-4 text-slate-400" />
                  </div>
                </button>
              </li>
            ))}
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

const LevelRow: React.FC<{
  label: string;
  level: number | null;
  allSources: readonly string[];
  presentSources: string[];
  measurements: SourceMeasurement[];
}> = ({ label, level, allSources, presentSources, measurements }) => {
  if (level == null) {
    return (
      <div className="rounded-lg border border-slate-100 bg-slate-50 p-4">
        <div className="flex items-center justify-between">
          <span className="text-sm font-semibold text-slate-500">{label}</span>
          <span className="text-[11px] font-medium text-slate-400">Not measured</span>
        </div>
        <p className="text-[11px] text-slate-400 mt-2">No {label.toLowerCase()} evidence reported for this line.</p>
      </div>
    );
  }

  const band = levelBand(level);
  const p = pct(level);

  return (
    <div className="rounded-lg border border-slate-200 bg-white p-4 space-y-3">
      <div className="flex items-center justify-between">
        <span className="text-sm font-semibold text-slate-700">{label}</span>
        <span className={`text-[11px] font-bold px-2 py-0.5 rounded-full border ${TIER_PILL_CLASS[band.tone]}`}>
          {band.label}
        </span>
      </div>

      {/* Low → High meter with a marker at the percentile */}
      <div>
        <div className="relative h-2.5 rounded-full bg-gradient-to-r from-slate-200 via-slate-300 to-indigo-500">
          <div
            className="absolute top-1/2 -translate-y-1/2 -translate-x-1/2 w-3 h-3 rounded-full bg-white border-2 border-indigo-600 shadow-sm"
            style={{ left: `${p}%` }}
            title={`${p}th percentile`}
          />
        </div>
        <div className="flex justify-between text-[10px] mt-1">
          <span className="text-slate-400">Low</span>
          <span className="font-mono font-semibold text-slate-600">{p}th percentile</span>
          <span className="text-slate-400">High</span>
        </div>
      </div>

      {/* Raw per-source measurements as a barplot; else which sources were present */}
      {measurements.length > 0 ? (
        <div>
          <span className="text-[10px] uppercase tracking-wide text-slate-400">Raw measurements by source</span>
          <div className="space-y-1.5 mt-1.5">
            {measurements.map((m) => {
              const w = Math.max(3, Math.min(100, (m.value / m.max) * 100));
              return (
                <div
                  key={m.source}
                  className="flex items-center gap-2"
                  title={`${m.source}: ${m.value} ${m.unit}`}
                >
                  <span className="w-16 shrink-0 text-[11px] font-medium text-slate-600">{m.source}</span>
                  <div className="flex-1 h-2 rounded-full bg-slate-100 overflow-hidden">
                    <div className="h-full rounded-full bg-indigo-500" style={{ width: `${w}%` }} />
                  </div>
                  <span className="w-28 shrink-0 text-right text-[10px] font-mono text-slate-500 tabular-nums">
                    {m.value} <span className="text-slate-400">{m.unit}</span>
                  </span>
                </div>
              );
            })}
          </div>
          <p className="text-[10px] text-slate-400 mt-1.5">Bars are scaled per source — units differ, so lengths aren’t comparable across rows.</p>
        </div>
      ) : (
        <div>
          <span className="text-[10px] uppercase tracking-wide text-slate-400">Sources</span>
          <div className="flex flex-wrap gap-1.5 mt-1">
            {allSources.map((s) => {
              const on = presentSources.includes(s);
              return (
                <span
                  key={s}
                  title={on ? `${s} contributed this level` : `${s} did not contribute for this line`}
                  className={`inline-flex items-center gap-1 text-[11px] font-medium px-2 py-0.5 rounded-full border ${
                    on
                      ? 'bg-indigo-50 text-indigo-700 border-indigo-200'
                      : 'bg-slate-50 text-slate-400 border-dashed border-slate-200'
                  }`}
                >
                  {on && <span aria-hidden>✓</span>}
                  {s}
                </span>
              );
            })}
          </div>
        </div>
      )}
    </div>
  );
};
