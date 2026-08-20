import React, { useEffect, useState } from 'react';
import { CellLineDetail, InspectTarget } from '../../types';
import { fetchCellLineDetail, toCellLineDetail } from '../../lib/api';
import { confTier, pct } from '../../lib/evidence';
import { TierPill } from '../common/TierPill';
import { EvidenceSpine } from '../common/EvidenceSpine';
import { ArrowLeft, AlertCircle, Loader2, Zap, ChevronRight } from 'lucide-react';

interface ProfileViewProps {
  target: InspectTarget;
  onBack: () => void;
  onInspect: (t: InspectTarget) => void; // re-inspect an RNA alternative
}

export const ProfileView: React.FC<ProfileViewProps> = ({ target, onBack, onInspect }) => {
  const [detail, setDetail] = useState<CellLineDetail | null>(null);
  const [loading, setLoading] = useState<boolean>(true);
  const [error, setError] = useState<string | null>(null);

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

      {/* Evidence layers reported by this endpoint */}
      <div className="bg-white rounded-xl shadow-sm border border-slate-200 p-6 space-y-4">
        <div className="flex items-center justify-between">
          <h3 className="text-sm font-bold text-slate-800">Alteration evidence</h3>
          <EvidenceSpine tracks={detail.tracks} />
        </div>
        <div className="grid grid-cols-1 sm:grid-cols-3 gap-3">
          <LayerCard
            label="Mutation"
            value={detail.pMutation != null ? `p = ${detail.pMutation.toFixed(3)}` : 'not reported'}
            active={detail.pMutation != null}
            badge={detail.driverAlteration ? 'driver' : undefined}
          />
          <LayerCard
            label="Fusion"
            value={detail.pFusion != null ? `p = ${detail.pFusion.toFixed(3)}` : 'none detected'}
            active={detail.pFusion != null}
          />
          <LayerCard
            label="Copy number"
            value={detail.hasCnaAlteration ? 'altered' : 'no alteration'}
            active={detail.hasCnaAlteration}
          />
        </div>
        <p className="text-[11px] text-slate-400">
          These are the alteration layers this endpoint reports. Expression and proteomics contribute
          to the score but aren't broken out per-layer here.
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
  );
};

const StatCard: React.FC<{ label: string; value: string; sub?: string }> = ({ label, value, sub }) => (
  <div className="bg-white border border-slate-200 rounded-xl p-4">
    <span className="text-xs text-slate-500 block">{label}</span>
    <span className="text-2xl font-black text-slate-900 font-mono">{value}</span>
    {sub && <span className="text-[11px] text-slate-400 block mt-0.5">{sub}</span>}
  </div>
);

const LayerCard: React.FC<{ label: string; value: string; active: boolean; badge?: string }> = ({
  label,
  value,
  active,
  badge,
}) => (
  <div className={`rounded-lg border p-3 ${active ? 'border-slate-200 bg-white' : 'border-slate-100 bg-slate-50'}`}>
    <div className="flex items-center justify-between">
      <span className={`text-xs font-semibold ${active ? 'text-slate-700' : 'text-slate-400'}`}>{label}</span>
      {badge && (
        <span className="inline-flex items-center gap-0.5 text-[9px] font-bold text-purple-700 bg-purple-50 border border-purple-200 px-1.5 py-0.5 rounded">
          <Zap className="w-2.5 h-2.5" /> {badge}
        </span>
      )}
    </div>
    <span className={`text-sm font-mono mt-1 block ${active ? 'text-slate-900' : 'text-slate-400'}`}>{value}</span>
  </div>
);
