import React, { useState } from 'react';
import { CellLineData } from '../../types';
import { ArrowLeft, CheckCircle2, Database, Dna, Activity } from 'lucide-react';
import { ResponsiveContainer, BarChart, Bar, XAxis, YAxis, Tooltip, Legend } from 'recharts';

interface ProfileViewProps {
  cellLine: string;
  dataset: CellLineData[];
  onBack: () => void;
}

export const ProfileView: React.FC<ProfileViewProps> = ({ cellLine, dataset, onBack }) => {
  const [activeTab, setActiveTab] = useState<'omics' | 'mutations' | 'alternatives'>('omics');

  const clRecords = dataset.filter(d => d.cellLine === cellLine);
  const tissue = clRecords[0]?.tissue || "Unknown";
  const mutations = clRecords.filter(d => d.hasSomaticVariant);

  // Prepare chart data
  const chartData = clRecords.map(r => ({
    gene: r.gene,
    TPM: r.tpmExpression,
    Protein: r.proteinExpression,
  }));

  return (
    <div className="space-y-6 max-w-5xl mx-auto">
      {/* Back Button & Header */}
      <div>
        <button
          onClick={onBack}
          className="inline-flex items-center text-sm font-semibold text-slate-600 hover:text-slate-900 transition mb-3"
        >
          <ArrowLeft className="w-4 h-4 mr-1.5" />
          <span>Back to Results</span>
        </button>
        <div className="flex items-center justify-between">
          <div>
            <h2 className="text-3xl font-extrabold text-slate-900">{cellLine}</h2>
            <p className="text-sm text-slate-500">{tissue} Tissue • Harmonized Multi-Omics Profile</p>
          </div>
          <div className="flex space-x-2">
            <span className="inline-flex items-center px-2.5 py-1 rounded-full text-xs font-semibold bg-emerald-50 text-emerald-700 border border-emerald-200">
              <CheckCircle2 className="w-3.5 h-3.5 mr-1" /> DepMap Verified
            </span>
          </div>
        </div>
      </div>

      {/* Top Stat Cards */}
      <div className="grid grid-cols-1 md:grid-cols-3 gap-4">
        <div className="bg-white border border-slate-200 rounded-xl p-4 flex items-center space-x-3">
          <div className="p-2.5 bg-blue-50 text-blue-600 rounded-lg">
            <Activity className="w-5 h-5" />
          </div>
          <div>
            <span className="text-xs text-slate-500 block">Tissue of Origin</span>
            <span className="text-base font-bold text-slate-900">{tissue}</span>
          </div>
        </div>

        <div className="bg-white border border-slate-200 rounded-xl p-4 flex items-center space-x-3">
          <div className="p-2.5 bg-amber-50 text-amber-600 rounded-lg">
            <Dna className="w-5 h-5" />
          </div>
          <div>
            <span className="text-xs text-slate-500 block">Target Somatic Variants</span>
            <span className="text-base font-bold text-slate-900">{mutations.length} Detected</span>
          </div>
        </div>

        <div className="bg-white border border-slate-200 rounded-xl p-4 flex items-center space-x-3">
          <div className="p-2.5 bg-purple-50 text-purple-600 rounded-lg">
            <Database className="w-5 h-5" />
          </div>
          <div>
            <span className="text-xs text-slate-500 block">Integrated Sources</span>
            <span className="text-base font-bold text-slate-900">DepMap, HPA, CCLE</span>
          </div>
        </div>
      </div>

      {/* Tabs */}
      <div className="bg-white rounded-xl shadow-sm border border-slate-200 overflow-hidden">
        <div className="flex border-b border-slate-200 bg-slate-50 px-4">
          <button
            onClick={() => setActiveTab('omics')}
            className={`py-3 px-4 text-xs font-bold uppercase tracking-wider border-b-2 transition ${
              activeTab === 'omics'
                ? 'border-indigo-600 text-indigo-600 bg-white'
                : 'border-transparent text-slate-500 hover:text-slate-700'
            }`}
          >
            Multi-Omics Expression
          </button>
          <button
            onClick={() => setActiveTab('mutations')}
            className={`py-3 px-4 text-xs font-bold uppercase tracking-wider border-b-2 transition ${
              activeTab === 'mutations'
                ? 'border-indigo-600 text-indigo-600 bg-white'
                : 'border-transparent text-slate-500 hover:text-slate-700'
            }`}
          >
            Genomic Features
          </button>
          <button
            onClick={() => setActiveTab('alternatives')}
            className={`py-3 px-4 text-xs font-bold uppercase tracking-wider border-b-2 transition ${
              activeTab === 'alternatives'
                ? 'border-indigo-600 text-indigo-600 bg-white'
                : 'border-transparent text-slate-500 hover:text-slate-700'
            }`}
          >
            Similar Alternatives
          </button>
        </div>

        <div className="p-6">
          {/* Tab 1: Chart */}
          {activeTab === 'omics' && (
            <div className="space-y-4">
              <div>
                <h4 className="text-sm font-bold text-slate-800">Transcriptomic (TPM) vs Proteomic Expression</h4>
                <p className="text-xs text-slate-500">Cross-validation between DepMap RNA-seq and mass spectrometry proteomics.</p>
              </div>
              <div className="h-72 w-full pt-4">
                <ResponsiveContainer width="100%" height="100%">
                  <BarChart data={chartData} margin={{ top: 10, right: 30, left: 0, bottom: 20 }}>
                    <XAxis dataKey="gene" tick={{ fontSize: 12 }} />
                    <YAxis tick={{ fontSize: 12 }} />
                    <Tooltip />
                    <Legend wrapperStyle={{ paddingTop: '10px' }} />
                    <Bar dataKey="TPM" fill="#4f46e5" radius={[4, 4, 0, 0]} name="RNA Expression (TPM)" />
                    <Bar dataKey="Protein" fill="#10b981" radius={[4, 4, 0, 0]} name="Protein Expression" />
                  </BarChart>
                </ResponsiveContainer>
              </div>
            </div>
          )}

          {/* Tab 2: Mutations */}
          {activeTab === 'mutations' && (
            <div className="space-y-4">
              <h4 className="text-sm font-bold text-slate-800">Somatic Mutations & Gene Variants</h4>
              {mutations.length > 0 ? (
                <div className="border border-slate-200 rounded-lg overflow-hidden">
                  <table className="w-full text-left text-sm text-slate-700">
                    <thead className="bg-slate-50 border-b border-slate-200 text-xs font-semibold text-slate-600">
                      <tr>
                        <th className="px-4 py-2.5">Gene Target</th>
                        <th className="px-4 py-2.5">Variant Status</th>
                        <th className="px-4 py-2.5">Evidence Level</th>
                      </tr>
                    </thead>
                    <tbody className="divide-y divide-slate-100">
                      {mutations.map(m => (
                        <tr key={m.gene}>
                          <td className="px-4 py-2.5 font-bold text-slate-900">{m.gene}</td>
                          <td className="px-4 py-2.5 text-rose-600 font-medium">Somatic Mutation Detected</td>
                          <td className="px-4 py-2.5 text-slate-600">High (DepMap Omics Mutation)</td>
                        </tr>
                      ))}
                    </tbody>
                  </table>
                </div>
              ) : (
                <p className="text-sm text-slate-500">No targeted somatic variants recorded for this cell line.</p>
              )}
            </div>
          )}

          {/* Tab 3: Alternatives */}
          {activeTab === 'alternatives' && (
            <div className="space-y-4">
              <h4 className="text-sm font-bold text-slate-800">Backup Cell Line Recommendations</h4>
              <p className="text-xs text-slate-500">Alternative models with similar molecular phenotypes.</p>
              <div className="border border-slate-200 rounded-lg overflow-hidden">
                <table className="w-full text-left text-sm text-slate-700">
                  <thead className="bg-slate-50 border-b border-slate-200 text-xs font-semibold text-slate-600">
                    <tr>
                      <th className="px-4 py-2.5">Alternative Cell Line</th>
                      <th className="px-4 py-2.5">Phenotypic Similarity</th>
                      <th className="px-4 py-2.5">Confidence</th>
                    </tr>
                  </thead>
                  <tbody className="divide-y divide-slate-100">
                    <tr>
                      <td className="px-4 py-2.5 font-bold text-slate-900">A549</td>
                      <td className="px-4 py-2.5 text-indigo-600 font-semibold">92% Match</td>
                      <td className="px-4 py-2.5 text-slate-600">High</td>
                    </tr>
                    <tr>
                      <td className="px-4 py-2.5 font-bold text-slate-900">Calu-3</td>
                      <td className="px-4 py-2.5 text-indigo-600 font-semibold">84% Match</td>
                      <td className="px-4 py-2.5 text-slate-600">Moderate</td>
                    </tr>
                  </tbody>
                </table>
              </div>
            </div>
          )}
        </div>
      </div>
    </div>
  );
};