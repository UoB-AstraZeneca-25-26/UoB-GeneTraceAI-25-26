import React from 'react';

export const AboutPage: React.FC = () => {
  return (
    <div className="max-w-3xl mx-auto space-y-6">
      <div>
        <h2 className="text-2xl font-bold text-slate-900">About Cell Line Finder</h2>
        <p className="text-slate-600 mt-1">
          A decision-support platform for selecting human cell lines using integrated multi-omics evidence.
        </p>
      </div>

      <div className="bg-white rounded-xl shadow-sm border border-slate-200 p-6 space-y-4">
        <h3 className="text-base font-bold text-slate-900">Purpose</h3>
        <p className="text-sm text-slate-600 leading-relaxed">
          Cell line choice influences assay relevance, experimental robustness, and reproducibility in biopharmaceutical research. The platform unifies molecular evidence from public datasets (like DepMap, CCLE, and HPA) into one harmonized, queryable system to provide ranked recommendations with confidence metrics.
        </p>

        <h3 className="text-base font-bold text-slate-900 pt-3">Planned Data Sources</h3>
        <ul className="grid grid-cols-2 gap-2 text-sm text-slate-600">
          <li className="flex items-center space-x-1.5"><span className="text-emerald-500">✓</span><span>DepMap TPM & Proteomics</span></li>
          <li className="flex items-center space-x-1.5"><span className="text-emerald-500">✓</span><span>Human Protein Atlas (HPA)</span></li>
          <li className="flex items-center space-x-1.5"><span className="text-emerald-500">✓</span><span>CCLE Metabolomics</span></li>
          <li className="flex items-center space-x-1.5"><span className="text-emerald-500">✓</span><span>Omics Somatic Mutations</span></li>
          <li className="flex items-center space-x-1.5"><span className="text-emerald-500">✓</span><span>Gene Fusion Filters</span></li>
          <li className="flex items-center space-x-1.5"><span className="text-emerald-500">✓</span><span>miRNA Expression Data</span></li>
        </ul>

        <h3 className="text-base font-bold text-slate-900 pt-3">End-to-End Workflow</h3>
        <div className="bg-slate-900 text-mulberry-300 font-mono text-xs p-4 rounded-lg leading-relaxed">
          User Input &rarr; Evidence Retrieval &rarr; Aggregation & Scoring &rarr; Ranked Recommendations &rarr; Profile Deep Dive & Export
        </div>
      </div>
    </div>
  );
};