import React from 'react';
import {
  Target,
  Layers,
  Ban,
  Table2,
  LayoutGrid,
  ChevronRight,
  Microscope,
  ArrowRight,
  Dna,
  ListOrdered,
} from 'lucide-react';

interface GuidePageProps {
  onStart: () => void;
}

/**
 * First-run walkthrough, structured around the user's mental model:
 * what it does → the three ways to search → how to read what comes back.
 */
export const GuidePage: React.FC<GuidePageProps> = ({ onStart }) => {
  return (
    <div className="max-w-3xl mx-auto space-y-8">
      <div>
        <h2 className="text-2xl font-bold text-slate-900 font-display tracking-tight">How to use Cell Line Finder</h2>
        <p className="text-slate-600 mt-1">A two-minute orientation — what it does, how to search, and how to read the results.</p>
      </div>

      {/* The idea in one line */}
      <section className="bg-plum text-white rounded-xl p-5">
        <div className="flex items-center justify-center gap-3 sm:gap-4 flex-wrap text-center">
          <span className="inline-flex items-center gap-2 font-semibold"><Dna className="w-5 h-5 text-mulberry-300" /> A gene you care about</span>
          <ArrowRight className="w-5 h-5 text-mulberry-300" />
          <span className="inline-flex items-center gap-2 font-semibold"><ListOrdered className="w-5 h-5 text-mulberry-300" /> A ranked list of cell lines</span>
          <ArrowRight className="w-5 h-5 text-mulberry-300" />
          <span className="inline-flex items-center gap-2 font-semibold"><Microscope className="w-5 h-5 text-mulberry-300" /> The evidence behind each one</span>
        </div>
        <p className="text-xs text-slate-300 text-center mt-3 leading-relaxed">
          You give it a target gene; it ranks human cancer cell lines by how well they model that gene, using
          integrated multi-omics evidence — then shows you why each line ranked where it did.
        </p>
      </section>

      {/* Three ways to search */}
      <section className="space-y-3">
        <h3 className="text-lg font-bold text-slate-900">Three ways to search</h3>
        <div className="grid grid-cols-1 sm:grid-cols-3 gap-3">
          <ModeCard
            icon={Target}
            title="Single gene"
            example="e.g. EGFR"
            body="One target gene. Ranks every line by how strongly the evidence supports that gene."
          />
          <ModeCard
            icon={Layers}
            title="Multiple targets"
            example="e.g. EGFR + KRAS"
            body="Two or more targets. Runs a joint ranking that combines the evidence across all of them."
          />
          <ModeCard
            icon={Ban}
            title="Selectivity"
            example="e.g. BRAF, not EGFR"
            body="A target plus an exclusion. Finds lines high in the target and low in the gene you exclude."
          />
        </div>
        <p className="text-xs text-slate-500">
          On the Find Cell Lines page, add target genes, optionally add an exclusion, set how many results with the
          slider, and hit <span className="font-semibold text-slate-700">Find Cell Lines</span>. The mode is chosen
          for you from what you enter.
        </p>
      </section>

      {/* The two scores */}
      <section className="bg-white rounded-xl shadow-sm border border-slate-200 p-6 space-y-3">
        <h3 className="text-base font-bold text-slate-900">Two scores to read</h3>
        <div className="grid grid-cols-1 sm:grid-cols-2 gap-4">
          <div>
            <h4 className="text-sm font-bold text-slate-800">Global score</h4>
            <p className="text-xs text-slate-600 mt-1 leading-relaxed">
              How the line ranks across <em>all</em> cell lines (0–1, higher is better). The top scores cluster
              together, so treat the leading lines as a shortlist rather than one winner.
            </p>
          </div>
          <div>
            <h4 className="text-sm font-bold text-slate-800">In-lineage score</h4>
            <p className="text-xs text-slate-600 mt-1 leading-relaxed">
              The same idea, but only among lines of the <em>same tissue lineage</em> — useful when you need a model
              from a particular tissue. The Lineage filter re-ranks the table by this score.
            </p>
          </div>
        </div>
      </section>

      {/* Reading the results — views */}
      <section className="space-y-3">
        <h3 className="text-lg font-bold text-slate-900">Three ways to view results</h3>
        <div className="grid grid-cols-1 sm:grid-cols-3 gap-3">
          <ViewCard icon={Table2} title="Table" body="The ranked list with the global and in-lineage scores. Filter by lineage, or Inspect any line." />
          <ViewCard icon={LayoutGrid} title="Evidence grid" body="Cell lines against their evidence layers — levels and alterations side by side. For multi and selectivity queries, layers are shown per gene." />
          <ViewCard icon={Layers} title="Lineages" body="Results grouped by tissue lineage, so you can compare candidates within and across tissues." />
        </div>
        <p className="text-xs text-slate-500">
          The panel on the right of the results always explains the view you’re looking at and what its numbers mean.
        </p>
      </section>

      {/* Reading the evidence */}
      <section className="bg-plum text-white rounded-xl p-6 space-y-4">
        <div className="flex items-center gap-2">
          <Microscope className="w-5 h-5 text-mulberry-300" />
          <h3 className="text-base font-bold">Reading the evidence</h3>
        </div>
        <p className="text-sm text-slate-300 leading-relaxed">
          Evidence comes in two kinds, read differently — in the grid and in a line’s profile:
        </p>
        <div className="grid grid-cols-1 sm:grid-cols-2 gap-3">
          <div className="bg-plum-accent rounded-lg p-4">
            <h4 className="text-sm font-bold text-white">Levels — Low / Moderate / High</h4>
            <p className="text-xs text-slate-300 mt-1 leading-relaxed">
              Expression and proteomics show a <span className="text-white font-semibold">relative level</span> — where
              the line sits for that gene compared with other cell lines, not an absolute amount. Chips show which
              datasets measured it (expression: DepMap / HPA / GEO · proteomics: ProCan / CCLE); more sources = more
              independent evidence.
            </p>
          </div>
          <div className="bg-plum-accent rounded-lg p-4">
            <h4 className="text-sm font-bold text-white">Alterations — present / absent</h4>
            <p className="text-xs text-slate-300 mt-1 leading-relaxed">
              Mutation, fusion and copy number are categorical: each is simply present or absent for that line. A{' '}
              <span className="font-mono text-mulberry-300">driver</span> tag means a known cancer-driver alteration,
              not just any variant.
            </p>
          </div>
        </div>
      </section>

      <div className="flex justify-center pt-2">
        <button
          onClick={onStart}
          className="inline-flex items-center gap-2 bg-mulberry-600 hover:bg-mulberry-700 text-white font-medium px-6 py-3 rounded-lg shadow-sm transition"
        >
          <span>Start your first search</span>
          <ChevronRight className="w-4 h-4" />
        </button>
      </div>
    </div>
  );
};

const ModeCard: React.FC<{ icon: React.ElementType; title: string; example: string; body: string }> = ({
  icon: Icon,
  title,
  example,
  body,
}) => (
  <div className="bg-white rounded-xl shadow-sm border border-slate-200 p-4">
    <div className="flex items-center gap-2">
      <span className="w-8 h-8 rounded-lg bg-mulberry-50 text-mulberry-600 flex items-center justify-center">
        <Icon className="w-4 h-4" />
      </span>
      <h4 className="text-sm font-bold text-slate-800">{title}</h4>
    </div>
    <p className="text-[11px] font-mono text-mulberry-700 mt-2">{example}</p>
    <p className="text-xs text-slate-600 mt-1 leading-relaxed">{body}</p>
  </div>
);

const ViewCard: React.FC<{ icon: React.ElementType; title: string; body: string }> = ({ icon: Icon, title, body }) => (
  <div className="bg-white rounded-xl shadow-sm border border-slate-200 p-4">
    <div className="flex items-center gap-2">
      <Icon className="w-4 h-4 text-mulberry-600" />
      <h4 className="text-sm font-bold text-slate-800">{title}</h4>
    </div>
    <p className="text-xs text-slate-600 mt-1.5 leading-relaxed">{body}</p>
  </div>
);
