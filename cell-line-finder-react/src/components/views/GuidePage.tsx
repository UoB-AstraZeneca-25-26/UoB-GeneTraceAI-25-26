import React from 'react';
import {
  Search,
  Target,
  Ban,
  SlidersHorizontal,
  Table2,
  Network,
  GitFork,
  Layers,
  ChevronRight,
  Microscope,
} from 'lucide-react';

interface GuidePageProps {
  onStart: () => void;
}

/**
 * A first-run walkthrough. It answers three questions a new user has: what this
 * tool does, how to run a search, and how to read what comes back — including
 * the two number types (percentiles vs. categorical calls) that trip people up.
 */
export const GuidePage: React.FC<GuidePageProps> = ({ onStart }) => {
  return (
    <div className="max-w-3xl mx-auto space-y-8">
      <div>
        <h2 className="text-2xl font-bold text-slate-900 font-display tracking-tight">How to use Cell Line Finder</h2>
        <p className="text-slate-600 mt-1">
          A two-minute orientation. Cell Line Finder ranks human cancer cell lines by how well they match a gene
          you care about, using integrated multi-omics evidence — then lets you inspect why each line ranked where it did.
        </p>
      </div>

      {/* Step-by-step */}
      <section className="bg-white rounded-xl shadow-sm border border-slate-200 divide-y divide-slate-100">
        <Step
          n={1}
          icon={Target}
          title="Pick your target gene(s)"
          body="On the Find Cell Lines page, add one or more target genes — type any symbol or tap a quick-pick. One target ranks lines for that gene; two or more runs a joint ranking across all of them."
        />
        <Step
          n={2}
          icon={Ban}
          title="Optionally, exclude a gene"
          body="Add an exclusion gene to rank for selectivity — lines where your target is high and the excluded gene is low. Useful when you want a background that lacks a confounder."
        />
        <Step
          n={3}
          icon={SlidersHorizontal}
          title="Choose how many results"
          body="Use the slider to set how many lines to return (3–30). You can change this any time and re-run. Clear filters resets everything to start over."
        />
        <Step
          n={4}
          icon={Search}
          title="Run the search"
          body="Hit Find Cell Lines. You'll land on Results, with the ranked table and a highlighted top pick."
          last
        />
      </section>

      {/* Reading the results */}
      <section className="space-y-3">
        <h3 className="text-lg font-bold text-slate-900">Reading your results</h3>
        <div className="grid grid-cols-1 sm:grid-cols-2 gap-3">
          <ViewCard
            icon={Table2}
            title="Table"
            body="The ranked list. The number is the score; the marker beside it shows a line's position within the shown set — lines in your top N all score highly, so treat them as comparable candidates."
          />
          <ViewCard
            icon={Network}
            title="Network"
            body="A radial map: distance from the centre follows rank (nearest = highest-ranked), and colour marks tissue lineage."
          />
          <ViewCard
            icon={GitFork}
            title="Omics map"
            body="Traces the evidence flow from the gene, through each measured layer, to the top lines. Ribbon thickness follows the per-layer signal."
          />
          <ViewCard
            icon={Layers}
            title="Lineages"
            body="Groups the results by tissue lineage so you can compare candidates within and across tissues."
          />
        </div>
      </section>

      {/* How to read the evidence — the part people find confusing */}
      <section className="bg-slate-900 text-slate-100 rounded-xl p-6 space-y-4">
        <div className="flex items-center gap-2">
          <Microscope className="w-5 h-5 text-indigo-300" />
          <h3 className="text-base font-bold">Reading the evidence in a profile</h3>
        </div>
        <p className="text-sm text-slate-300 leading-relaxed">
          Click <span className="font-semibold text-white">Inspect</span> on any line to open its profile. The
          evidence comes in two kinds — they're read differently:
        </p>
        <div className="grid grid-cols-1 sm:grid-cols-2 gap-3">
          <div className="bg-slate-800/60 rounded-lg p-4">
            <h4 className="text-sm font-bold text-white">Alteration layers — yes / no</h4>
            <p className="text-xs text-slate-300 mt-1 leading-relaxed">
              Mutation, fusion and copy number are categorical: each is simply present or absent for that line.
              A <span className="font-mono text-indigo-300">driver</span> tag means the mutation is a known
              cancer driver, not just any variant.
            </p>
          </div>
          <div className="bg-slate-800/60 rounded-lg p-4">
            <h4 className="text-sm font-bold text-white">Expression &amp; proteomics — levels</h4>
            <p className="text-xs text-slate-300 mt-1 leading-relaxed">
              Shown as a within-gene percentile: where this line sits among all cell lines for the gene.
              <span className="text-white font-semibold"> 0 = lowest, 100 = highest.</span> It's a relative
              level, not an absolute amount or a probability. Chips show which datasets each level came from
              (expression: DepMap / HPA / GEO · proteomics: ProCan / CCLE).
            </p>
          </div>
        </div>
      </section>

      <div className="flex justify-center pt-2">
        <button
          onClick={onStart}
          className="inline-flex items-center gap-2 bg-indigo-600 hover:bg-indigo-700 text-white font-medium px-6 py-3 rounded-lg shadow-sm transition"
        >
          <span>Start your first search</span>
          <ChevronRight className="w-4 h-4" />
        </button>
      </div>
    </div>
  );
};

const Step: React.FC<{
  n: number;
  icon: React.ElementType;
  title: string;
  body: string;
  last?: boolean;
}> = ({ n, icon: Icon, title, body }) => (
  <div className="flex gap-4 p-5">
    <div className="shrink-0">
      <div className="w-9 h-9 rounded-lg bg-indigo-50 text-indigo-600 flex items-center justify-center font-bold font-mono">
        {n}
      </div>
    </div>
    <div className="min-w-0">
      <div className="flex items-center gap-2">
        <Icon className="w-4 h-4 text-slate-400" />
        <h4 className="text-sm font-bold text-slate-800">{title}</h4>
      </div>
      <p className="text-sm text-slate-600 mt-1 leading-relaxed">{body}</p>
    </div>
  </div>
);

const ViewCard: React.FC<{ icon: React.ElementType; title: string; body: string }> = ({
  icon: Icon,
  title,
  body,
}) => (
  <div className="bg-white rounded-xl shadow-sm border border-slate-200 p-4">
    <div className="flex items-center gap-2">
      <Icon className="w-4 h-4 text-indigo-600" />
      <h4 className="text-sm font-bold text-slate-800">{title}</h4>
    </div>
    <p className="text-xs text-slate-600 mt-1.5 leading-relaxed">{body}</p>
  </div>
);
