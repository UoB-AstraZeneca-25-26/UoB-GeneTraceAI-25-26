import React from 'react';
import { ExternalLink } from 'lucide-react';

export const AboutPage: React.FC = () => {
  return (
    <div className="max-w-3xl mx-auto space-y-6">
      <div>
        <h2 className="text-2xl font-bold text-slate-900">About Cell Line Finder</h2>
        <p className="text-slate-600 mt-1">
          A decision-support tool for choosing human cancer cell-line models using integrated multi-omics evidence.
        </p>
      </div>

      {/* Purpose */}
      <div className="bg-white rounded-xl shadow-sm border border-slate-200 p-6 space-y-3">
        <h3 className="text-base font-bold text-slate-900">Purpose</h3>
        <p className="text-sm text-slate-600 leading-relaxed">
          The cell line you pick shapes how relevant, robust and reproducible an experiment is. Yet the evidence
          needed to choose well is scattered across many public datasets, each with its own format and conventions.
          Cell Line Finder harmonises that evidence into one queryable system: you give it a target gene (or a
          combination), and it returns a ranked shortlist of cell lines together with the molecular evidence behind
          each recommendation, so the choice is transparent rather than a guess.
        </p>
        <p className="text-xs text-slate-500 leading-relaxed">
          Rankings are model predictions derived from public data — a starting point for selection, to be confirmed
          experimentally, not a substitute for it.
        </p>
      </div>

      {/* Layers of evidence */}
      <div className="bg-white rounded-xl shadow-sm border border-slate-200 p-6 space-y-3">
        <h3 className="text-base font-bold text-slate-900">How it combines evidence</h3>
        <p className="text-sm text-slate-600 leading-relaxed">
          Transcription and translation are regulated independently, so no single measurement decides whether a line
          is a good model — each molecular layer answers a different question, and they are strongest together.
          The tool integrates complementary layers:
        </p>
        <div className="divide-y divide-slate-100 border border-slate-100 rounded-lg">
          {LAYERS.map((l) => (
            <div key={l.layer} className="flex gap-4 px-4 py-2.5">
              <span className="w-32 shrink-0 text-sm font-semibold text-slate-800">{l.layer}</span>
              <span className="text-sm text-slate-600">{l.answers}</span>
            </div>
          ))}
        </div>
        <p className="text-xs text-slate-500 leading-relaxed">
          Transcriptomic, proteomic and genomic evidence are scored; copy number acts as a directional modifier on
          confidence. Expression and proteomics are reported as relative levels (Low / Moderate / High); mutations,
          fusions and copy number as present or absent.
        </p>
      </div>

      {/* Data sources */}
      <div className="bg-white rounded-xl shadow-sm border border-slate-200 p-6 space-y-4">
        <h3 className="text-base font-bold text-slate-900">Data sources</h3>
        <p className="text-sm text-slate-600 leading-relaxed">
          The corpus draws on established public programmes, each contributing a different piece of the picture:
        </p>
        <div className="grid grid-cols-1 sm:grid-cols-2 gap-3">
          {SOURCES.map((s) => (
            <div key={s.name} className="border border-slate-200 rounded-lg p-4">
              <div className="flex items-start justify-between gap-2">
                <h4 className="text-sm font-bold text-slate-800">{s.name}</h4>
                <a
                  href={s.url}
                  target="_blank"
                  rel="noopener noreferrer"
                  className="text-slate-400 hover:text-mulberry-600 shrink-0"
                  aria-label={`Open ${s.name}`}
                >
                  <ExternalLink className="w-3.5 h-3.5" />
                </a>
              </div>
              <p className="text-xs text-slate-600 mt-1.5 leading-relaxed">{s.body}</p>
              <p className="text-[11px] text-slate-400 mt-2">{s.cite}</p>
            </div>
          ))}
        </div>
        <p className="text-[11px] text-slate-400 leading-relaxed pt-1">
          Cell-line identities are resolved through Cellosaurus to guard against misidentified or contaminated lines.
          COSMIC-derived driver annotations are used under its academic licence and are not redistributed.
        </p>
      </div>
    </div>
  );
};

const LAYERS: { layer: string; answers: string }[] = [
  { layer: 'Transcriptomics', answers: 'Is the gene transcribed, and how strongly? (broad coverage, but a weak predictor of protein)' },
  { layer: 'Proteomics', answers: 'Is the protein actually present, and roughly how much? (biased toward abundant proteins)' },
  { layer: 'Somatic mutation', answers: 'Is the gene altered, and is the alteration a known driver?' },
  { layer: 'Gene fusion', answers: 'Does a fusion involving the gene exist in this line?' },
  { layer: 'Copy number', answers: 'Is the gene amplified or deleted? (used to modify confidence)' },
];

const SOURCES: { name: string; url: string; body: string; cite: string }[] = [
  {
    name: 'DepMap / CCLE',
    url: 'https://depmap.org/portal/',
    body: 'The Cancer Cell Line Encyclopedia, delivered through the Cancer Dependency Map (Broad Institute) on a quarterly cadence. Provides RNA-seq expression (log2 TPM+1), somatic mutations, gene fusions, copy number and derived signatures across ~1,100 cancer cell lines.',
    cite: 'Ghandi et al., Nature 2019 · Tsherniak et al., Cell 2017',
  },
  {
    name: 'Human Protein Atlas (HPA)',
    url: 'https://www.proteinatlas.org/',
    body: 'RNA and protein expression across human tissues and cell lines, from antibody-based imaging and mass spectrometry. Contributes cell-line transcriptomics (nTPM) and independently benchmarks lines against matched tumour cohorts.',
    cite: 'Uhlén et al., Science 2015',
  },
  {
    name: 'Cell-line proteomics (CCLE & ProCan)',
    url: 'https://depmap.sanger.ac.uk/documentation/datasets/proteomics/',
    body: 'Two complementary mass-spectrometry maps of the cell-line proteome — CCLE (TMT, ~375 lines; Nusinow) and ProCan-DepMapSanger (DIA, 949 lines; Gonçalves). Detection is biased toward abundant proteins, so absence signals weak evidence rather than certain absence.',
    cite: 'Nusinow et al., Cell 2020 · Gonçalves et al., Cancer Cell 2022',
  },
  {
    name: 'Gene Expression Omnibus (GEO)',
    url: 'https://www.ncbi.nlm.nih.gov/geo/',
    body: 'NCBI’s public repository of functional-genomics data submitted by independent labs. Adds an expression compendium; because submissions are heterogeneous, values are standardised within method rather than assumed comparable across sources.',
    cite: 'Edgar et al., Nucleic Acids Res 2002',
  },
  {
    name: 'COSMIC — Cancer Gene Census',
    url: 'https://cancer.sanger.ac.uk/cosmic',
    body: 'A conservative, expert-curated catalogue of genes causally implicated in cancer (with driver / tumour-suppressor tiers), maintained at the Wellcome Sanger Institute. Supplies the driver annotations behind the “driver” tag.',
    cite: 'Sondka et al., Nature Reviews Cancer 2018',
  },
  {
    name: 'Cellosaurus',
    url: 'https://www.cellosaurus.org/',
    body: 'A knowledge resource cataloguing >100,000 cell lines with standardised identifiers, synonyms and provenance. Used to resolve cell-line identity (ACH ids) and to flag potentially misidentified or contaminated lines.',
    cite: 'Bairoch, J Biomol Tech 2018',
  },
  {
    name: 'Ensembl',
    url: 'https://www.ensembl.org/',
    body: 'Genome-annotation database providing gene models and stable identifiers (ENSG ids). Anchors genes to a consistent reference so evidence from different sources maps to the same gene.',
    cite: 'Cunningham et al., Nucleic Acids Res 2022',
  },
  {
    name: 'HGNC',
    url: 'https://www.genenames.org/',
    body: 'The HUGO Gene Nomenclature Committee — the authority for approved human gene symbols, full names and aliases. Powers gene-symbol search, validation and the “did you mean?” suggestions.',
    cite: 'Seal et al., Nucleic Acids Res 2023',
  },
];
