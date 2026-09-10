import { CellLineData } from '../types';

export const ALL_GENES = [
  "EGFR",
  "KRAS",
  "BRAF",
  "TP53",
  "PIK3CA",
  "ALK",
  "MET",
  "ERBB2",
  "PTEN",
  "MYC",
  "BCR::ABL1"
];

/** Demo gene reference (symbol ↔ full name ↔ Ensembl gene id). Powers the
 *  searchable dropdown and the Reference lookup. To be replaced by the full gene
 *  table later. */
export const GENE_REF: { symbol: string; name: string; ensg: string }[] = [
  { symbol: 'EGFR', name: 'Epidermal growth factor receptor', ensg: 'ENSG00000146648' },
  { symbol: 'KRAS', name: 'KRAS proto-oncogene, GTPase', ensg: 'ENSG00000133703' },
  { symbol: 'NRAS', name: 'NRAS proto-oncogene, GTPase', ensg: 'ENSG00000213281' },
  { symbol: 'HRAS', name: 'HRas proto-oncogene, GTPase', ensg: 'ENSG00000174775' },
  { symbol: 'BRAF', name: 'B-Raf proto-oncogene, serine/threonine kinase', ensg: 'ENSG00000157764' },
  { symbol: 'TP53', name: 'Tumor protein p53', ensg: 'ENSG00000141510' },
  { symbol: 'PIK3CA', name: 'Phosphatidylinositol-4,5-bisphosphate 3-kinase catalytic subunit alpha', ensg: 'ENSG00000121879' },
  { symbol: 'PTEN', name: 'Phosphatase and tensin homolog', ensg: 'ENSG00000171862' },
  { symbol: 'ALK', name: 'ALK receptor tyrosine kinase', ensg: 'ENSG00000171094' },
  { symbol: 'MET', name: 'MET proto-oncogene, receptor tyrosine kinase', ensg: 'ENSG00000105976' },
  { symbol: 'ERBB2', name: 'Erb-b2 receptor tyrosine kinase 2 (HER2)', ensg: 'ENSG00000141736' },
  { symbol: 'ERBB3', name: 'Erb-b2 receptor tyrosine kinase 3', ensg: 'ENSG00000065361' },
  { symbol: 'MYC', name: 'MYC proto-oncogene, bHLH transcription factor', ensg: 'ENSG00000136997' },
  { symbol: 'MYCN', name: 'MYCN proto-oncogene, bHLH transcription factor', ensg: 'ENSG00000134323' },
  { symbol: 'RB1', name: 'RB transcriptional corepressor 1', ensg: 'ENSG00000139687' },
  { symbol: 'CDK4', name: 'Cyclin dependent kinase 4', ensg: 'ENSG00000135446' },
  { symbol: 'CDK6', name: 'Cyclin dependent kinase 6', ensg: 'ENSG00000105810' },
  { symbol: 'CDKN2A', name: 'Cyclin dependent kinase inhibitor 2A', ensg: 'ENSG00000147889' },
  { symbol: 'CCND1', name: 'Cyclin D1', ensg: 'ENSG00000110092' },
  { symbol: 'CCNE1', name: 'Cyclin E1', ensg: 'ENSG00000105173' },
  { symbol: 'MDM2', name: 'MDM2 proto-oncogene', ensg: 'ENSG00000135679' },
  { symbol: 'MDM4', name: 'MDM4 regulator of p53', ensg: 'ENSG00000198625' },
  { symbol: 'AKT1', name: 'AKT serine/threonine kinase 1', ensg: 'ENSG00000142208' },
  { symbol: 'MTOR', name: 'Mechanistic target of rapamycin kinase', ensg: 'ENSG00000198793' },
  { symbol: 'STK11', name: 'Serine/threonine kinase 11 (LKB1)', ensg: 'ENSG00000118046' },
  { symbol: 'KEAP1', name: 'Kelch like ECH associated protein 1', ensg: 'ENSG00000079999' },
  { symbol: 'NFE2L2', name: 'Nuclear factor erythroid 2 like 2 (NRF2)', ensg: 'ENSG00000116044' },
  { symbol: 'SMARCA4', name: 'SWI/SNF related, matrix associated, actin dependent regulator of chromatin, subfamily a, member 4', ensg: 'ENSG00000127616' },
  { symbol: 'ARID1A', name: 'AT-rich interaction domain 1A', ensg: 'ENSG00000117713' },
  { symbol: 'VHL', name: 'Von Hippel-Lindau tumor suppressor', ensg: 'ENSG00000134086' },
  { symbol: 'NF1', name: 'Neurofibromin 1', ensg: 'ENSG00000196712' },
  { symbol: 'NF2', name: 'Neurofibromin 2 (merlin)', ensg: 'ENSG00000186575' },
  { symbol: 'APC', name: 'APC regulator of WNT signaling pathway', ensg: 'ENSG00000134982' },
  { symbol: 'CTNNB1', name: 'Catenin beta 1', ensg: 'ENSG00000168036' },
  { symbol: 'SMAD4', name: 'SMAD family member 4', ensg: 'ENSG00000141646' },
  { symbol: 'GNAS', name: 'GNAS complex locus', ensg: 'ENSG00000087460' },
  { symbol: 'IDH1', name: 'Isocitrate dehydrogenase (NADP(+)) 1', ensg: 'ENSG00000138413' },
  { symbol: 'IDH2', name: 'Isocitrate dehydrogenase (NADP(+)) 2', ensg: 'ENSG00000182054' },
  { symbol: 'FLT3', name: 'Fms related receptor tyrosine kinase 3', ensg: 'ENSG00000122025' },
  { symbol: 'KIT', name: 'KIT proto-oncogene, receptor tyrosine kinase', ensg: 'ENSG00000157404' },
  { symbol: 'PDGFRA', name: 'Platelet derived growth factor receptor alpha', ensg: 'ENSG00000134853' },
  { symbol: 'FGFR1', name: 'Fibroblast growth factor receptor 1', ensg: 'ENSG00000077782' },
  { symbol: 'FGFR2', name: 'Fibroblast growth factor receptor 2', ensg: 'ENSG00000066468' },
  { symbol: 'FGFR3', name: 'Fibroblast growth factor receptor 3', ensg: 'ENSG00000068078' },
  { symbol: 'ROS1', name: 'ROS proto-oncogene 1, receptor tyrosine kinase', ensg: 'ENSG00000047936' },
  { symbol: 'RET', name: 'Ret proto-oncogene', ensg: 'ENSG00000165731' },
  { symbol: 'NTRK1', name: 'Neurotrophic receptor tyrosine kinase 1', ensg: 'ENSG00000198400' },
  { symbol: 'ABL1', name: 'ABL proto-oncogene 1, non-receptor tyrosine kinase', ensg: 'ENSG00000097007' },
  { symbol: 'JAK2', name: 'Janus kinase 2', ensg: 'ENSG00000096968' },
  { symbol: 'BRCA1', name: 'BRCA1 DNA repair associated', ensg: 'ENSG00000012048' },
  { symbol: 'BRCA2', name: 'BRCA2 DNA repair associated', ensg: 'ENSG00000139618' },
  { symbol: 'ATM', name: 'ATM serine/threonine kinase', ensg: 'ENSG00000149311' },
  { symbol: 'CDH1', name: 'Cadherin 1', ensg: 'ENSG00000039068' },
  { symbol: 'VEGFA', name: 'Vascular endothelial growth factor A', ensg: 'ENSG00000112715' },
  { symbol: 'ESR1', name: 'Estrogen receptor 1', ensg: 'ENSG00000091831' },
  { symbol: 'AR', name: 'Androgen receptor', ensg: 'ENSG00000169083' },
  { symbol: 'BCL2', name: 'BCL2 apoptosis regulator', ensg: 'ENSG00000171791' },
  { symbol: 'MCL1', name: 'MCL1 apoptosis regulator, BCL2 family member', ensg: 'ENSG00000143384' },
  { symbol: 'CCND3', name: 'Cyclin D3', ensg: 'ENSG00000112576' },
];

/** Full gene universe for the searchable dropdown (symbols only). */
export const GENE_UNIVERSE: string[] = GENE_REF.map((g) => g.symbol);

export const ALL_LINEAGES = [
  "Any",
  "Lung",
  "Breast",
  "Colon",
  "Liver",
  "Blood",
  "Cervix",
  "Prostate",
  "Brain"
];

export function generateMockData(): CellLineData[] {
  const cellLines = [
    { name: "HCC827", tissue: "Lung" },
    { name: "A549", tissue: "Lung" },
    { name: "NCI-H1975", tissue: "Lung" },
    { name: "Calu-3", tissue: "Lung" },
    { name: "MCF7", tissue: "Breast" },
    { name: "HeLa", tissue: "Cervix" },
    { name: "PC3", tissue: "Prostate" },
    { name: "K562", tissue: "Blood" },
    { name: "Jurkat", tissue: "Blood" },
    { name: "HepG2", tissue: "Liver" },
    { name: "U87", tissue: "Brain" },
    { name: "HT29", tissue: "Colon" }
  ];

  const data: CellLineData[] = [];

  cellLines.forEach(({ name, tissue }) => {
    ALL_GENES.forEach((gene) => {
      let tpm = Math.floor(Math.random() * 140) + 10;
      
      // Known biological patterns for realism
      if (gene === "BCR::ABL1") {
        tpm = name === "K562" ? 220 : 0;
      } else if (gene === "EGFR" && (name === "HCC827" || name === "A549")) {
        tpm = Math.floor(Math.random() * 40) + 160;
      }

      const protein = Math.round(tpm * (0.6 + Math.random() * 0.5));
      const hasMutation = Math.random() < 0.25;

      data.push({
        cellLine: name,
        tissue,
        gene,
        tpmExpression: tpm,
        proteinExpression: protein,
        hasSomaticVariant: hasMutation
      });
    });
  });

  return data;
}

export const MOCK_DATA = generateMockData();