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