import { CellLineData, QueryParams, ScoredResult } from '../types';

export function runScoringPipeline(
  dataset: CellLineData[],
  params: QueryParams
): ScoredResult[] {
  const { targets, exclusions, lineage, topK } = params;

  // 1. Filter by Lineage
  const lineageFiltered = lineage !== "Any" 
    ? dataset.filter(d => d.tissue.toLowerCase() === lineage.toLowerCase())
    : dataset;

  // 2. Identify unique cell lines in this subset
  const uniqueCellLines = Array.from(new Set(lineageFiltered.map(d => d.cellLine)));

  const scored: ScoredResult[] = [];

  uniqueCellLines.forEach(cl => {
    const clData = lineageFiltered.filter(d => d.cellLine === cl);
    const tissue = clData[0]?.tissue || "Unknown";

    // Target metrics
    const targetRecords = clData.filter(d => targets.includes(d.gene));
    const meanTpm = targetRecords.length > 0
      ? targetRecords.reduce((sum, r) => sum + r.tpmExpression, 0) / targetRecords.length
      : 0;

    const mutations = targetRecords.filter(r => r.hasSomaticVariant).length;

    // Exclusion penalty calculation
    const exclusionRecords = clData.filter(d => exclusions.includes(d.gene));
    const meanExclusionTpm = exclusionRecords.length > 0
      ? exclusionRecords.reduce((sum, r) => sum + r.tpmExpression, 0) / exclusionRecords.length
      : 0;

    // Score calculation
    const rawScore = (meanTpm / 1.5) - (meanExclusionTpm / 2);
    const confidenceScore = Math.min(100, Math.max(0, Math.round(rawScore * 10) / 10));

    scored.push({
      cellLine: cl,
      tissue,
      meanTpm: Math.round(meanTpm * 10) / 10,
      mutations,
      confidenceScore,
      exclusionPenalty: Math.round(meanExclusionTpm * 10) / 10
    });
  });

  // 3. Sort by score descending and return Top K
  return scored
    .sort((a, b) => b.confidenceScore - a.confidenceScore)
    .slice(0, topK);
}