/**
 * TypeScript port of chainofcustody/three_prime/filtering_on_target.py
 * — greedy_mirna_cover().
 *
 * Given a set of target cell types (to protect) and off-target cell types
 * (to suppress), finds a minimal set of miRNAs that are:
 *   - silent in ALL target cells (mean RPM < targetThreshold)
 *   - collectively expressed (≥ coverThreshold) in every off-target cell
 *
 * When multiple targets are specified the intersection of their candidate
 * sets is used (a miRNA must be silent in all targets simultaneously).
 */

import type { MirnaData, GreedyParams, GreedyResult, GreedyStepResult } from "./types";

function getMeanRpm(
  meanMatrix: Record<string, Record<string, number>>,
  mirna: string,
  cellType: string
): number {
  return meanMatrix[mirna]?.[cellType] ?? 0;
}

export function greedyCover(data: MirnaData, params: GreedyParams): GreedyResult {
  const { targets, offTargets, targetThreshold, coverThreshold, maxMirnas } = params;

  if (targets.length === 0 || offTargets.length === 0) {
    return {
      success: false,
      selectedMirnas: [],
      steps: [],
      uncovered: [...offTargets],
      allOffTargets: [...offTargets],
    };
  }

  // Candidates: miRNAs silent in ALL target cells
  const candidateMirnas = data.mirnas.filter((mirna) =>
    targets.every(
      (target) => getMeanRpm(data.mean_matrix, mirna, target) < targetThreshold
    )
  );

  // Pre-compute coverage sets: for each candidate, which off-targets does it cover?
  const mirnaCoverage = new Map<string, Set<string>>();
  for (const mirna of candidateMirnas) {
    const covered = new Set<string>();
    for (const ot of offTargets) {
      if (getMeanRpm(data.mean_matrix, mirna, ot) >= coverThreshold) {
        covered.add(ot);
      }
    }
    mirnaCoverage.set(mirna, covered);
  }

  const uncovered = new Set(offTargets);
  const selected: string[] = [];
  const steps: GreedyStepResult[] = [];
  const selectedSet = new Set<string>();

  while (uncovered.size > 0 && selected.length < maxMirnas) {
    let bestMirna: string | null = null;
    let bestNewCover = new Set<string>();

    for (const mirna of candidateMirnas) {
      if (selectedSet.has(mirna)) continue;
      const newCover = new Set<string>();
      for (const ot of mirnaCoverage.get(mirna)!) {
        if (uncovered.has(ot)) newCover.add(ot);
      }
      if (newCover.size > bestNewCover.size) {
        bestMirna = mirna;
        bestNewCover = newCover;
      }
    }

    if (bestMirna === null || bestNewCover.size === 0) break;

    // Compute mean RPM across all target cells
    const targetRpm =
      targets.reduce(
        (sum, t) => sum + getMeanRpm(data.mean_matrix, bestMirna!, t),
        0
      ) / targets.length;

    steps.push({
      mirnaId: bestMirna,
      seed: data.mir_to_seed[bestMirna] ?? "",
      matureSeq: data.mature_seqs[bestMirna] ?? "",
      targetRpm,
      newlyCovered: [...bestNewCover],
    });

    selected.push(bestMirna);
    selectedSet.add(bestMirna);
    for (const ot of bestNewCover) uncovered.delete(ot);
  }

  return {
    success: uncovered.size === 0,
    selectedMirnas: selected,
    steps,
    uncovered: [...uncovered],
    allOffTargets: [...offTargets],
  };
}
