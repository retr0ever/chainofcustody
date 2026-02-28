export interface MirnaData {
  cell_types: string[];
  mirnas: string[];
  /** mirna_id → { cell_type → mean_rpm } (sparse: zero entries omitted) */
  mean_matrix: Record<string, Record<string, number>>;
  mir_to_seed: Record<string, string>;
  mature_seqs: Record<string, string>;
}

export interface GreedyParams {
  targets: string[];
  offTargets: string[];
  targetThreshold: number;
  coverThreshold: number;
  maxMirnas: number;
}

export interface GreedyStepResult {
  mirnaId: string;
  seed: string;
  matureSeq: string;
  /** Mean RPM across all target cells (average) */
  targetRpm: number;
  /** Off-target cells newly covered by this miRNA */
  newlyCovered: string[];
}

export interface GreedyResult {
  success: boolean;
  selectedMirnas: string[];
  steps: GreedyStepResult[];
  uncovered: string[];
  allOffTargets: string[];
}
