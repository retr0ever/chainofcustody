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

/* ──────────────────────────────────────────────────────────────────────
   Evaluation pipeline types (mirrors Python score_parsed + compute_fitness)
   ────────────────────────────────────────────────────────────────────── */

export type MetricStatus = "GREEN" | "AMBER" | "RED" | "GREY";

export interface SequenceInfo {
  total_length: number;
  full_length: number;
  cap5_length: number;
  utr5_length: number;
  cds_length: number;
  utr3_length: number;
  num_codons: number;
}

export interface Utr5Accessibility {
  mfe: number | null;
  mfe_per_nt: number | null;
  utr5_length: number;
  fold_window: number;
  status: string;
  message: string;
}

export interface GlobalMfe {
  mfe: number;
  mfe_per_nt: number;
  length: number;
  method: string;
}

export interface StructureScores {
  utr5_accessibility: Utr5Accessibility;
  global_mfe: GlobalMfe;
}

export interface ManufacturingScores {
  gc_windows: {
    pass: boolean;
    violations: Array<{ position: number; gc_content: number; issue: string }>;
    windows_checked: number;
  };
  homopolymers: {
    pass: boolean;
    violations: Array<{ position: number; nucleotide: string; length: number; sequence: string }>;
    max_allowed: number;
  };
  restriction_sites: {
    pass: boolean;
    violations: Array<{ position: number; enzyme: string; sequence: string; strand: string }>;
    enzymes_checked: string[];
  };
  uorfs: {
    pass: boolean;
    count: number;
    positions: number[];
    violations: Array<{ position: number; sequence: string }>;
  };
  total_violations: number;
  utr5_violations: number;
  overall_pass: boolean;
}

export interface StabilityScores {
  gc3: number;
  mfe_per_nt: number;
  stability_score: number;
  status: MetricStatus;
}

export interface RibonnScores {
  mean_te: number;
  target_cell_type: string;
  target_te: number;
  mean_off_target_te: number;
  per_tissue: Record<string, number> | null;
  status: MetricStatus;
  message: string;
}

export interface ScoreReport {
  sequence_info: SequenceInfo;
  structure_scores: StructureScores;
  manufacturing_scores: ManufacturingScores;
  stability_scores: StabilityScores;
  ribonn_scores: RibonnScores;
  summary: Record<string, MetricStatus>;
}

export interface FitnessMetricScore {
  value: number;
  weight: number;
  weighted: number;
  status: MetricStatus;
}

export interface FitnessResult {
  scores: Record<string, FitnessMetricScore>;
  overall: number;
  suggestions: Array<{ metric: string; priority: string; action: string }>;
}

export interface CandidateSequence {
  utr5: string;
  cds: string;
  utr3: string;
}

export interface CandidateEvaluation {
  id: string;
  name: string;
  gene: string;
  target_cell_type: string;
  created_at: string;
  sequence: CandidateSequence;
  report: ScoreReport;
  fitness: FitnessResult;
}
