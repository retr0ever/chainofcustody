"use client";

import { useState } from "react";
import type { MetricStatus, FitnessMetricScore, ScoreReport } from "@/lib/types";
import { formatMetricName, formatScore, formatPercentage } from "@/lib/format";
import StatusBadge from "./StatusBadge";
import MetricBar from "./MetricBar";

interface MetricDetailProps {
  metricKey: string;
  score: FitnessMetricScore;
  report: ScoreReport;
  suggestions?: Array<{ priority: string; action: string }>;
}

function SubRow({ label, value, mono }: { label: string; value: string; mono?: boolean }) {
  return (
    <div className="flex items-center justify-between py-1.5">
      <span className="text-xs" style={{ color: "var(--text-secondary)" }}>
        {label}
      </span>
      <span
        className={`text-xs ${mono ? "font-mono tabular-nums" : ""}`}
        style={{ color: "var(--text-primary)" }}
      >
        {value}
      </span>
    </div>
  );
}

function Utr5Details({ report }: { report: ScoreReport }) {
  const a = report.structure_scores.utr5_accessibility;
  const g = report.structure_scores.global_mfe;
  return (
    <>
      <SubRow label="MFE" value={a.mfe != null ? `${a.mfe.toFixed(1)} kcal/mol` : "\u2014"} mono />
      <SubRow label="MFE/nt" value={a.mfe_per_nt != null ? `${a.mfe_per_nt.toFixed(3)} kcal/mol/nt` : "\u2014"} mono />
      <SubRow label="5\u2032UTR length" value={`${a.utr5_length} nt`} mono />
      <SubRow label="Fold window" value={`${a.fold_window} nt`} mono />
      <SubRow label="Global MFE" value={`${g.mfe.toFixed(1)} kcal/mol`} mono />
      <SubRow label="Global MFE/nt" value={`${g.mfe_per_nt.toFixed(3)} kcal/mol/nt`} mono />
      <div className="mt-2 text-xs" style={{ color: "var(--text-tertiary)" }}>
        {a.message}
      </div>
    </>
  );
}

function ManufacturingDetails({ report }: { report: ScoreReport }) {
  const m = report.manufacturing_scores;
  return (
    <>
      <SubRow label="Total violations" value={String(m.total_violations)} mono />
      <SubRow label="5\u2032UTR violations" value={String(m.utr5_violations)} mono />
      <SubRow label="GC windows" value={m.gc_windows.pass ? "Pass" : `${m.gc_windows.violations.length} violations`} />
      <SubRow label="Homopolymers" value={m.homopolymers.pass ? "Pass" : `${m.homopolymers.violations.length} runs`} />
      <SubRow label="Restriction sites" value={m.restriction_sites.pass ? "Pass" : `${m.restriction_sites.violations.length} sites`} />
      <SubRow label="Upstream ORFs" value={m.uorfs.pass ? "Pass (0)" : `${m.uorfs.count} uORFs found`} />
      {m.uorfs.violations.length > 0 && (
        <div className="mt-1 text-xs" style={{ color: "var(--amber)" }}>
          uORF positions: {m.uorfs.positions.join(", ")}
        </div>
      )}
    </>
  );
}

function StabilityDetails({ report }: { report: ScoreReport }) {
  const s = report.stability_scores;
  return (
    <>
      <SubRow label="GC3 (wobble)" value={formatPercentage(s.gc3)} mono />
      <SubRow label="MFE/nt" value={`${s.mfe_per_nt.toFixed(3)} kcal/mol/nt`} mono />
      <SubRow label="Combined score" value={formatScore(s.stability_score)} mono />
    </>
  );
}

function SpecificityDetails({ report }: { report: ScoreReport }) {
  const r = report.ribonn_scores;
  return (
    <>
      <SubRow label="Target cell type" value={r.target_cell_type} />
      <SubRow label="Target TE" value={formatScore(r.target_te, 4)} mono />
      <SubRow label="Mean TE (all tissues)" value={formatScore(r.mean_te, 4)} mono />
      <SubRow label="Mean off-target TE" value={formatScore(r.mean_off_target_te, 4)} mono />
      <div className="mt-2 text-xs" style={{ color: "var(--text-tertiary)" }}>
        {r.message}
      </div>
    </>
  );
}

const DETAIL_COMPONENTS: Record<string, React.FC<{ report: ScoreReport }>> = {
  utr5_accessibility: Utr5Details,
  manufacturability: ManufacturingDetails,
  stability: StabilityDetails,
  specificity: SpecificityDetails,
};

const PRIORITY_COLOURS: Record<string, string> = {
  high: "var(--red)",
  medium: "var(--amber)",
  low: "var(--text-secondary)",
};

export default function MetricDetail({ metricKey, score, report, suggestions }: MetricDetailProps) {
  const [expanded, setExpanded] = useState(false);
  const DetailComponent = DETAIL_COMPONENTS[metricKey];

  return (
    <div
      className="rounded-xl border overflow-hidden"
      style={{ background: "var(--bg-surface)", borderColor: "var(--border)" }}
    >
      {/* Header */}
      <button
        onClick={() => setExpanded(!expanded)}
        className="w-full flex items-center gap-3 px-5 py-4 transition-colors cursor-pointer"
        onMouseEnter={(e) => (e.currentTarget.style.background = "var(--bg-hover)")}
        onMouseLeave={(e) => (e.currentTarget.style.background = "transparent")}
      >
        <svg
          viewBox="0 0 12 12"
          className={`w-3 h-3 shrink-0 transition-transform duration-200 ${expanded ? "rotate-90" : ""}`}
          fill="currentColor"
          style={{ color: "var(--text-tertiary)" }}
        >
          <path d="M4.5 2l4 4-4 4" />
        </svg>

        <StatusBadge status={score.status as MetricStatus} size="md" />

        <span className="text-sm font-medium flex-1 text-left" style={{ color: "var(--text-primary)" }}>
          {formatMetricName(metricKey)}
        </span>

        <div className="flex items-center gap-3 shrink-0">
          <MetricBar value={score.value} status={score.status as MetricStatus} showValue />
          <span className="text-xs" style={{ color: "var(--text-tertiary)" }}>
            {formatPercentage(score.weight)}
          </span>
        </div>
      </button>

      {/* Expanded detail */}
      {expanded && (
        <div
          className="px-5 pb-4 pt-0"
          style={{ borderTop: "1px solid var(--border)" }}
        >
          <div className="pt-3">
            {DetailComponent && <DetailComponent report={report} />}
          </div>

          {/* Suggestions */}
          {suggestions && suggestions.length > 0 && (
            <div className="mt-3 pt-3" style={{ borderTop: "1px solid var(--border)" }}>
              <span className="text-xs font-medium block mb-2" style={{ color: "var(--text-secondary)" }}>
                Suggestions
              </span>
              {suggestions.map((s, i) => (
                <div key={i} className="flex items-start gap-2 mb-1.5">
                  <span
                    className="text-xs font-medium uppercase shrink-0 mt-0.5"
                    style={{ color: PRIORITY_COLOURS[s.priority] }}
                  >
                    {s.priority}
                  </span>
                  <span className="text-xs" style={{ color: "var(--text-primary)" }}>
                    {s.action}
                  </span>
                </div>
              ))}
            </div>
          )}
        </div>
      )}
    </div>
  );
}
