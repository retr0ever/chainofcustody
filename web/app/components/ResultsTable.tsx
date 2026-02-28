"use client";

import type { GreedyResult } from "@/lib/types";

interface ResultsTableProps {
  result: GreedyResult;
  targets: string[];
}

export default function ResultsTable({ result, targets }: ResultsTableProps) {
  const { selectedMirnas, steps, uncovered, allOffTargets, success } = result;
  const coverPct =
    allOffTargets.length > 0
      ? (((allOffTargets.length - uncovered.length) / allOffTargets.length) * 100).toFixed(1)
      : "0";

  return (
    <div className="flex flex-col gap-4">
      {/* Summary banner */}
      <div
        className="rounded-lg px-4 py-3 text-sm flex flex-wrap gap-4 items-center border"
        style={{
          background: success ? "var(--green-bg)" : "var(--amber-bg)",
          borderColor: success ? "var(--green)" : "var(--amber)",
          color: success ? "var(--green)" : "var(--amber)",
        }}
      >
        <span>
          <span className="font-semibold">{success ? "Full coverage" : "Partial coverage"}</span>
          {" \u2014 "}
          {allOffTargets.length - uncovered.length}/{allOffTargets.length} off-targets covered (
          {coverPct}%)
        </span>
        <span>
          <span className="font-semibold">{selectedMirnas.length}</span> miRNA
          {selectedMirnas.length !== 1 ? "s" : ""} selected
        </span>
        {targets.length > 0 && (
          <span>
            Targeting{" "}
            <span className="font-semibold">
              {targets.map((t) => t.replace(/_/g, " ")).join(", ")}
            </span>
          </span>
        )}
        {uncovered.length > 0 && (
          <details className="w-full">
            <summary className="cursor-pointer font-medium">
              {uncovered.length} uncovered off-target{uncovered.length !== 1 ? "s" : ""}
            </summary>
            <p className="mt-1 text-xs font-mono" style={{ color: "var(--text-secondary)" }}>
              {uncovered.map((c) => c.replace(/_/g, " ")).join(", ")}
            </p>
          </details>
        )}
      </div>

      {/* miRNA table */}
      {steps.length > 0 ? (
        <div
          className="overflow-x-auto rounded-lg border"
          style={{ borderColor: "var(--border)" }}
        >
          <table className="w-full text-sm" style={{ minWidth: 480 }}>
            <thead style={{ background: "var(--bg-inset)" }}>
              <tr style={{ borderBottom: "1px solid var(--border)" }}>
                <th className="px-2 sm:px-3 py-2 text-left font-medium text-xs" style={{ color: "var(--text-secondary)" }}>#</th>
                <th className="px-2 sm:px-3 py-2 text-left font-medium text-xs" style={{ color: "var(--text-secondary)" }}>MiRBase ID</th>
                <th className="px-2 sm:px-3 py-2 text-left font-medium text-xs hidden sm:table-cell" style={{ color: "var(--text-secondary)" }}>Seed</th>
                <th className="px-2 sm:px-3 py-2 text-right font-medium text-xs" style={{ color: "var(--text-secondary)" }}>Target RPM</th>
                <th className="px-2 sm:px-3 py-2 text-right font-medium text-xs" style={{ color: "var(--text-secondary)" }}>Covered</th>
                <th className="px-2 sm:px-3 py-2 text-left font-medium text-xs hidden lg:table-cell" style={{ color: "var(--text-secondary)" }}>Mature sequence</th>
              </tr>
            </thead>
            <tbody>
              {steps.map((step, i) => (
                <tr
                  key={step.mirnaId}
                  className="transition-colors"
                  style={{ borderBottom: "1px solid var(--border)" }}
                  onMouseEnter={(e) => (e.currentTarget.style.background = "var(--bg-hover)")}
                  onMouseLeave={(e) => (e.currentTarget.style.background = "transparent")}
                >
                  <td className="px-2 sm:px-3 py-2.5 sm:py-2 font-mono text-xs" style={{ color: "var(--text-tertiary)" }}>{i + 1}</td>
                  <td className="px-2 sm:px-3 py-2.5 sm:py-2 font-mono font-medium text-xs whitespace-nowrap" style={{ color: "var(--text-primary)" }}>{step.mirnaId}</td>
                  <td className="px-2 sm:px-3 py-2.5 sm:py-2 font-mono text-xs hidden sm:table-cell" style={{ color: "var(--text-secondary)" }}>{step.seed || "\u2014"}</td>
                  <td className="px-2 sm:px-3 py-2.5 sm:py-2 text-right font-mono text-xs" style={{ color: "var(--text-secondary)" }}>{step.targetRpm.toFixed(1)}</td>
                  <td className="px-2 sm:px-3 py-2.5 sm:py-2 text-right">
                    <span
                      className="inline-block rounded-full text-xs px-2 py-0.5 font-medium"
                      style={{ background: "var(--primary-bg)", color: "var(--primary)" }}
                    >
                      {step.newlyCovered.length}
                    </span>
                  </td>
                  <td className="px-2 sm:px-3 py-2.5 sm:py-2 font-mono text-xs hidden lg:table-cell max-w-xs truncate" style={{ color: "var(--text-tertiary)" }}>
                    {step.matureSeq || "\u2014"}
                  </td>
                </tr>
              ))}
            </tbody>
          </table>
        </div>
      ) : (
        <p className="text-sm text-center py-6" style={{ color: "var(--text-tertiary)" }}>
          No miRNAs found — try lowering the silence or coverage thresholds.
        </p>
      )}
    </div>
  );
}
