"use client";

import type { CandidateEvaluation, MetricStatus } from "@/lib/types";
import { METRIC_ORDER, formatMetricName, formatScore } from "@/lib/format";
import StatusBadge from "./StatusBadge";

interface Props {
  candidates: CandidateEvaluation[];
}

export default function ComparisonTable({ candidates }: Props) {
  if (candidates.length === 0) return null;

  const rows = [
    {
      key: "overall",
      label: "Overall Fitness",
      getValue: (c: CandidateEvaluation) => c.fitness.overall,
      getStatus: (c: CandidateEvaluation): MetricStatus =>
        c.fitness.overall >= 0.7 ? "GREEN" : c.fitness.overall >= 0.4 ? "AMBER" : "RED",
      bold: true,
    },
    ...METRIC_ORDER.map((m) => ({
      key: m,
      label: formatMetricName(m),
      getValue: (c: CandidateEvaluation) => c.fitness.scores[m]?.value ?? 0,
      getStatus: (c: CandidateEvaluation): MetricStatus =>
        (c.fitness.scores[m]?.status as MetricStatus) ?? "GREY",
      bold: false,
    })),
  ];

  return (
    <div
      className="rounded-xl border overflow-hidden"
      style={{ background: "var(--bg-surface)", borderColor: "var(--border)" }}
    >
      <div className="overflow-x-auto">
        <table className="w-full">
          <thead>
            <tr style={{ background: "var(--bg-inset)" }}>
              <th
                className="px-5 py-3 text-left text-xs font-medium uppercase tracking-wider"
                style={{ color: "var(--text-secondary)" }}
              >
                Metric
              </th>
              {candidates.map((c) => (
                <th
                  key={c.id}
                  className="px-5 py-3 text-left text-xs font-medium uppercase tracking-wider"
                  style={{ color: "var(--text-primary)" }}
                >
                  {c.name}
                </th>
              ))}
            </tr>
          </thead>
          <tbody>
            {rows.map((row) => {
              const values = candidates.map(row.getValue);
              const best = Math.max(...values);

              return (
                <tr
                  key={row.key}
                  style={{
                    borderTop: "1px solid var(--border)",
                    background: row.bold ? "var(--bg-raised)" : undefined,
                  }}
                >
                  <td className="px-5 py-3">
                    <span
                      className={`text-sm ${row.bold ? "font-semibold" : "font-medium"}`}
                      style={{ color: "var(--text-primary)" }}
                    >
                      {row.label}
                    </span>
                  </td>
                  {candidates.map((c, i) => {
                    const val = values[i];
                    const isWinner = val === best && candidates.length > 1;
                    const status = row.getStatus(c);

                    return (
                      <td key={c.id} className="px-5 py-3">
                        <div className="flex items-center gap-2">
                          <StatusBadge status={status} size="sm" />
                          <span
                            className={`text-sm font-mono tabular-nums ${
                              row.bold ? "font-semibold" : ""
                            }`}
                            style={{
                              color: isWinner ? "var(--green)" : "var(--text-primary)",
                            }}
                          >
                            {formatScore(val)}
                          </span>
                          {isWinner && (
                            <svg
                              viewBox="0 0 12 12"
                              className="w-3 h-3"
                              fill="var(--green)"
                            >
                              <path d="M10 3L4.5 8.5 2 6" stroke="currentColor" strokeWidth="1.5" fill="none" strokeLinecap="round" strokeLinejoin="round" />
                            </svg>
                          )}
                        </div>
                      </td>
                    );
                  })}
                </tr>
              );
            })}
          </tbody>
        </table>
      </div>
    </div>
  );
}
