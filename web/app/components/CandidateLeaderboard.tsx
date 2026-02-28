"use client";

import { useState, useMemo } from "react";
import Link from "next/link";
import type { CandidateEvaluation, MetricStatus } from "@/lib/types";
import { METRIC_ORDER, formatMetricShort, formatScore } from "@/lib/format";
import StatusBadge from "./StatusBadge";
import MetricBar from "./MetricBar";

type SortKey = "overall" | typeof METRIC_ORDER[number];

interface Props {
  candidates: CandidateEvaluation[];
  compareSet: Set<string>;
  onToggleCompare: (id: string) => void;
}

export default function CandidateLeaderboard({ candidates, compareSet, onToggleCompare }: Props) {
  const [sortKey, setSortKey] = useState<SortKey>("overall");
  const [sortDesc, setSortDesc] = useState(true);

  const sorted = useMemo(() => {
    const arr = [...candidates];
    arr.sort((a, b) => {
      let va: number, vb: number;
      if (sortKey === "overall") {
        va = a.fitness.overall;
        vb = b.fitness.overall;
      } else {
        va = a.fitness.scores[sortKey]?.value ?? 0;
        vb = b.fitness.scores[sortKey]?.value ?? 0;
      }
      return sortDesc ? vb - va : va - vb;
    });
    return arr;
  }, [candidates, sortKey, sortDesc]);

  const handleSort = (key: SortKey) => {
    if (sortKey === key) {
      setSortDesc(!sortDesc);
    } else {
      setSortKey(key);
      setSortDesc(true);
    }
  };

  const SortHeader = ({ k, label, className }: { k: SortKey; label: string; className?: string }) => (
    <th
      className={`px-3 py-2.5 text-left text-xs font-medium uppercase tracking-wider cursor-pointer select-none transition-colors ${className ?? ""}`}
      style={{ color: sortKey === k ? "var(--primary)" : "var(--text-secondary)" }}
      onClick={() => handleSort(k)}
    >
      <span className="flex items-center gap-1">
        {label}
        {sortKey === k && (
          <svg viewBox="0 0 10 10" className="w-2.5 h-2.5" fill="currentColor">
            {sortDesc ? (
              <path d="M5 7L1 3h8z" />
            ) : (
              <path d="M5 3L1 7h8z" />
            )}
          </svg>
        )}
      </span>
    </th>
  );

  return (
    <div
      className="rounded-xl border overflow-hidden"
      style={{ background: "var(--bg-surface)", borderColor: "var(--border)" }}
    >
      <div className="overflow-x-auto">
        <table className="w-full">
          <thead>
            <tr style={{ background: "var(--bg-inset)" }}>
              <th className="w-10 px-3 py-2.5" />
              <th className="w-8 px-2 py-2.5">
                <span className="text-xs" style={{ color: "var(--text-tertiary)" }}>#</span>
              </th>
              <SortHeader k="overall" label="Candidate" className="min-w-[140px]" />
              <SortHeader k="overall" label="Fitness" className="min-w-[120px]" />
              {METRIC_ORDER.map((m) => (
                <SortHeader key={m} k={m} label={formatMetricShort(m)} className="min-w-[120px]" />
              ))}
            </tr>
          </thead>
          <tbody>
            {sorted.map((c, i) => (
              <tr
                key={c.id}
                className="transition-colors"
                style={{ borderTop: "1px solid var(--border)" }}
                onMouseEnter={(e) => (e.currentTarget.style.background = "var(--bg-hover)")}
                onMouseLeave={(e) => (e.currentTarget.style.background = "transparent")}
              >
                {/* Checkbox */}
                <td className="px-3 py-3">
                  <label className="flex items-center justify-center cursor-pointer">
                    <input
                      type="checkbox"
                      checked={compareSet.has(c.id)}
                      onChange={() => onToggleCompare(c.id)}
                      className="sr-only peer"
                    />
                    <span
                      className="w-4 h-4 rounded border flex items-center justify-center transition-colors peer-checked:border-transparent"
                      style={{
                        borderColor: compareSet.has(c.id) ? "var(--primary)" : "var(--border-strong)",
                        background: compareSet.has(c.id) ? "var(--primary)" : "transparent",
                      }}
                    >
                      {compareSet.has(c.id) && (
                        <svg viewBox="0 0 12 12" className="w-2.5 h-2.5" fill="var(--bg-base)">
                          <path d="M10 3L4.5 8.5 2 6" stroke="currentColor" strokeWidth="2" fill="none" strokeLinecap="round" strokeLinejoin="round" />
                        </svg>
                      )}
                    </span>
                  </label>
                </td>

                {/* Rank */}
                <td className="px-2 py-3">
                  <span
                    className="text-xs font-mono tabular-nums"
                    style={{ color: "var(--text-tertiary)" }}
                  >
                    {i + 1}
                  </span>
                </td>

                {/* Name */}
                <td className="px-3 py-3">
                  <Link
                    href={`/candidate/${c.id}`}
                    className="text-sm font-medium transition-colors"
                    style={{ color: "var(--text-primary)" }}
                    onMouseEnter={(e) => (e.currentTarget.style.color = "var(--primary)")}
                    onMouseLeave={(e) => (e.currentTarget.style.color = "var(--text-primary)")}
                  >
                    {c.name}
                  </Link>
                  <div className="text-xs mt-0.5" style={{ color: "var(--text-tertiary)" }}>
                    {c.gene} / {c.target_cell_type}
                  </div>
                </td>

                {/* Overall fitness */}
                <td className="px-3 py-3">
                  <div className="flex items-center gap-2">
                    <span
                      className="text-sm font-mono font-semibold tabular-nums"
                      style={{ color: "var(--text-primary)" }}
                    >
                      {formatScore(c.fitness.overall)}
                    </span>
                    <MetricBar
                      value={c.fitness.overall}
                      status={c.fitness.overall >= 0.7 ? "GREEN" : c.fitness.overall >= 0.4 ? "AMBER" : "RED"}
                      compact
                    />
                  </div>
                </td>

                {/* Per-metric columns */}
                {METRIC_ORDER.map((m) => {
                  const s = c.fitness.scores[m];
                  if (!s) return <td key={m} />;
                  return (
                    <td key={m} className="px-3 py-3">
                      <div className="flex items-center gap-2">
                        <StatusBadge status={s.status as MetricStatus} size="sm" />
                        <MetricBar value={s.value} status={s.status as MetricStatus} compact />
                        <span
                          className="text-xs font-mono tabular-nums"
                          style={{ color: "var(--text-secondary)" }}
                        >
                          {formatScore(s.value)}
                        </span>
                      </div>
                    </td>
                  );
                })}
              </tr>
            ))}
          </tbody>
        </table>
      </div>
    </div>
  );
}
