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
        className={`rounded-lg px-4 py-3 text-sm flex flex-wrap gap-4 items-center ${
          success
            ? "bg-emerald-50 border border-emerald-200 text-emerald-800"
            : "bg-amber-50 border border-amber-200 text-amber-800"
        }`}
      >
        <span>
          <span className="font-semibold">{success ? "Full coverage" : "Partial coverage"}</span>
          {" — "}
          {allOffTargets.length - uncovered.length}/{allOffTargets.length} off-targets covered (
          {coverPct}%)
        </span>
        <span>
          <span className="font-semibold">{selectedMirnas.length}</span> miRNA
          {selectedMirnas.length !== 1 ? "s" : ""} selected
        </span>
        {targets.length > 0 && (
          <span>
            Protecting{" "}
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
            <p className="mt-1 text-xs font-mono">
              {uncovered.map((c) => c.replace(/_/g, " ")).join(", ")}
            </p>
          </details>
        )}
      </div>

      {/* miRNA table */}
      {steps.length > 0 ? (
        <div className="overflow-x-auto rounded-lg border border-zinc-200">
          <table className="w-full text-sm">
            <thead className="bg-zinc-50 border-b border-zinc-200">
              <tr>
                <th className="px-3 py-2 text-left font-medium text-zinc-500 text-xs">#</th>
                <th className="px-3 py-2 text-left font-medium text-zinc-500 text-xs">MiRBase ID</th>
                <th className="px-3 py-2 text-left font-medium text-zinc-500 text-xs">Seed</th>
                <th className="px-3 py-2 text-right font-medium text-zinc-500 text-xs">
                  Target RPM (avg)
                </th>
                <th className="px-3 py-2 text-right font-medium text-zinc-500 text-xs">
                  Newly covered
                </th>
                <th className="px-3 py-2 text-left font-medium text-zinc-500 text-xs hidden lg:table-cell">
                  Mature sequence
                </th>
              </tr>
            </thead>
            <tbody className="divide-y divide-zinc-100">
              {steps.map((step, i) => (
                <tr key={step.mirnaId} className="hover:bg-zinc-50 transition-colors">
                  <td className="px-3 py-2 text-zinc-400 font-mono text-xs">{i + 1}</td>
                  <td className="px-3 py-2 font-mono font-medium text-zinc-800 text-xs">
                    {step.mirnaId}
                  </td>
                  <td className="px-3 py-2 font-mono text-xs text-zinc-600">{step.seed || "—"}</td>
                  <td className="px-3 py-2 text-right font-mono text-xs text-zinc-600">
                    {step.targetRpm.toFixed(1)}
                  </td>
                  <td className="px-3 py-2 text-right">
                    <span className="inline-block rounded-full bg-sky-100 text-sky-700 text-xs px-2 py-0.5 font-medium">
                      {step.newlyCovered.length}
                    </span>
                  </td>
                  <td className="px-3 py-2 font-mono text-xs text-zinc-400 hidden lg:table-cell max-w-xs truncate">
                    {step.matureSeq || "—"}
                  </td>
                </tr>
              ))}
            </tbody>
          </table>
        </div>
      ) : (
        <p className="text-sm text-zinc-400 text-center py-6">
          No miRNAs found — try lowering the silence or coverage thresholds.
        </p>
      )}
    </div>
  );
}
