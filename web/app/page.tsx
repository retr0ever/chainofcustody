"use client";

import { useState, useMemo, useEffect } from "react";
import { useMirnaData } from "@/lib/data";
import { greedyCover } from "@/lib/greedy";
import type { GreedyResult, GreedyParams } from "@/lib/types";
import CellTypeSelector from "./components/CellTypeSelector";
import ParameterControls from "./components/ParameterControls";
import ResultsTable from "./components/ResultsTable";
import MirnaChart from "./components/MirnaChart";

function useDebounced<T>(value: T, delayMs: number): T {
  const [debounced, setDebounced] = useState(value);
  useEffect(() => {
    const id = setTimeout(() => setDebounced(value), delayMs);
    return () => clearTimeout(id);
  }, [value, delayMs]);
  return debounced;
}

export default function Home() {
  const { data, loading, error } = useMirnaData();

  const [targets, setTargets] = useState<string[]>([]);
  const [offTargets, setOffTargets] = useState<string[]>([]);
  const [targetThreshold, setTargetThreshold] = useState(10);
  const [coverThreshold, setCoverThreshold] = useState(1000);
  const [maxMirnas, setMaxMirnas] = useState(20);

  // Debounce the algorithm params so slider dragging doesn't stall the UI
  const params = useMemo<GreedyParams | null>(() => {
    if (!data || targets.length === 0 || offTargets.length === 0) return null;
    return { targets, offTargets, targetThreshold, coverThreshold, maxMirnas };
  }, [data, targets, offTargets, targetThreshold, coverThreshold, maxMirnas]);

  const debouncedParams = useDebounced(params, 150);

  const result = useMemo<GreedyResult | null>(() => {
    if (!data || !debouncedParams) return null;
    return greedyCover(data, debouncedParams);
  }, [data, debouncedParams]);

  return (
    <div className="min-h-screen bg-zinc-50 font-sans">
      {/* Header */}
      <header className="bg-white border-b border-zinc-200 px-6 py-4 sticky top-0 z-10">
        <div className="max-w-screen-xl mx-auto flex items-center gap-3">
          <div className="flex items-center gap-2">
            <span className="text-lg font-bold tracking-tight text-zinc-900">miRNA</span>
            <span className="text-lg font-light text-zinc-400">cell-type selector</span>
          </div>
          {data && (
            <span className="text-xs text-zinc-400 bg-zinc-100 rounded-full px-2.5 py-1 font-mono">
              {data.mirnas.length} miRNAs · {data.cell_types.length} cell types
            </span>
          )}
        </div>
      </header>

      <main className="max-w-screen-xl mx-auto px-6 py-8">
        {loading && (
          <div className="flex items-center justify-center h-64 text-zinc-500 text-sm gap-3">
            <svg
              className="animate-spin w-5 h-5 text-sky-500"
              viewBox="0 0 24 24"
              fill="none"
            >
              <circle
                className="opacity-25"
                cx="12"
                cy="12"
                r="10"
                stroke="currentColor"
                strokeWidth="4"
              />
              <path
                className="opacity-75"
                fill="currentColor"
                d="M4 12a8 8 0 018-8v4a4 4 0 00-4 4H4z"
              />
            </svg>
            Loading expression data…
          </div>
        )}

        {error && (
          <div className="rounded-lg bg-red-50 border border-red-200 text-red-700 px-5 py-4 text-sm">
            <strong>Error loading data:</strong> {error}
          </div>
        )}

        {data && (
          <div className="grid grid-cols-1 lg:grid-cols-[320px_1fr] gap-8 items-start">
            {/* Left sidebar: controls */}
            <aside className="flex flex-col gap-6">
              <section className="bg-white rounded-xl border border-zinc-200 p-5 shadow-sm">
                <h2 className="text-sm font-semibold text-zinc-800 mb-4">Cell type selection</h2>
                <CellTypeSelector
                  allCellTypes={data.cell_types}
                  targets={targets}
                  offTargets={offTargets}
                  onTargetsChange={setTargets}
                  onOffTargetsChange={setOffTargets}
                />
              </section>

              <section className="bg-white rounded-xl border border-zinc-200 p-5 shadow-sm">
                <h2 className="text-sm font-semibold text-zinc-800 mb-4">Algorithm parameters</h2>
                <ParameterControls
                  targetThreshold={targetThreshold}
                  coverThreshold={coverThreshold}
                  maxMirnas={maxMirnas}
                  onTargetThresholdChange={setTargetThreshold}
                  onCoverThresholdChange={setCoverThreshold}
                  onMaxMirnasChange={setMaxMirnas}
                />
              </section>

              <section className="bg-zinc-50 rounded-xl border border-zinc-200 p-4 text-xs text-zinc-500 space-y-1.5">
                <p className="font-medium text-zinc-600">How it works</p>
                <p>
                  Select <span className="text-rose-600 font-medium">target</span> cell types (to
                  protect — miRNAs must be silent here) and{" "}
                  <span className="text-sky-600 font-medium">off-target</span> cell types (to
                  suppress).
                </p>
                <p>
                  The greedy set-cover algorithm finds the fewest miRNAs that collectively silence
                  every off-target while remaining below the threshold in all targets.
                </p>
              </section>
            </aside>

            {/* Right panel: results */}
            <div className="flex flex-col gap-6">
              {/* Prompt if nothing selected */}
              {targets.length === 0 || offTargets.length === 0 ? (
                <div className="flex flex-col items-center justify-center h-64 rounded-xl border border-dashed border-zinc-200 bg-white text-zinc-400 text-sm gap-2">
                  <svg
                    className="w-8 h-8 text-zinc-300"
                    fill="none"
                    viewBox="0 0 24 24"
                    stroke="currentColor"
                  >
                    <path
                      strokeLinecap="round"
                      strokeLinejoin="round"
                      strokeWidth={1.5}
                      d="M9 17V7m0 10a2 2 0 01-2 2H5a2 2 0 01-2-2V7a2 2 0 012-2h2a2 2 0 012 2m0 10a2 2 0 002 2h2a2 2 0 002-2M9 7a2 2 0 012-2h2a2 2 0 012 2m0 10V7m0 10a2 2 0 002 2h2a2 2 0 002-2V7a2 2 0 00-2-2h-2a2 2 0 00-2 2"
                    />
                  </svg>
                  <p>Select at least one target and one off-target cell type to begin.</p>
                </div>
              ) : result ? (
                <>
                  {/* Chart */}
                  <section className="bg-white rounded-xl border border-zinc-200 p-5 shadow-sm">
                    <h2 className="text-sm font-semibold text-zinc-800 mb-4">
                      miRNA expression across selected cell types
                    </h2>
                    <MirnaChart
                      data={data}
                      selectedMirnas={result.selectedMirnas}
                      targets={targets}
                      offTargets={offTargets}
                    />
                  </section>

                  {/* Results table */}
                  <section className="bg-white rounded-xl border border-zinc-200 p-5 shadow-sm">
                    <h2 className="text-sm font-semibold text-zinc-800 mb-4">
                      Selected miRNAs
                    </h2>
                    <ResultsTable result={result} targets={targets} />
                  </section>
                </>
              ) : null}
            </div>
          </div>
        )}
      </main>
    </div>
  );
}
