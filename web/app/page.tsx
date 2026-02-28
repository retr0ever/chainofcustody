"use client";

import { useState, useMemo, useEffect } from "react";
import { useMirnaData } from "@/lib/mirna-data";
import { greedyCover } from "@/lib/greedy";
import { generateSpongeUtr } from "@/lib/sponge";
import type { GreedyResult, GreedyParams } from "@/lib/types";
import CellTypeSelector from "./components/CellTypeSelector";
import ParameterControls from "./components/ParameterControls";
import ResultsTable from "./components/ResultsTable";
import MirnaChart from "./components/MirnaChart";
import SequencePanel from "./components/SequencePanel";
import StructurePlots from "./components/StructurePlots";

function useDebounced<T>(value: T, delayMs: number): T {
  const [debounced, setDebounced] = useState(value);
  useEffect(() => {
    const id = setTimeout(() => setDebounced(value), delayMs);
    return () => clearTimeout(id);
  }, [value, delayMs]);
  return debounced;
}

export default function HomePage() {
  const { data, loading, error } = useMirnaData();

  const [targets, setTargets] = useState<string[]>([]);
  const [offTargets, setOffTargets] = useState<string[]>([]);
  const [targetThreshold, setTargetThreshold] = useState(10);
  const [coverThreshold, setCoverThreshold] = useState(1000);
  const [maxMirnas, setMaxMirnas] = useState(20);

  const params = useMemo<GreedyParams | null>(() => {
    if (!data || targets.length === 0 || offTargets.length === 0) return null;
    return { targets, offTargets, targetThreshold, coverThreshold, maxMirnas };
  }, [data, targets, offTargets, targetThreshold, coverThreshold, maxMirnas]);

  const debouncedParams = useDebounced(params, 150);

  const result = useMemo<GreedyResult | null>(() => {
    if (!data || !debouncedParams) return null;
    return greedyCover(data, debouncedParams);
  }, [data, debouncedParams]);

  // Generate mRNA sponge 3'UTR from selected miRNAs
  const sponge = useMemo(() => {
    if (!result || result.steps.length === 0) return null;
    const stepsWithSeq = result.steps.filter((s) => s.matureSeq.length > 0);
    if (stepsWithSeq.length === 0) return null;
    return generateSpongeUtr(
      stepsWithSeq.map((s) => s.matureSeq),
      stepsWithSeq.map((s) => s.mirnaId),
    );
  }, [result]);

  return (
    <div className="px-3 py-5 sm:px-6 sm:py-6 lg:p-8 max-w-[1600px]">
      {/* Header */}
      <div className="mb-5 sm:mb-6 lg:mb-8 flex flex-col sm:flex-row sm:items-center gap-2 sm:gap-3">
        <div>
          <h1
            className="text-lg sm:text-xl lg:text-2xl font-semibold tracking-tight"
            style={{ color: "var(--text-primary)" }}
          >
            The Optimizer
          </h1>
          <p className="text-xs sm:text-sm mt-1" style={{ color: "var(--text-secondary)" }}>
            Greedy set-cover algorithm for tissue-specific 3&apos;UTR sponge design
          </p>
        </div>
        {data && (
          <span
            className="text-xs rounded-full px-2.5 py-1 font-mono"
            style={{ background: "var(--bg-raised)", color: "var(--text-secondary)" }}
          >
            {data.mirnas.length} miRNAs &middot; {data.cell_types.length} cell types
          </span>
        )}
      </div>

      {loading && (
        <div className="flex items-center justify-center h-64 text-sm gap-3" style={{ color: "var(--text-secondary)" }}>
          <svg className="animate-spin w-5 h-5" style={{ color: "var(--primary)" }} viewBox="0 0 24 24" fill="none">
            <circle className="opacity-25" cx="12" cy="12" r="10" stroke="currentColor" strokeWidth="4" />
            <path className="opacity-75" fill="currentColor" d="M4 12a8 8 0 018-8v4a4 4 0 00-4 4H4z" />
          </svg>
          Loading expression data...
        </div>
      )}

      {error && (
        <div
          className="rounded-lg px-5 py-4 text-sm"
          style={{ background: "var(--red-bg)", color: "var(--red)", border: "1px solid var(--red)" }}
        >
          <strong>Error loading data:</strong> {error}
        </div>
      )}

      {data && (
        <div className="grid grid-cols-1 lg:grid-cols-[320px_1fr] gap-4 sm:gap-6 lg:gap-8 items-start">
          {/* Left sidebar */}
          <aside className="flex flex-col gap-4 sm:gap-6">
            <section
              className="rounded-xl border p-4 sm:p-5"
              style={{ background: "var(--bg-surface)", borderColor: "var(--border)" }}
            >
              <h2
                className="text-sm font-semibold mb-3 sm:mb-4"
                style={{ color: "var(--text-primary)" }}
              >
                Cell type selection
              </h2>
              <CellTypeSelector
                allCellTypes={data.cell_types}
                targets={targets}
                offTargets={offTargets}
                onTargetsChange={setTargets}
                onOffTargetsChange={setOffTargets}
              />
            </section>

            <section
              className="rounded-xl border p-5"
              style={{ background: "var(--bg-surface)", borderColor: "var(--border)" }}
            >
              <h2
                className="text-sm font-semibold mb-4"
                style={{ color: "var(--text-primary)" }}
              >
                Algorithm parameters
              </h2>
              <ParameterControls
                targetThreshold={targetThreshold}
                coverThreshold={coverThreshold}
                maxMirnas={maxMirnas}
                onTargetThresholdChange={setTargetThreshold}
                onCoverThresholdChange={setCoverThreshold}
                onMaxMirnasChange={setMaxMirnas}
              />
            </section>

            <section
              className="rounded-xl border p-4 text-xs space-y-3"
              style={{
                background: "var(--bg-raised)",
                borderColor: "var(--border)",
                color: "var(--text-secondary)",
              }}
            >
              <p className="font-medium" style={{ color: "var(--text-primary)" }}>How it works</p>
              <p>
                Select <span style={{ color: "var(--red)" }}>target</span> cell types
                (to protect — miRNAs must be silent here) and{" "}
                <span style={{ color: "var(--primary)" }}>off-target</span> cell types
                (to suppress).
              </p>
              <p>
                The greedy set-cover algorithm finds the fewest miRNAs that collectively
                silence every off-target while remaining below the threshold in all targets.
                The resulting mRNA sponge 3&apos;UTR contains bulged binding sites for the
                selected miRNAs.
              </p>

              <div style={{ borderTop: "1px solid var(--border)" }} className="pt-3">
                <p className="font-medium mb-2" style={{ color: "var(--text-primary)" }}>Quick start</p>
                <ol className="space-y-1.5 list-none">
                  <li className="flex gap-2">
                    <span className="shrink-0 w-4 h-4 rounded-full flex items-center justify-center text-[10px] font-bold" style={{ background: "var(--red-bg)", color: "var(--red)" }}>1</span>
                    <span>Switch to <span style={{ color: "var(--red)" }}>Targets</span> and click the cell types your therapeutic must spare.</span>
                  </li>
                  <li className="flex gap-2">
                    <span className="shrink-0 w-4 h-4 rounded-full flex items-center justify-center text-[10px] font-bold" style={{ background: "var(--primary-bg)", color: "var(--primary)" }}>2</span>
                    <span>Switch to <span style={{ color: "var(--primary)" }}>Off-targets</span> and select cell types to suppress — or use the quick-select button to add all remaining.</span>
                  </li>
                  <li className="flex gap-2">
                    <span className="shrink-0 w-4 h-4 rounded-full flex items-center justify-center text-[10px] font-bold" style={{ background: "var(--bg-inset)", color: "var(--text-secondary)" }}>3</span>
                    <span>Tune the thresholds and max miRNAs if needed. Results update automatically.</span>
                  </li>
                  <li className="flex gap-2">
                    <span className="shrink-0 w-4 h-4 rounded-full flex items-center justify-center text-[10px] font-bold" style={{ background: "var(--bg-inset)", color: "var(--text-secondary)" }}>4</span>
                    <span>Review the mRNA sponge sequence and structure plots on the right.</span>
                  </li>
                </ol>
              </div>
            </section>
          </aside>

          {/* Right panel */}
          <div className="flex flex-col gap-4 sm:gap-6">
            {targets.length === 0 || offTargets.length === 0 ? (
              <div
                className="flex flex-col items-center justify-center h-48 sm:h-64 rounded-xl border border-dashed px-4 text-center"
                style={{ borderColor: "var(--border)", color: "var(--text-tertiary)" }}
              >
                <svg className="w-8 h-8 mb-2" fill="none" viewBox="0 0 24 24" stroke="currentColor" opacity={0.3}>
                  <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={1.5} d="M9 17V7m0 10a2 2 0 01-2 2H5a2 2 0 01-2-2V7a2 2 0 012-2h2a2 2 0 012 2m0 10a2 2 0 002 2h2a2 2 0 002-2M9 7a2 2 0 012-2h2a2 2 0 012 2m0 10V7m0 10a2 2 0 002 2h2a2 2 0 002-2V7a2 2 0 00-2-2h-2a2 2 0 00-2 2" />
                </svg>
                <p className="text-sm">Select at least one target and one off-target cell type to begin.</p>
              </div>
            ) : result ? (
              <>
                {/* Expression chart — first */}
                <section
                  className="rounded-xl border p-4 sm:p-5"
                  style={{ background: "var(--bg-surface)", borderColor: "var(--border)" }}
                >
                  <h2
                    className="text-sm font-semibold mb-3 sm:mb-4"
                    style={{ color: "var(--text-primary)" }}
                  >
                    Sponge target expression across cell types
                  </h2>
                  <MirnaChart
                    data={data}
                    selectedMirnas={result.selectedMirnas}
                    targets={targets}
                    offTargets={offTargets}
                  />
                </section>

                {/* Selected targets table */}
                <section
                  className="rounded-xl border p-4 sm:p-5"
                  style={{ background: "var(--bg-surface)", borderColor: "var(--border)" }}
                >
                  <h2
                    className="text-sm font-semibold mb-3 sm:mb-4"
                    style={{ color: "var(--text-primary)" }}
                  >
                    Selected sponge targets
                  </h2>
                  <ResultsTable result={result} targets={targets} />
                </section>

                {/* mRNA sponge sequence */}
                {sponge && (
                  <section
                    className="rounded-xl border p-4 sm:p-5"
                    style={{ background: "var(--bg-surface)", borderColor: "var(--border)" }}
                  >
                    <h2
                      className="text-sm font-semibold mb-3 sm:mb-4"
                      style={{ color: "var(--text-primary)" }}
                    >
                      mRNA sponge 3&apos;UTR sequence
                    </h2>
                    <SequencePanel sponge={sponge} />
                  </section>
                )}

                {/* Structure visualisations (client-side) */}
                {sponge && (
                  <StructurePlots
                    sponge={sponge}
                    mirnaNames={result.steps
                      .filter((s) => s.matureSeq.length > 0)
                      .map((s) => s.mirnaId)}
                  />
                )}
              </>
            ) : null}
          </div>
        </div>
      )}
    </div>
  );
}
