"use client";

import { useState, useMemo, useEffect } from "react";
import { useMirnaData } from "@/lib/mirna-data";
import { greedyCover } from "@/lib/greedy";
import { generateUtr } from "@/lib/utr-design";
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
  const [utr5Seq, setUtr5Seq] = useState("");
  const [cdsSeq, setCdsSeq] = useState("");
  const [geneSymbol, setGeneSymbol] = useState("");
  const [fetchingCds, setFetchingCds] = useState(false);
  const [optimizationMode, setOptimizationMode] = useState<"manual" | "optimize">("manual");
  const [optimizing, setOptimizing] = useState(false);

  const handleFetchCds = async () => {
    if (!geneSymbol) return;
    setFetchingCds(true);
    try {
      const apiBase = process.env.NEXT_PUBLIC_STRUCTURE_API ?? "http://localhost:8000";
      const response = await fetch(`${apiBase}/api/gene-cds/${geneSymbol}`);
      const data = await response.json();
      if (data.ok) {
        setCdsSeq(data.cds);
      } else {
        alert(`Error fetching CDS: ${data.error}`);
      }
    } catch (err) {
      alert("Failed to connect to structure API server.");
    } finally {
      setFetchingCds(false);
    }
  };

  const handleRunOptimization = async () => {
    if (!geneSymbol) return;
    setOptimizing(true);
    try {
      const apiBase = process.env.NEXT_PUBLIC_STRUCTURE_API ?? "http://localhost:8000";
      const response = await fetch(`${apiBase}/api/optimize`, {
        method: "POST",
        headers: { "Content-Type": "application/json" },
        body: JSON.stringify({ 
          gene: geneSymbol,
          target_cell_type: targets[0] || "Dendritic_cell",
          n_gen: 10
        }),
      });
      const data = await response.json();
      if (data.ok && data.best?.utr5) {
        setUtr5Seq(data.best.utr5);
        if (data.best.cds) setCdsSeq(data.best.cds);
      } else {
        alert(`Optimisation failed: ${data.error || "Unknown error"}`);
      }
    } catch (err) {
      alert("Failed to connect to structure API server.");
    } finally {
      setOptimizing(false);
    }
  };

  const params = useMemo<GreedyParams | null>(() => {
    if (!data || targets.length === 0 || offTargets.length === 0) return null;
    return { targets, offTargets, targetThreshold, coverThreshold, maxMirnas };
  }, [data, targets, offTargets, targetThreshold, coverThreshold, maxMirnas]);

  const debouncedParams = useDebounced(params, 150);

  const result = useMemo<GreedyResult | null>(() => {
    if (!data || !debouncedParams) return null;
    return greedyCover(data, debouncedParams);
  }, [data, debouncedParams]);

  // Generate 3'UTR design from selected miRNAs
  const design = useMemo(() => {
    if (!result || result.steps.length === 0) return null;
    const stepsWithSeq = result.steps.filter((s) => s.matureSeq.length > 0);
    if (stepsWithSeq.length === 0) return null;
    return generateUtr(
      stepsWithSeq.map((s) => s.matureSeq),
      stepsWithSeq.map((s) => s.mirnaId),
      16,
      utr5Seq,
      cdsSeq
    );
  }, [result, utr5Seq, cdsSeq]);

  return (
    <div className="px-3 py-5 sm:px-6 sm:py-6 lg:p-8 max-w-[1800px]">
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
            Greedy set-cover algorithm for tissue-specific 3&apos;UTR design
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
        <div className="grid grid-cols-1 lg:grid-cols-[320px_1fr] xl:grid-cols-[320px_1fr_480px] gap-4 sm:gap-6 lg:gap-8 items-start">
          {/* Column 1: Selection & Controls */}
          <aside className="flex flex-col gap-4 sm:gap-6 lg:sticky lg:top-8 max-h-[calc(100vh-4rem)] overflow-y-auto pr-1">
            <section
              className="rounded-xl border p-4 sm:p-5 shadow-sm transition-shadow hover:shadow-md"
              style={{ background: "var(--bg-surface)", borderColor: "var(--border)" }}
            >
              <h2
                className="text-sm font-semibold mb-3 sm:mb-4 flex items-center gap-2"
                style={{ color: "var(--text-primary)" }}
              >
                <div className="w-1 h-4 bg-[var(--primary)] rounded-full" />
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
              className="rounded-xl border p-5 shadow-sm transition-shadow hover:shadow-md"
              style={{ background: "var(--bg-surface)", borderColor: "var(--border)" }}
            >
              <h2
                className="text-sm font-semibold mb-4 flex items-center gap-2"
                style={{ color: "var(--text-primary)" }}
              >
                <div className="w-1 h-4 bg-[var(--primary)] rounded-full" />
                Construct design
              </h2>
              <ParameterControls
                targetThreshold={targetThreshold}
                coverThreshold={coverThreshold}
                maxMirnas={maxMirnas}
                utr5Seq={utr5Seq}
                cdsSeq={cdsSeq}
                geneSymbol={geneSymbol}
                fetchingCds={fetchingCds}
                optimizationMode={optimizationMode}
                optimizing={optimizing}
                onTargetThresholdChange={setTargetThreshold}
                onCoverThresholdChange={setCoverThreshold}
                onMaxMirnasChange={setMaxMirnas}
                onUtr5SeqChange={setUtr5Seq}
                onCdsSeqChange={setCdsSeq}
                onGeneSymbolChange={setGeneSymbol}
                onFetchCds={handleFetchCds}
                onOptimizationModeChange={setOptimizationMode}
                onRunOptimization={handleRunOptimization}
              />
            </section>
          </aside>

          {/* Column 2: Expression & miRNA Discovery */}
          <div className="flex flex-col gap-4 sm:gap-6 min-w-0">
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
                  className="rounded-xl border p-4 sm:p-5 shadow-sm"
                  style={{ background: "var(--bg-surface)", borderColor: "var(--border)" }}
                >
                  <h2
                    className="text-sm font-semibold mb-3 sm:mb-4 flex items-center justify-between"
                    style={{ color: "var(--text-primary)" }}
                  >
                    <div className="flex items-center gap-2">
                      <div className="w-1 h-4 bg-[var(--primary)] rounded-full" />
                      Differential expression
                    </div>
                    <span className="text-[10px] font-mono text-[var(--text-tertiary)] bg-[var(--bg-inset)] px-2 py-0.5 rounded">
                      Top candidates
                    </span>
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
                  className="rounded-xl border p-4 sm:p-5 shadow-sm"
                  style={{ background: "var(--bg-surface)", borderColor: "var(--border)" }}
                >
                  <h2
                    className="text-sm font-semibold mb-3 sm:mb-4 flex items-center gap-2"
                    style={{ color: "var(--text-primary)" }}
                  >
                    <div className="w-1 h-4 bg-[var(--primary)] rounded-full" />
                    miRNA sponge logic
                  </h2>
                  <ResultsTable result={result} targets={targets} />
                </section>

                {/* Mobile views for sequence/structure */}
                <div className="xl:hidden flex flex-col gap-4 sm:gap-6">
                  {design && (
                    <section
                      className="rounded-xl border p-4 sm:p-5"
                      style={{ background: "var(--bg-surface)", borderColor: "var(--border)" }}
                    >
                      <h2 className="text-sm font-semibold mb-3 sm:mb-4" style={{ color: "var(--text-primary)" }}>Construct sequence</h2>
                      <SequencePanel design={design} />
                    </section>
                  )}
                  {design && (
                    <StructurePlots
                      design={design}
                      mirnaNames={result.steps.filter((s) => s.matureSeq.length > 0).map((s) => s.mirnaId)}
                    />
                  )}
                </div>
              </>
            ) : null}
          </div>

          {/* Column 3: Sequence & Folding (Desktop Sticky) */}
          <aside className="hidden xl:flex flex-col gap-4 sm:gap-6 lg:sticky lg:top-8 max-h-[calc(100vh-4rem)] overflow-y-auto pr-1 custom-scrollbar">
            {design && (
              <>
                <section
                  className="rounded-xl border p-4 sm:p-5 shadow-sm"
                  style={{ background: "var(--bg-surface)", borderColor: "var(--border)" }}
                >
                  <h2
                    className="text-sm font-semibold mb-3 sm:mb-4 flex items-center gap-2"
                    style={{ color: "var(--text-primary)" }}
                  >
                    <div className="w-1 h-4 bg-[var(--primary)] rounded-full" />
                    mRNA sequence
                  </h2>
                  <SequencePanel design={design} />
                </section>

                <StructurePlots
                  design={design}
                  mirnaNames={result.steps
                    .filter((s) => s.matureSeq.length > 0)
                    .map((s) => s.mirnaId)}
                />
              </>
            )}
          </aside>
        </div>
      )}
    </div>
  );
}
