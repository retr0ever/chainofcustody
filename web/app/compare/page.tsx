"use client";

import { Suspense, useMemo } from "react";
import { useSearchParams } from "next/navigation";
import Link from "next/link";
import { useCandidates } from "@/lib/candidates";
import { METRIC_ORDER } from "@/lib/format";
import { CHART_PALETTE } from "@/lib/colours";
import RadarChart from "../components/RadarChart";
import ComparisonTable from "../components/ComparisonTable";

function CompareContent() {
  const searchParams = useSearchParams();
  const { candidates, loading, compareSet, toggleCompare } = useCandidates();

  const idsFromUrl = useMemo(() => {
    const raw = searchParams.get("ids");
    return raw ? raw.split(",").filter(Boolean) : [];
  }, [searchParams]);

  const selectedIds = useMemo(() => {
    return idsFromUrl.length > 0 ? new Set(idsFromUrl) : compareSet;
  }, [idsFromUrl, compareSet]);

  const selected = useMemo(
    () => candidates.filter((c) => selectedIds.has(c.id)),
    [candidates, selectedIds]
  );

  const radarCandidates = useMemo(
    () =>
      selected.map((c, i) => ({
        id: c.id,
        name: c.name,
        scores: Object.fromEntries(
          METRIC_ORDER.map((m) => [m, c.fitness.scores[m]?.value ?? 0])
        ),
        colour: CHART_PALETTE[i % CHART_PALETTE.length],
      })),
    [selected]
  );

  if (loading) {
    return (
      <div className="flex items-center justify-center h-screen">
        <div className="flex items-center gap-3" style={{ color: "var(--text-secondary)" }}>
          <svg className="animate-spin w-5 h-5" style={{ color: "var(--primary)" }} viewBox="0 0 24 24" fill="none">
            <circle className="opacity-25" cx="12" cy="12" r="10" stroke="currentColor" strokeWidth="4" />
            <path className="opacity-75" fill="currentColor" d="M4 12a8 8 0 018-8v4a4 4 0 00-4 4H4z" />
          </svg>
          <span className="text-sm">Loading...</span>
        </div>
      </div>
    );
  }

  return (
    <div className="p-8 max-w-7xl mx-auto">
      {/* Header */}
      <div className="mb-8">
        <div className="flex items-center gap-2 mb-2 text-xs" style={{ color: "var(--text-tertiary)" }}>
          <Link href="/" className="transition-colors" style={{ color: "var(--text-secondary)" }}
            onMouseEnter={(e) => (e.currentTarget.style.color = "var(--primary)")}
            onMouseLeave={(e) => (e.currentTarget.style.color = "var(--text-secondary)")}
          >
            Dashboard
          </Link>
          <span>/</span>
          <span style={{ color: "var(--text-primary)" }}>Compare</span>
        </div>
        <h1
          className="text-2xl font-semibold tracking-tight"
          style={{ color: "var(--text-primary)" }}
        >
          Candidate Comparison
        </h1>
      </div>

      {/* Candidate selector */}
      <div className="mb-8">
        <div className="flex flex-wrap gap-2">
          {candidates.map((c) => {
            const isSelected = selectedIds.has(c.id);
            return (
              <button
                key={c.id}
                onClick={() => toggleCompare(c.id)}
                className="text-xs font-medium px-3 py-1.5 rounded-lg border transition-all cursor-pointer"
                style={{
                  borderColor: isSelected
                    ? CHART_PALETTE[selected.findIndex((s) => s.id === c.id) % CHART_PALETTE.length]
                    : "var(--border)",
                  background: isSelected ? "var(--bg-raised)" : "transparent",
                  color: isSelected
                    ? CHART_PALETTE[selected.findIndex((s) => s.id === c.id) % CHART_PALETTE.length]
                    : "var(--text-secondary)",
                }}
              >
                {c.name}
              </button>
            );
          })}
        </div>
        <p className="text-xs mt-2" style={{ color: "var(--text-tertiary)" }}>
          Select 2-4 candidates to compare. {selected.length} selected.
        </p>
      </div>

      {selected.length < 2 ? (
        <div
          className="flex flex-col items-center justify-center h-48 rounded-xl border border-dashed"
          style={{ borderColor: "var(--border)", color: "var(--text-tertiary)" }}
        >
          <svg viewBox="0 0 24 24" className="w-8 h-8 mb-2" fill="currentColor" opacity={0.3}>
            <path d="M9 17V7m0 10a2 2 0 01-2 2H5a2 2 0 01-2-2V7a2 2 0 012-2h2a2 2 0 012 2m0 10a2 2 0 002 2h2a2 2 0 002-2M9 7a2 2 0 012-2h2a2 2 0 012 2m0 10V7m0 10a2 2 0 002 2h2a2 2 0 002-2V7a2 2 0 00-2-2h-2a2 2 0 00-2 2" />
          </svg>
          <p className="text-sm">Select at least 2 candidates to compare.</p>
        </div>
      ) : (
        <div className="flex flex-col gap-8">
          {/* Radar chart */}
          <div
            className="rounded-xl border p-6"
            style={{ background: "var(--bg-surface)", borderColor: "var(--border)" }}
          >
            <h2
              className="text-sm font-semibold mb-2 uppercase tracking-wider"
              style={{ color: "var(--text-secondary)" }}
            >
              Metric Overlay
            </h2>
            <RadarChart candidates={radarCandidates} height={400} />
          </div>

          {/* Comparison table */}
          <div>
            <h2
              className="text-sm font-semibold mb-4 uppercase tracking-wider"
              style={{ color: "var(--text-secondary)" }}
            >
              Metric Comparison
            </h2>
            <ComparisonTable candidates={selected} />
          </div>
        </div>
      )}
    </div>
  );
}

export default function ComparePage() {
  return (
    <Suspense
      fallback={
        <div className="flex items-center justify-center h-screen">
          <div className="flex items-center gap-3" style={{ color: "var(--text-secondary)" }}>
            <svg className="animate-spin w-5 h-5" style={{ color: "var(--primary)" }} viewBox="0 0 24 24" fill="none">
              <circle className="opacity-25" cx="12" cy="12" r="10" stroke="currentColor" strokeWidth="4" />
              <path className="opacity-75" fill="currentColor" d="M4 12a8 8 0 018-8v4a4 4 0 00-4 4H4z" />
            </svg>
            <span className="text-sm">Loading...</span>
          </div>
        </div>
      }
    >
      <CompareContent />
    </Suspense>
  );
}
