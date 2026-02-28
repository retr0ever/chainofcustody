"use client";

import { useMemo, useState } from "react";
import dynamic from "next/dynamic";
import type { MirnaData } from "@/lib/types";
import { colours, CHART_PALETTE, PLOTLY_DARK_LAYOUT } from "@/lib/colours";

const Plot = dynamic(() => import("react-plotly.js"), {
  ssr: false,
  loading: () => (
    <div
      className="flex items-center justify-center h-64 text-sm gap-2 animate-pulse"
      style={{ color: "var(--text-tertiary)" }}
    >
      <svg className="w-4 h-4" viewBox="0 0 24 24" fill="none" stroke="currentColor">
        <rect x="3" y="3" width="4" height="18" rx="1" strokeWidth={2} />
        <rect x="10" y="8" width="4" height="13" rx="1" strokeWidth={2} />
        <rect x="17" y="5" width="4" height="16" rx="1" strokeWidth={2} />
      </svg>
      Rendering chart...
    </div>
  ),
});

interface MirnaChartProps {
  data: MirnaData;
  selectedMirnas: string[];
  targets: string[];
  offTargets: string[];
}

export default function MirnaChart({
  data,
  selectedMirnas,
  targets,
  offTargets,
}: MirnaChartProps) {
  const [logScale, setLogScale] = useState(true);

  // Sort cell types by total expression across selected miRNAs (descending)
  // within each group (targets, then off-targets)
  const orderedCellTypes = useMemo(() => {
    function totalExpr(ct: string): number {
      let sum = 0;
      for (const mirna of selectedMirnas) {
        sum += data.mean_matrix[mirna]?.[ct] ?? 0;
      }
      return sum;
    }
    const targetsSorted = [...targets].sort((a, b) => totalExpr(b) - totalExpr(a));
    const offTargetsSorted = [...offTargets].sort((a, b) => totalExpr(b) - totalExpr(a));
    return [...targetsSorted, ...offTargetsSorted];
  }, [targets, offTargets, selectedMirnas, data.mean_matrix]);

  const traces = useMemo(() => {
    return selectedMirnas.map((mirna, idx) => {
      const color = CHART_PALETTE[idx % CHART_PALETTE.length];
      const xLabels = orderedCellTypes.map((ct) => ct.replace(/_/g, " "));
      const yValues = orderedCellTypes.map(
        (ct) => data.mean_matrix[mirna]?.[ct] ?? 0
      );
      return {
        type: "bar" as const,
        name: mirna,
        x: xLabels,
        y: yValues,
        marker: { color },
        hovertemplate:
          `<b>${mirna}</b><br>%{x}<br>Mean RPM: %{y:.1f}<extra></extra>`,
      };
    });
  }, [selectedMirnas, orderedCellTypes, data.mean_matrix]);

  const nTargets = targets.length;

  const shapes = useMemo(() => {
    const result = [];
    if (nTargets > 0) {
      result.push({
        type: "rect" as const,
        x0: -0.5,
        x1: nTargets - 0.5,
        y0: 0,
        y1: 1,
        yref: "paper" as const,
        fillcolor: "rgba(248, 113, 113, 0.06)",
        line: { width: 0 },
      });
    }
    if (nTargets > 0 && nTargets < orderedCellTypes.length) {
      result.push({
        type: "line" as const,
        x0: nTargets - 0.5,
        x1: nTargets - 0.5,
        y0: 0,
        y1: 1,
        yref: "paper" as const,
        line: { color: colours.red, width: 2, dash: "dot" as const },
      });
    }
    return result;
  }, [nTargets, orderedCellTypes.length]);

  const annotations = useMemo(() => {
    const result = [];
    if (nTargets > 0) {
      result.push({
        x: (nTargets - 1) / 2,
        y: 1.04,
        xref: "x" as const,
        yref: "paper" as const,
        text: "Target cells",
        showarrow: false,
        font: { size: 11, color: colours.red },
      });
    }
    if (nTargets < orderedCellTypes.length) {
      result.push({
        x: nTargets + (orderedCellTypes.length - nTargets - 1) / 2,
        y: 1.04,
        xref: "x" as const,
        yref: "paper" as const,
        text: "Off-target cells",
        showarrow: false,
        font: { size: 11, color: colours.primary },
      });
    }
    return result;
  }, [nTargets, orderedCellTypes.length]);

  if (orderedCellTypes.length === 0 || selectedMirnas.length === 0) {
    return (
      <div
        className="flex items-center justify-center h-64 rounded-xl border border-dashed text-sm"
        style={{ borderColor: "var(--border)", color: "var(--text-tertiary)" }}
      >
        {orderedCellTypes.length === 0
          ? "Select target and off-target cell types to see the chart."
          : "Run the algorithm to see expression data."}
      </div>
    );
  }

  // Ensure a minimum width per bar so the chart remains readable on mobile
  const minChartWidth = Math.max(600, orderedCellTypes.length * 18);

  return (
    <div>
      <div className="flex justify-end mb-2">
        <div
          className="inline-flex rounded-md border overflow-hidden text-xs font-medium"
          style={{ borderColor: "var(--border)" }}
        >
          <button
            onClick={() => setLogScale(false)}
            className="px-3 py-2 transition-colors cursor-pointer"
            style={{
              background: !logScale ? "var(--primary)" : "transparent",
              color: !logScale ? "var(--bg-base)" : "var(--text-secondary)",
            }}
          >
            Linear
          </button>
          <button
            onClick={() => setLogScale(true)}
            className="px-3 py-2 transition-colors cursor-pointer"
            style={{
              background: logScale ? "var(--primary)" : "transparent",
              color: logScale ? "var(--bg-base)" : "var(--text-secondary)",
              borderLeft: "1px solid var(--border)",
            }}
          >
            Log
          </button>
        </div>
      </div>

      <div className="overflow-x-auto -mx-5 px-5 sm:mx-0 sm:px-0">
        <div style={{ minWidth: minChartWidth }}>
          <Plot
            data={traces}
            layout={{
              ...PLOTLY_DARK_LAYOUT,
              barmode: "stack",
              margin: { t: 48, b: 120, l: 72, r: 16 },
              xaxis: {
                tickangle: -45,
                tickfont: { size: 11, color: colours.textSecondary },
                gridcolor: colours.border,
                automargin: true,
              },
              yaxis: {
                title: { text: "Mean expression (RPM)", font: { size: 12, color: colours.textSecondary } },
                tickfont: { size: 11, color: colours.textSecondary },
                type: logScale ? "log" : "linear",
                gridcolor: colours.border,
              },
              legend: {
                orientation: "h" as const,
                y: -0.35,
                x: 0,
                font: { size: 11, color: colours.textSecondary },
              },
              shapes,
              annotations,
            }}
            config={{ responsive: true, displayModeBar: false, displaylogo: false }}
            style={{ width: "100%", minHeight: 380 }}
            useResizeHandler
          />
        </div>
      </div>
    </div>
  );
}
