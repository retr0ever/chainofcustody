"use client";

import { useMemo, useState } from "react";
import dynamic from "next/dynamic";
import type { MirnaData } from "@/lib/types";

const Plot = dynamic(() => import("react-plotly.js"), {
  ssr: false,
  loading: () => (
    <div className="flex items-center justify-center h-64 text-zinc-400 text-sm gap-2 animate-pulse">
      <svg className="w-4 h-4" viewBox="0 0 24 24" fill="none" stroke="currentColor">
        <rect x="3" y="3" width="4" height="18" rx="1" strokeWidth={2} />
        <rect x="10" y="8" width="4" height="13" rx="1" strokeWidth={2} />
        <rect x="17" y="5" width="4" height="16" rx="1" strokeWidth={2} />
      </svg>
      Rendering chart…
    </div>
  ),
});

const PALETTE = [
  "#3b82f6",
  "#10b981",
  "#f59e0b",
  "#8b5cf6",
  "#ef4444",
  "#06b6d4",
  "#f97316",
  "#84cc16",
  "#ec4899",
  "#14b8a6",
  "#a855f7",
  "#eab308",
  "#6366f1",
  "#22c55e",
  "#f43f5e",
];

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

  const orderedCellTypes = useMemo(() => {
    const targetsSorted = [...targets].sort();
    const offTargetsSorted = [...offTargets].sort();
    return [...targetsSorted, ...offTargetsSorted];
  }, [targets, offTargets]);

  const traces = useMemo(() => {
    return selectedMirnas.map((mirna, idx) => {
      const color = PALETTE[idx % PALETTE.length];
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
          `<b>${mirna}</b><br>` +
          "%{x}<br>" +
          "Mean RPM: %{y:.1f}<extra></extra>",
      };
    });
  }, [selectedMirnas, orderedCellTypes, data.mean_matrix]);

  const nTargets = targets.length;

  const shapes = useMemo(() => {
    if (nTargets === 0 || nTargets >= orderedCellTypes.length) return [];
    return [
      {
        type: "line" as const,
        x0: nTargets - 0.5,
        x1: nTargets - 0.5,
        y0: 0,
        y1: 1,
        yref: "paper" as const,
        line: { color: "#e11d48", width: 2, dash: "dot" as const },
      },
    ];
  }, [nTargets, orderedCellTypes.length]);

  const targetShapes = useMemo(() => {
    if (nTargets === 0) return [];
    return [
      {
        type: "rect" as const,
        x0: -0.5,
        x1: nTargets - 0.5,
        y0: 0,
        y1: 1,
        yref: "paper" as const,
        fillcolor: "rgba(239,68,68,0.06)",
        line: { width: 0 },
      },
    ];
  }, [nTargets]);

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
        font: { size: 11, color: "#e11d48" },
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
        font: { size: 11, color: "#0284c7" },
      });
    }
    return result;
  }, [nTargets, orderedCellTypes.length]);

  if (orderedCellTypes.length === 0 || selectedMirnas.length === 0) {
    return (
      <div className="flex items-center justify-center h-64 rounded-xl border border-dashed border-zinc-200 text-zinc-400 text-sm">
        {orderedCellTypes.length === 0
          ? "Select target and off-target cell types to see the chart."
          : "Run the algorithm to see expression data."}
      </div>
    );
  }

  return (
    <div>
      <div className="flex justify-end mb-2">
        <div className="inline-flex rounded-md border border-zinc-200 overflow-hidden text-xs font-medium">
          <button
            onClick={() => setLogScale(false)}
            className={`px-3 py-1.5 transition-colors ${
              !logScale
                ? "bg-zinc-800 text-white"
                : "bg-white text-zinc-500 hover:bg-zinc-50"
            }`}
          >
            Linear
          </button>
          <button
            onClick={() => setLogScale(true)}
            className={`px-3 py-1.5 border-l border-zinc-200 transition-colors ${
              logScale
                ? "bg-zinc-800 text-white"
                : "bg-white text-zinc-500 hover:bg-zinc-50"
            }`}
          >
            Log
          </button>
        </div>
      </div>

      <Plot
        data={traces}
        layout={{
          barmode: "stack",
          margin: { t: 48, b: 120, l: 72, r: 16 },
          xaxis: {
            tickangle: -45,
            tickfont: { size: 11 },
            automargin: true,
          },
          yaxis: {
            title: { text: "Mean expression (RPM)", font: { size: 12 } },
            tickfont: { size: 11 },
            type: logScale ? "log" : "linear",
          },
          legend: {
            orientation: "h",
            y: -0.35,
            x: 0,
            font: { size: 11 },
          },
          shapes: [...targetShapes, ...shapes],
          annotations,
          plot_bgcolor: "#ffffff",
          paper_bgcolor: "#ffffff",
          font: { family: "ui-sans-serif, system-ui, sans-serif" },
        }}
        config={{ responsive: true, displayModeBar: true, displaylogo: false }}
        style={{ width: "100%", minHeight: 420 }}
        useResizeHandler
      />
    </div>
  );
}
