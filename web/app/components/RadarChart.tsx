"use client";

import dynamic from "next/dynamic";
import { useMemo } from "react";
import { colours, PLOTLY_DARK_LAYOUT, CHART_PALETTE } from "@/lib/colours";
import { formatMetricName } from "@/lib/format";

const Plot = dynamic(() => import("react-plotly.js"), { ssr: false });

interface RadarCandidate {
  id: string;
  name: string;
  scores: Record<string, number>;
  colour?: string;
}

interface RadarChartProps {
  candidates: RadarCandidate[];
  metrics?: string[];
  height?: number;
}

const DEFAULT_METRICS = [
  "utr5_accessibility",
  "manufacturability",
  "stability",
  "specificity",
];

export default function RadarChart({
  candidates,
  metrics = DEFAULT_METRICS,
  height = 340,
}: RadarChartProps) {
  const labels = useMemo(() => metrics.map(formatMetricName), [metrics]);

  const traces = useMemo(
    () =>
      candidates.map((c, i) => ({
        type: "scatterpolar" as const,
        r: [...metrics.map((m) => c.scores[m] ?? 0), c.scores[metrics[0]] ?? 0],
        theta: [...labels, labels[0]],
        fill: "toself" as const,
        fillcolor: `${c.colour ?? CHART_PALETTE[i % CHART_PALETTE.length]}18`,
        line: {
          color: c.colour ?? CHART_PALETTE[i % CHART_PALETTE.length],
          width: 2,
        },
        marker: {
          size: 5,
          color: c.colour ?? CHART_PALETTE[i % CHART_PALETTE.length],
        },
        name: c.name,
        hovertemplate: "%{theta}: %{r:.2f}<extra>%{fullData.name}</extra>",
      })),
    [candidates, metrics, labels]
  );

  const layout = useMemo(
    () => ({
      ...PLOTLY_DARK_LAYOUT,
      height,
      showlegend: candidates.length > 1,
      legend: {
        font: { color: colours.textSecondary, size: 11 },
        x: 0.5,
        xanchor: "center" as const,
        y: -0.15,
        orientation: "h" as const,
      },
      polar: {
        bgcolor: colours.bgSurface,
        radialaxis: {
          visible: true,
          range: [0, 1],
          tickvals: [0.2, 0.4, 0.6, 0.8, 1.0],
          ticktext: ["0.2", "0.4", "0.6", "0.8", "1.0"],
          gridcolor: colours.border,
          linecolor: colours.border,
          tickfont: { color: colours.textTertiary, size: 9 },
        },
        angularaxis: {
          gridcolor: colours.border,
          linecolor: colours.border,
          tickfont: { color: colours.textSecondary, size: 11 },
        },
      },
      margin: { t: 24, r: 48, b: candidates.length > 1 ? 60 : 24, l: 48 },
    }),
    [height, candidates.length]
  );

  return (
    <Plot
      data={traces}
      layout={layout}
      config={{ displayModeBar: false, responsive: true }}
      useResizeHandler
      style={{ width: "100%", height }}
    />
  );
}
