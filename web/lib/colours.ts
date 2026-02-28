/** Colour tokens and Plotly config — theme-aware via CSS variables. */

/**
 * Static dark palette kept for backwards compat / non-DOM contexts.
 * Prefer reading CSS variables at runtime for theme awareness.
 */
export const colours = {
  bgBase: "#0C1017",
  bgSurface: "#131924",
  bgRaised: "#1A2233",
  bgHover: "#212D42",
  bgInset: "#0A0E14",

  textPrimary: "#E8ECF2",
  textSecondary: "#8896AB",
  textTertiary: "#556378",

  border: "#1E2A3A",
  borderStrong: "#2A3A50",

  primary: "#5B8FA8",
  primaryHover: "#4A7A94",
  accent: "#E8766A",

  green: "#34D399",
  greenBg: "#0D2A1F",
  amber: "#FBBF24",
  amberBg: "#2A2210",
  red: "#F87171",
  redBg: "#2A1414",
  grey: "#6B7280",

  seqUtr5: "#8CB4C8",
  seqCds: "#5B8FA8",
  seqUtr3: "#E8766A",
} as const;

/** Status colour mapping (works on both themes — green/amber/red adapt via CSS vars). */
export const statusColour: Record<string, string> = {
  GREEN: "var(--green)",
  AMBER: "var(--amber)",
  RED: "var(--red)",
  GREY: "var(--grey)",
};

export const statusBg: Record<string, string> = {
  GREEN: "var(--green-bg)",
  AMBER: "var(--amber-bg)",
  RED: "var(--red-bg)",
  GREY: "var(--grey-bg)",
};

/** Plotly trace palette for multi-series charts. */
export const CHART_PALETTE = [
  "#5B8FA8",
  "#E8766A",
  "#8CB4C8",
  "#34D399",
  "#FBBF24",
  "#A78BFA",
  "#F472B6",
  "#38BDF8",
  "#FB923C",
  "#6EE7B7",
];

/** Read a CSS variable value from the document root. */
function getCssVar(name: string, fallback: string): string {
  if (typeof window === "undefined") return fallback;
  return getComputedStyle(document.documentElement).getPropertyValue(name).trim() || fallback;
}

/** Build Plotly layout matching the current CSS theme. */
export function getPlotlyLayout(): Record<string, unknown> {
  const bgSurface = getCssVar("--bg-surface", colours.bgSurface);
  const textSec = getCssVar("--text-secondary", colours.textSecondary);
  const border = getCssVar("--border", colours.border);

  return {
    plot_bgcolor: bgSurface,
    paper_bgcolor: bgSurface,
    font: {
      color: textSec,
      family: "var(--font-geist-sans), system-ui, sans-serif",
    },
    xaxis: {
      gridcolor: border,
      zerolinecolor: border,
      tickfont: { color: textSec },
    },
    yaxis: {
      gridcolor: border,
      zerolinecolor: border,
      tickfont: { color: textSec },
    },
    legend: {
      font: { color: textSec },
    },
    margin: { t: 32, r: 16, b: 40, l: 48 },
  };
}

/** Legacy export — prefer getPlotlyLayout() for theme awareness. */
export const PLOTLY_DARK_LAYOUT: Record<string, unknown> = {
  plot_bgcolor: colours.bgSurface,
  paper_bgcolor: colours.bgSurface,
  font: {
    color: colours.textSecondary,
    family: "var(--font-geist-sans), system-ui, sans-serif",
  },
  xaxis: {
    gridcolor: colours.border,
    zerolinecolor: colours.border,
    tickfont: { color: colours.textSecondary },
  },
  yaxis: {
    gridcolor: colours.border,
    zerolinecolor: colours.border,
    tickfont: { color: colours.textSecondary },
  },
  legend: {
    font: { color: colours.textSecondary },
  },
  margin: { t: 32, r: 16, b: 40, l: 48 },
};
