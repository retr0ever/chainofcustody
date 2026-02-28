/** Formatting helpers for metric display. */

const METRIC_LABELS: Record<string, string> = {
  utr5_accessibility: "5\u2032UTR Accessibility",
  manufacturability: "Manufacturability",
  stability: "Stability",
  specificity: "Specificity",
};

const METRIC_SHORT: Record<string, string> = {
  utr5_accessibility: "UTR5",
  manufacturability: "Mfg",
  stability: "Stab",
  specificity: "Spec",
};

export function formatMetricName(key: string): string {
  return METRIC_LABELS[key] ?? key;
}

export function formatMetricShort(key: string): string {
  return METRIC_SHORT[key] ?? key;
}

export function formatScore(value: number | null | undefined, decimals = 2): string {
  if (value == null) return "\u2014";
  return value.toFixed(decimals);
}

export function formatPercentage(value: number | null | undefined): string {
  if (value == null) return "\u2014";
  return `${(value * 100).toFixed(0)}%`;
}

export function formatTimestamp(iso: string): string {
  const d = new Date(iso);
  return d.toLocaleDateString("en-GB", {
    day: "numeric",
    month: "short",
    year: "numeric",
  });
}

/** Metric keys in display order. */
export const METRIC_ORDER = [
  "utr5_accessibility",
  "manufacturability",
  "stability",
  "specificity",
] as const;
