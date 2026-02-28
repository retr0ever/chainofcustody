"use client";

import type { MetricStatus } from "@/lib/types";

const DOT_COLOURS: Record<string, string> = {
  GREEN: "var(--green)",
  AMBER: "var(--amber)",
  RED: "var(--red)",
  GREY: "var(--grey)",
};

const BG_COLOURS: Record<string, string> = {
  GREEN: "var(--green-bg)",
  AMBER: "var(--amber-bg)",
  RED: "var(--red-bg)",
  GREY: "var(--grey-bg)",
};

const STATUS_LABELS: Record<string, string> = {
  GREEN: "Pass",
  AMBER: "Warn",
  RED: "Fail",
  GREY: "N/A",
};

interface StatusBadgeProps {
  status: MetricStatus;
  size?: "sm" | "md";
  label?: string;
  showLabel?: boolean;
}

export default function StatusBadge({
  status,
  size = "sm",
  label,
  showLabel = false,
}: StatusBadgeProps) {
  const dotSize = size === "sm" ? 6 : 8;
  const text = label ?? (showLabel ? STATUS_LABELS[status] : undefined);

  return (
    <span
      className={`inline-flex items-center gap-1.5 rounded-full ${
        size === "sm" ? "text-xs px-2 py-0.5" : "text-sm px-2.5 py-1"
      }`}
      style={{ background: BG_COLOURS[status] }}
    >
      <span
        className="rounded-full shrink-0"
        style={{
          width: dotSize,
          height: dotSize,
          background: DOT_COLOURS[status],
          boxShadow: `0 0 6px ${DOT_COLOURS[status]}40`,
        }}
      />
      {text && (
        <span style={{ color: DOT_COLOURS[status] }}>{text}</span>
      )}
    </span>
  );
}
