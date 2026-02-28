"use client";

import type { MetricStatus } from "@/lib/types";

const FILL_COLOURS: Record<string, string> = {
  GREEN: "var(--green)",
  AMBER: "var(--amber)",
  RED: "var(--red)",
  GREY: "var(--grey)",
};

interface MetricBarProps {
  value: number;
  status: MetricStatus;
  compact?: boolean;
  showValue?: boolean;
}

export default function MetricBar({
  value,
  status,
  compact = false,
  showValue = false,
}: MetricBarProps) {
  const clamped = Math.max(0, Math.min(1, value));
  const height = compact ? 4 : 8;

  return (
    <div className={`flex items-center gap-2 ${compact ? "min-w-[48px]" : "min-w-[80px]"}`}>
      <div
        className="flex-1 rounded-full overflow-hidden"
        style={{ height, background: "var(--bg-inset)" }}
      >
        <div
          className="h-full rounded-full transition-all duration-300"
          style={{
            width: `${clamped * 100}%`,
            background: FILL_COLOURS[status],
            opacity: 0.85,
          }}
        />
      </div>
      {showValue && (
        <span
          className="text-xs font-mono tabular-nums shrink-0"
          style={{ color: "var(--text-secondary)", minWidth: 32, textAlign: "right" }}
        >
          {clamped.toFixed(2)}
        </span>
      )}
    </div>
  );
}
