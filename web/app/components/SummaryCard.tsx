"use client";

interface SummaryCardProps {
  title: string;
  value: string | number;
  subtitle?: string;
  icon?: React.ReactNode;
}

export default function SummaryCard({ title, value, subtitle, icon }: SummaryCardProps) {
  return (
    <div
      className="rounded-xl p-5 border"
      style={{
        background: "var(--bg-surface)",
        borderColor: "var(--border)",
      }}
    >
      <div className="flex items-start justify-between">
        <div className="flex flex-col gap-1">
          <span
            className="text-xs font-medium uppercase tracking-wider"
            style={{ color: "var(--text-secondary)" }}
          >
            {title}
          </span>
          <span
            className="text-2xl font-semibold font-mono tabular-nums"
            style={{ color: "var(--text-primary)" }}
          >
            {value}
          </span>
          {subtitle && (
            <span className="text-xs" style={{ color: "var(--text-tertiary)" }}>
              {subtitle}
            </span>
          )}
        </div>
        {icon && (
          <div
            className="w-9 h-9 rounded-lg flex items-center justify-center shrink-0"
            style={{ background: "var(--primary-bg)", color: "var(--primary)" }}
          >
            {icon}
          </div>
        )}
      </div>
    </div>
  );
}
