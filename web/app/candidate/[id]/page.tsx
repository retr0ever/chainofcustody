"use client";

import { use } from "react";
import Link from "next/link";
import { useCandidates } from "@/lib/candidates";
import type { MetricStatus } from "@/lib/types";
import { METRIC_ORDER, formatMetricName, formatScore } from "@/lib/format";
import RadarChart from "../../components/RadarChart";
import MetricDetail from "../../components/MetricDetail";
import SequenceViewer from "../../components/SequenceViewer";
import StatusBadge from "../../components/StatusBadge";

export default function CandidateDetailPage({
  params,
}: {
  params: Promise<{ id: string }>;
}) {
  const { id } = use(params);
  const { getCandidate, loading } = useCandidates();
  const candidate = getCandidate(id);

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

  if (!candidate) {
    return (
      <div className="p-8">
        <div
          className="rounded-lg px-5 py-4 text-sm"
          style={{ background: "var(--red-bg)", color: "var(--red)", border: "1px solid var(--red)" }}
        >
          Candidate <code className="font-mono">{id}</code> not found.
        </div>
        <Link
          href="/"
          className="inline-block mt-4 text-sm"
          style={{ color: "var(--primary)" }}
        >
          &larr; Back to dashboard
        </Link>
      </div>
    );
  }

  const c = candidate;
  const radarScores: Record<string, number> = {};
  for (const m of METRIC_ORDER) {
    radarScores[m] = c.fitness.scores[m]?.value ?? 0;
  }

  const overallStatus: MetricStatus =
    c.fitness.overall >= 0.7 ? "GREEN" : c.fitness.overall >= 0.4 ? "AMBER" : "RED";

  return (
    <div className="p-8 max-w-7xl mx-auto">
      {/* Breadcrumb */}
      <div className="flex items-center gap-2 mb-6 text-xs" style={{ color: "var(--text-tertiary)" }}>
        <Link href="/" className="transition-colors" style={{ color: "var(--text-secondary)" }}
          onMouseEnter={(e) => (e.currentTarget.style.color = "var(--primary)")}
          onMouseLeave={(e) => (e.currentTarget.style.color = "var(--text-secondary)")}
        >
          Dashboard
        </Link>
        <span>/</span>
        <span style={{ color: "var(--text-primary)" }}>{c.name}</span>
      </div>

      {/* Header */}
      <div className="flex flex-col sm:flex-row sm:items-center gap-4 mb-8">
        <div className="flex-1">
          <h1 className="text-2xl font-semibold tracking-tight" style={{ color: "var(--text-primary)" }}>
            {c.name}
          </h1>
          <p className="text-sm mt-1" style={{ color: "var(--text-secondary)" }}>
            {c.gene} &middot; {c.target_cell_type}
          </p>
        </div>

        <div className="flex items-center gap-4">
          {/* Overall score */}
          <div
            className="flex items-center gap-3 rounded-xl px-5 py-3 border"
            style={{ background: "var(--bg-surface)", borderColor: "var(--border)" }}
          >
            <span className="text-xs uppercase tracking-wider" style={{ color: "var(--text-secondary)" }}>
              Fitness
            </span>
            <span
              className="text-2xl font-bold font-mono tabular-nums"
              style={{ color: "var(--text-primary)" }}
            >
              {formatScore(c.fitness.overall)}
            </span>
            <StatusBadge status={overallStatus} size="md" showLabel />
          </div>
        </div>
      </div>

      {/* Status badges row */}
      <div className="flex flex-wrap gap-2 mb-8">
        {METRIC_ORDER.map((m) => (
          <div
            key={m}
            className="flex items-center gap-2 rounded-lg px-3 py-2 border"
            style={{ background: "var(--bg-surface)", borderColor: "var(--border)" }}
          >
            <StatusBadge status={c.fitness.scores[m]?.status as MetricStatus} size="sm" />
            <span className="text-xs" style={{ color: "var(--text-secondary)" }}>
              {formatMetricName(m)}
            </span>
            <span className="text-xs font-mono" style={{ color: "var(--text-primary)" }}>
              {formatScore(c.fitness.scores[m]?.value)}
            </span>
          </div>
        ))}
      </div>

      {/* Two-column layout */}
      <div className="grid grid-cols-1 lg:grid-cols-[1fr_380px] gap-8 items-start">
        {/* Left: radar + metric details */}
        <div className="flex flex-col gap-4">
          {/* Radar chart */}
          <div
            className="rounded-xl border p-5"
            style={{ background: "var(--bg-surface)", borderColor: "var(--border)" }}
          >
            <h2
              className="text-sm font-semibold mb-2 uppercase tracking-wider"
              style={{ color: "var(--text-secondary)" }}
            >
              Metric Overview
            </h2>
            <RadarChart
              candidates={[{ id: c.id, name: c.name, scores: radarScores }]}
              height={340}
            />
          </div>

          {/* Metric detail sections */}
          <div className="flex flex-col gap-3">
            {METRIC_ORDER.map((m) => {
              const metricSuggestions = c.fitness.suggestions.filter((s) => s.metric === m);
              return (
                <MetricDetail
                  key={m}
                  metricKey={m}
                  score={c.fitness.scores[m]}
                  report={c.report}
                  suggestions={metricSuggestions.length > 0 ? metricSuggestions : undefined}
                />
              );
            })}
          </div>
        </div>

        {/* Right: sequence info + viewer */}
        <div className="flex flex-col gap-4 lg:sticky lg:top-8">
          {/* Sequence info */}
          <div
            className="rounded-xl border p-5"
            style={{ background: "var(--bg-surface)", borderColor: "var(--border)" }}
          >
            <h2
              className="text-sm font-semibold mb-4 uppercase tracking-wider"
              style={{ color: "var(--text-secondary)" }}
            >
              Sequence Info
            </h2>
            <div className="flex flex-col gap-2">
              {[
                ["Total length", `${c.report.sequence_info.total_length} nt`],
                ["5\u2032UTR", `${c.report.sequence_info.utr5_length} nt`],
                ["CDS", `${c.report.sequence_info.cds_length} nt (${c.report.sequence_info.num_codons} codons)`],
                ["3\u2032UTR", `${c.report.sequence_info.utr3_length} nt`],
              ].map(([label, value]) => (
                <div key={label} className="flex items-center justify-between">
                  <span className="text-xs" style={{ color: "var(--text-secondary)" }}>{label}</span>
                  <span className="text-xs font-mono tabular-nums" style={{ color: "var(--text-primary)" }}>{value}</span>
                </div>
              ))}
            </div>
          </div>

          {/* Sequence viewer */}
          <div
            className="rounded-xl border p-5"
            style={{ background: "var(--bg-surface)", borderColor: "var(--border)" }}
          >
            <h2
              className="text-sm font-semibold mb-4 uppercase tracking-wider"
              style={{ color: "var(--text-secondary)" }}
            >
              Construct
            </h2>
            <SequenceViewer
              utr5={c.sequence.utr5}
              cds={c.sequence.cds}
              utr3={c.sequence.utr3}
            />
          </div>
        </div>
      </div>
    </div>
  );
}
