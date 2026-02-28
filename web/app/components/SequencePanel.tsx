"use client";

import { useState, useCallback } from "react";
import type { SpongeResult } from "@/lib/sponge";

interface SequencePanelProps {
  sponge: SpongeResult;
}

function CopyButton({ text, label }: { text: string; label?: string }) {
  const [copied, setCopied] = useState(false);

  const copy = useCallback(() => {
    navigator.clipboard.writeText(text).then(() => {
      setCopied(true);
      setTimeout(() => setCopied(false), 1500);
    });
  }, [text]);

  return (
    <button
      onClick={copy}
      className="inline-flex items-center gap-1 text-xs px-2 py-1 rounded-md transition-colors cursor-pointer shrink-0"
      style={{
        background: copied ? "var(--green-bg)" : "var(--bg-raised)",
        color: copied ? "var(--green)" : "var(--text-secondary)",
        border: `1px solid ${copied ? "var(--green)" : "var(--border)"}`,
      }}
    >
      {copied ? (
        <svg viewBox="0 0 16 16" fill="currentColor" className="w-3 h-3">
          <path d="M13.78 4.22a.75.75 0 010 1.06l-7.25 7.25a.75.75 0 01-1.06 0L2.22 9.28a.75.75 0 011.06-1.06L6 10.94l6.72-6.72a.75.75 0 011.06 0z" />
        </svg>
      ) : (
        <svg viewBox="0 0 16 16" fill="currentColor" className="w-3 h-3">
          <path d="M0 6.75C0 5.784.784 5 1.75 5h1.5a.75.75 0 010 1.5h-1.5a.25.25 0 00-.25.25v7.5c0 .138.112.25.25.25h7.5a.25.25 0 00.25-.25v-1.5a.75.75 0 011.5 0v1.5A1.75 1.75 0 019.25 16h-7.5A1.75 1.75 0 010 14.25v-7.5z" />
          <path d="M5 1.75C5 .784 5.784 0 6.75 0h7.5C15.216 0 16 .784 16 1.75v7.5A1.75 1.75 0 0114.25 11h-7.5A1.75 1.75 0 015 9.25v-7.5zm1.75-.25a.25.25 0 00-.25.25v7.5c0 .138.112.25.25.25h7.5a.25.25 0 00.25-.25v-7.5a.25.25 0 00-.25-.25h-7.5z" />
        </svg>
      )}
      {label ?? "Copy"}
    </button>
  );
}

const NUC_COLOURS: Record<string, string> = {
  A: "var(--green)",
  U: "var(--red)",
  G: "var(--primary)",
  C: "var(--amber)",
};

export default function SequencePanel({ sponge }: SequencePanelProps) {
  const [view, setView] = useState<"full" | "cassette">("full");
  const seq = view === "full" ? sponge.fullUtr3 : sponge.cassette;

  return (
    <div className="flex flex-col gap-3">
      {/* Header with stats and copy */}
      <div className="flex flex-wrap items-center justify-between gap-2">
        <div className="flex items-center gap-2">
          <button
            onClick={() => setView("full")}
            className="text-xs px-3 py-1.5 rounded-md transition-colors cursor-pointer"
            style={{
              background: view === "full" ? "var(--primary-bg)" : "var(--bg-raised)",
              color: view === "full" ? "var(--primary)" : "var(--text-secondary)",
              border: `1px solid ${view === "full" ? "var(--primary)" : "var(--border)"}`,
            }}
          >
            Full 3&apos;UTR
          </button>
          <button
            onClick={() => setView("cassette")}
            className="text-xs px-3 py-1.5 rounded-md transition-colors cursor-pointer"
            style={{
              background: view === "cassette" ? "var(--primary-bg)" : "var(--bg-raised)",
              color: view === "cassette" ? "var(--primary)" : "var(--text-secondary)",
              border: `1px solid ${view === "cassette" ? "var(--primary)" : "var(--border)"}`,
            }}
          >
            Cassette only
          </button>
        </div>
        <CopyButton text={seq} label="Copy" />
      </div>

      {/* Stats */}
      <div className="flex flex-wrap items-center gap-2 sm:gap-3 text-[11px] sm:text-[10px]" style={{ color: "var(--text-tertiary)" }}>
        <span>{seq.length} nt</span>
        <span>{sponge.sites.length} binding site{sponge.sites.length !== 1 ? " types" : ""}</span>
        <span>{sponge.numSites} repeats</span>
      </div>

      {/* Sequence display */}
      <div
        className="rounded-lg p-3 overflow-x-auto max-h-48 overflow-y-auto font-mono text-xs leading-relaxed break-all"
        style={{ background: "var(--bg-inset)" }}
      >
        {seq.toUpperCase().split("").map((nt, i) => (
          <span
            key={i}
            style={{ color: NUC_COLOURS[nt] ?? "var(--text-secondary)" }}
          >
            {nt}
          </span>
        ))}
      </div>
    </div>
  );
}
