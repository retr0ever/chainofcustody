"use client";

import { useState, useCallback, useMemo } from "react";

interface SequenceViewerProps {
  utr5: string;
  cds: string;
  utr3: string;
  chunkSize?: number;
  showPositions?: boolean;
}

const KOZAK = "GCCACC";

export default function SequenceViewer({
  utr5,
  cds,
  utr3,
  chunkSize = 60,
  showPositions = true,
}: SequenceViewerProps) {
  const [copied, setCopied] = useState(false);
  const full = utr5 + cds + utr3;

  const handleCopy = useCallback(async () => {
    await navigator.clipboard.writeText(full);
    setCopied(true);
    setTimeout(() => setCopied(false), 2000);
  }, [full]);

  const regions = useMemo(() => {
    const spans: Array<{ char: string; region: "utr5" | "cds" | "utr3" }> = [];
    for (let i = 0; i < utr5.length; i++) spans.push({ char: utr5[i], region: "utr5" });
    for (let i = 0; i < cds.length; i++) spans.push({ char: cds[i], region: "cds" });
    for (let i = 0; i < utr3.length; i++) spans.push({ char: utr3[i], region: "utr3" });
    return spans;
  }, [utr5, cds, utr3]);

  const lines = useMemo(() => {
    const result: Array<{ pos: number; spans: typeof regions }> = [];
    for (let i = 0; i < regions.length; i += chunkSize) {
      result.push({
        pos: i + 1,
        spans: regions.slice(i, i + chunkSize),
      });
    }
    return result;
  }, [regions, chunkSize]);

  const colourMap = {
    utr5: "var(--seq-utr5)",
    cds: "var(--seq-cds)",
    utr3: "var(--seq-utr3)",
  };

  return (
    <div className="flex flex-col gap-3">
      {/* Legend + copy button */}
      <div className="flex items-center justify-between">
        <div className="flex items-center gap-4">
          {(["utr5", "cds", "utr3"] as const).map((r) => (
            <div key={r} className="flex items-center gap-1.5">
              <span
                className="w-2.5 h-2.5 rounded-sm"
                style={{ background: colourMap[r] }}
              />
              <span className="text-xs" style={{ color: "var(--text-secondary)" }}>
                {r === "utr5" ? "5\u2032UTR" : r === "cds" ? "CDS" : "3\u2032UTR"}
              </span>
            </div>
          ))}
        </div>
        <button
          onClick={handleCopy}
          className="text-xs px-2.5 py-1 rounded-md border transition-colors cursor-pointer"
          style={{
            borderColor: "var(--border)",
            color: copied ? "var(--green)" : "var(--text-secondary)",
            background: copied ? "var(--green-bg)" : "transparent",
          }}
        >
          {copied ? "Copied" : "Copy sequence"}
        </button>
      </div>

      {/* Sequence display */}
      <div
        className="rounded-lg p-4 overflow-x-auto"
        style={{ background: "var(--bg-inset)", border: "1px solid var(--border)" }}
      >
        <pre className="font-mono text-xs leading-6">
          {lines.map((line) => (
            <div key={line.pos} className="flex">
              {showPositions && (
                <span
                  className="select-none tabular-nums mr-4 shrink-0"
                  style={{ color: "var(--text-tertiary)", minWidth: 48, textAlign: "right" }}
                >
                  {line.pos}
                </span>
              )}
              <span className="break-all">
                {line.spans.map((s, j) => (
                  <span key={j} style={{ color: colourMap[s.region] }}>
                    {s.char}
                  </span>
                ))}
              </span>
            </div>
          ))}
        </pre>
      </div>

      {/* Info line */}
      <div className="flex gap-4 text-xs" style={{ color: "var(--text-tertiary)" }}>
        <span>{full.length} nt total</span>
        <span>{utr5.length} nt 5&prime;UTR</span>
        <span>{cds.length} nt CDS</span>
        <span>{utr3.length} nt 3&prime;UTR</span>
      </div>
    </div>
  );
}
