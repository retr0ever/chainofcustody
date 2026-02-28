"use client";

import { useMemo, useState } from "react";
import type { SpongeResult, SpongeSite } from "@/lib/sponge";

interface StructurePlotsProps {
  sponge: SpongeResult;
  mirnaNames: string[];
}

const SITE_COLOURS = [
  "#5B8FA8", "#E8766A", "#34D399", "#FBBF24",
  "#A78BFA", "#F472B6", "#67E8F9", "#FCA5A1",
];

const SPACER_COLOUR = "var(--border-strong)";

function colourForMirna(index: number): string {
  return SITE_COLOURS[index % SITE_COLOURS.length];
}

// ── 1D Linear Binding Site Map ──────────────────────────────

function LinearSpongeMap({ sponge, mirnaNames }: { sponge: SpongeResult; mirnaNames: string[] }) {
  const totalLen = sponge.fullUtr3.length;
  const siteRegions = sponge.regions.filter((r) => r.type === "site");

  // SVG dimensions
  const width = 900;
  const height = 260;
  const barY = 140;
  const barH = 28;
  const padX = 40;
  const usable = width - padX * 2;

  const toX = (pos: number) => padX + (pos / totalLen) * usable;

  // Unique miRNAs for legend
  const uniqueMirnas = useMemo(() => {
    const seen = new Map<string, number>();
    for (const r of siteRegions) {
      if (r.mirnaId && !seen.has(r.mirnaId)) {
        seen.set(r.mirnaId, r.mirnaIndex ?? 0);
      }
    }
    return Array.from(seen.entries()).map(([id, idx]) => ({ id, colour: colourForMirna(idx) }));
  }, [siteRegions]);

  return (
    <div className="overflow-x-auto -mx-4 px-4 sm:mx-0 sm:px-0">
      <svg
        viewBox={`0 0 ${width} ${height}`}
        className="w-full"
        style={{ minWidth: 600, maxHeight: 280 }}
      >
        {/* Title */}
        <text
          x={width / 2} y={20}
          textAnchor="middle"
          fill="var(--text-primary)"
          fontSize={13}
          fontWeight={600}
        >
          mRNA Sponge 3&#x2032;UTR &middot; {totalLen} nt &middot; {siteRegions.length} binding sites
        </text>

        {/* 5' and 3' labels */}
        <text x={padX - 8} y={barY + barH / 2 + 4} textAnchor="end" fill="var(--text-tertiary)" fontSize={11} fontFamily="var(--font-mono)">5&#x2032;</text>
        <text x={width - padX + 8} y={barY + barH / 2 + 4} textAnchor="start" fill="var(--text-tertiary)" fontSize={11} fontFamily="var(--font-mono)">3&#x2032;</text>

        {/* Backbone line */}
        <line x1={padX} y1={barY + barH / 2} x2={width - padX} y2={barY + barH / 2} stroke="var(--border-strong)" strokeWidth={2} />

        {/* Region blocks */}
        {sponge.regions.map((region, i) => {
          const x = toX(region.start);
          const w = Math.max(1, toX(region.end) - toX(region.start));

          let fill: string;
          let opacity = 1;
          switch (region.type) {
            case "stop":
              fill = "var(--red)";
              opacity = 0.5;
              break;
            case "lead_in":
            case "lead_out":
              fill = "var(--text-tertiary)";
              opacity = 0.4;
              break;
            case "spacer":
              fill = SPACER_COLOUR;
              opacity = 0.3;
              break;
            case "poly_a":
              fill = "var(--text-secondary)";
              opacity = 0.3;
              break;
            case "site":
              fill = colourForMirna(region.mirnaIndex ?? 0);
              opacity = 0.85;
              break;
            default:
              fill = "var(--border)";
          }

          return (
            <rect
              key={i}
              x={x}
              y={barY}
              width={w}
              height={barH}
              rx={2}
              fill={fill}
              opacity={opacity}
            />
          );
        })}

        {/* Site labels with staggered connectors */}
        {siteRegions.map((region, i) => {
          const cx = toX(region.start) + (toX(region.end) - toX(region.start)) / 2;
          const yText = i % 2 === 0 ? 60 : 85;
          const colour = colourForMirna(region.mirnaIndex ?? 0);
          const label = region.mirnaId ?? `Site ${i + 1}`;

          return (
            <g key={`label-${i}`}>
              <line
                x1={cx} y1={barY - 2}
                x2={cx} y2={yText + 12}
                stroke="var(--border-strong)"
                strokeWidth={0.5}
                opacity={0.5}
              />
              <text
                x={cx} y={yText}
                textAnchor="middle"
                fill={colour}
                fontSize={8}
                fontFamily="var(--font-mono)"
                transform={`rotate(-30 ${cx} ${yText})`}
              >
                {label} #{i + 1}
              </text>
            </g>
          );
        })}

        {/* Position scale */}
        {[0, 0.25, 0.5, 0.75, 1].map((frac) => {
          const x = toX(frac * totalLen);
          return (
            <g key={frac}>
              <line x1={x} y1={barY + barH + 4} x2={x} y2={barY + barH + 10} stroke="var(--text-tertiary)" strokeWidth={0.5} />
              <text x={x} y={barY + barH + 22} textAnchor="middle" fill="var(--text-tertiary)" fontSize={9} fontFamily="var(--font-mono)">
                {Math.round(frac * totalLen)}
              </text>
            </g>
          );
        })}

        {/* Legend */}
        {uniqueMirnas.map((m, i) => {
          const lx = padX + i * 130;
          const ly = height - 16;
          return (
            <g key={m.id}>
              <rect x={lx} y={ly - 6} width={10} height={10} rx={2} fill={m.colour} opacity={0.85} />
              <text x={lx + 14} y={ly + 3} fill="var(--text-secondary)" fontSize={9} fontFamily="var(--font-mono)">
                {m.id}
              </text>
            </g>
          );
        })}
      </svg>
    </div>
  );
}


// ── Nussinov Algorithm (basic RNA secondary structure prediction) ─────

function nussinovFold(seq: string): string {
  const n = seq.length;
  if (n < 5) return ".".repeat(n);

  const canPair = (a: string, b: string): boolean => {
    const pair = a + b;
    return pair === "AU" || pair === "UA" || pair === "GC" || pair === "CG" || pair === "GU" || pair === "UG";
  };

  // DP table
  const dp: number[][] = Array.from({ length: n }, () => new Array(n).fill(0));
  const MIN_LOOP = 3;

  for (let span = MIN_LOOP + 1; span < n; span++) {
    for (let i = 0; i < n - span; i++) {
      const j = i + span;
      // Case 1: j unpaired
      dp[i][j] = dp[i][j - 1];
      // Case 2: j paired with some k
      for (let k = i; k < j - MIN_LOOP; k++) {
        if (canPair(seq[k], seq[j])) {
          const score = (k > i ? dp[i][k - 1] : 0) + 1 + (k + 1 <= j - 1 ? dp[k + 1][j - 1] : 0);
          dp[i][j] = Math.max(dp[i][j], score);
        }
      }
    }
  }

  // Traceback
  const pairs = new Array(n).fill(-1);
  const traceback = (i: number, j: number) => {
    if (i >= j) return;
    if (dp[i][j] === dp[i][j - 1]) {
      traceback(i, j - 1);
    } else {
      for (let k = i; k < j - MIN_LOOP; k++) {
        if (canPair(seq[k], seq[j])) {
          const score = (k > i ? dp[i][k - 1] : 0) + 1 + (k + 1 <= j - 1 ? dp[k + 1][j - 1] : 0);
          if (dp[i][j] === score) {
            pairs[k] = j;
            pairs[j] = k;
            if (k > i) traceback(i, k - 1);
            traceback(k + 1, j - 1);
            return;
          }
        }
      }
    }
  };
  traceback(0, n - 1);

  const result = new Array(n).fill(".");
  for (let i = 0; i < n; i++) {
    if (pairs[i] > i) {
      result[i] = "(";
      result[pairs[i]] = ")";
    }
  }
  return result.join("");
}


// ── Arc Diagram (secondary structure visualisation) ──────────

function ArcDiagram({ sponge, mirnaNames }: { sponge: SpongeResult; mirnaNames: string[] }) {
  const [computing, setComputing] = useState(false);
  const seq = sponge.cassette.toUpperCase();

  const { structure, pairs } = useMemo(() => {
    const struct = nussinovFold(seq);
    const p: [number, number][] = [];
    const stack: number[] = [];
    for (let i = 0; i < struct.length; i++) {
      if (struct[i] === "(") stack.push(i);
      else if (struct[i] === ")") {
        const j = stack.pop();
        if (j !== undefined) p.push([j, i]);
      }
    }
    return { structure: struct, pairs: p };
  }, [seq]);

  // Build per-nucleotide colour map from sponge regions
  const ntColours = useMemo(() => {
    const colours = new Array(sponge.fullUtr3.length).fill(SPACER_COLOUR);
    for (const region of sponge.regions) {
      if (region.type === "site") {
        const c = colourForMirna(region.mirnaIndex ?? 0);
        for (let i = region.start; i < region.end; i++) {
          colours[i] = c;
        }
      }
    }
    // Offset to cassette region: find first site region
    const cassetteStart = sponge.regions.find((r) => r.type === "site")?.start ?? 0;
    return colours.slice(cassetteStart, cassetteStart + seq.length);
  }, [sponge, seq.length]);

  const n = seq.length;
  const width = 900;
  const lineY = 280;
  const padX = 20;
  const usable = width - padX * 2;
  const maxArcH = 240;

  const toX = (pos: number) => padX + (pos / n) * usable;

  // Filter to show only significant arcs (span > 4)
  const significantPairs = pairs.filter(([i, j]) => j - i > 4);
  const maxSpan = Math.max(1, ...significantPairs.map(([i, j]) => j - i));

  const numPairs = pairs.length;
  const gcContent = ((seq.match(/[GC]/g)?.length ?? 0) / n * 100).toFixed(1);

  return (
    <div>
      <div className="flex flex-wrap items-center gap-3 mb-3 text-xs" style={{ color: "var(--text-secondary)" }}>
        <span>{n} nt cassette</span>
        <span>{numPairs} predicted base pairs</span>
        <span>{gcContent}% GC</span>
      </div>

      <div className="overflow-x-auto -mx-4 px-4 sm:mx-0 sm:px-0">
        <svg
          viewBox={`0 0 ${width} ${lineY + 30}`}
          className="w-full"
          style={{ minWidth: 600 }}
        >
          {/* Arc pairs */}
          {significantPairs.map(([i, j], idx) => {
            const x1 = toX(i);
            const x2 = toX(j);
            const cx = (x1 + x2) / 2;
            const span = j - i;
            const ry = (span / maxSpan) * maxArcH * 0.85;
            const rx = (x2 - x1) / 2;

            return (
              <ellipse
                key={idx}
                cx={cx}
                cy={lineY}
                rx={rx}
                ry={Math.min(ry, maxArcH)}
                fill="none"
                stroke={ntColours[i] || "var(--border-strong)"}
                strokeWidth={0.6}
                opacity={0.35}
              />
            );
          })}

          {/* Sequence baseline - coloured blocks per region */}
          {sponge.regions.filter(r => r.type === "site" || r.type === "spacer").map((region, i) => {
            const cassetteStart = sponge.regions.find((r) => r.type === "site")?.start ?? 0;
            const relStart = region.start - cassetteStart;
            const relEnd = region.end - cassetteStart;
            if (relStart < 0 || relEnd > n) return null;

            const x = toX(relStart);
            const w = Math.max(1, toX(relEnd) - toX(relStart));
            const fill = region.type === "site"
              ? colourForMirna(region.mirnaIndex ?? 0)
              : SPACER_COLOUR;

            return (
              <rect
                key={i}
                x={x}
                y={lineY - 3}
                width={w}
                height={6}
                rx={1}
                fill={fill}
                opacity={region.type === "site" ? 0.8 : 0.3}
              />
            );
          })}

          {/* Position scale */}
          {[0, 0.25, 0.5, 0.75, 1].map((frac) => {
            const x = toX(frac * n);
            return (
              <g key={frac}>
                <line x1={x} y1={lineY + 8} x2={x} y2={lineY + 14} stroke="var(--text-tertiary)" strokeWidth={0.5} />
                <text x={x} y={lineY + 26} textAnchor="middle" fill="var(--text-tertiary)" fontSize={9} fontFamily="var(--font-mono)">
                  {Math.round(frac * n)}
                </text>
              </g>
            );
          })}
        </svg>
      </div>
    </div>
  );
}


// ── Duplex Detail (miRNA-site pairing) ───────────────────────

const WC = new Set(["AU", "UA", "GC", "CG"]);
const WOBBLE = new Set(["GU", "UG"]);

function DuplexView({ site, index }: { site: SpongeSite; index: number }) {
  const mirna = site.mirnaSeq;
  const siteRev = site.siteSeq.split("").reverse().join("");
  const colour = colourForMirna(index);

  const len = Math.max(mirna.length, siteRev.length);
  const lines = { mirna: "", bonds: "", site: "" };
  for (let i = 0; i < len; i++) {
    const m = mirna[i] ?? " ";
    const s = siteRev[i] ?? " ";
    lines.mirna += m;
    lines.site += s;
    const pair = m + s;
    if (WC.has(pair)) lines.bonds += "|";
    else if (WOBBLE.has(pair)) lines.bonds += ":";
    else lines.bonds += " ";
  }

  return (
    <div
      className="rounded-lg p-3 font-mono text-xs leading-relaxed overflow-x-auto"
      style={{ background: "var(--bg-inset)" }}
    >
      <div className="flex items-center gap-2 whitespace-nowrap">
        <span className="text-[10px] w-8 shrink-0 text-right" style={{ color: "var(--text-tertiary)" }}>5&#x2032;</span>
        <span style={{ color: "var(--text-primary)" }}>{lines.mirna}</span>
        <span className="text-[10px] shrink-0" style={{ color: "var(--text-tertiary)" }}>3&#x2032; miRNA</span>
      </div>
      <div className="flex items-center gap-2 whitespace-nowrap">
        <span className="w-8 shrink-0" />
        <span style={{ color: "var(--text-tertiary)" }}>{lines.bonds}</span>
      </div>
      <div className="flex items-center gap-2 whitespace-nowrap">
        <span className="text-[10px] w-8 shrink-0 text-right" style={{ color: "var(--text-tertiary)" }}>3&#x2032;</span>
        <span style={{ color: colour }}>{lines.site}</span>
        <span className="text-[10px] shrink-0" style={{ color: "var(--text-tertiary)" }}>5&#x2032; site</span>
      </div>
    </div>
  );
}


// ── Main Export ──────────────────────────────────────────────

export default function StructurePlots({ sponge, mirnaNames }: StructurePlotsProps) {
  const [selectedSite, setSelectedSite] = useState(0);
  const site = sponge.sites[selectedSite];

  return (
    <div className="flex flex-col gap-4 sm:gap-6">
      {/* 1D Linear binding site map */}
      <section
        className="rounded-xl border p-4 sm:p-5"
        style={{ background: "var(--bg-surface)", borderColor: "var(--border)" }}
      >
        <h2 className="text-sm font-semibold mb-3 sm:mb-4" style={{ color: "var(--text-primary)" }}>
          mRNA sponge — linear binding site map
        </h2>
        <LinearSpongeMap sponge={sponge} mirnaNames={mirnaNames} />
      </section>

      {/* Arc diagram — predicted secondary structure */}
      <section
        className="rounded-xl border p-4 sm:p-5"
        style={{ background: "var(--bg-surface)", borderColor: "var(--border)" }}
      >
        <h2 className="text-sm font-semibold mb-3 sm:mb-4" style={{ color: "var(--text-primary)" }}>
          mRNA sponge — predicted secondary structure
        </h2>
        <ArcDiagram sponge={sponge} mirnaNames={mirnaNames} />
      </section>

      {/* Binding site duplex detail */}
      {site && (
        <section
          className="rounded-xl border p-4 sm:p-5"
          style={{ background: "var(--bg-surface)", borderColor: "var(--border)" }}
        >
          <h2 className="text-sm font-semibold mb-3 sm:mb-4" style={{ color: "var(--text-primary)" }}>
            Binding site duplex detail
          </h2>

          {/* Site tabs */}
          {sponge.sites.length > 1 && (
            <div className="flex flex-wrap gap-1.5 mb-4">
              {sponge.sites.map((s, i) => (
                <button
                  key={s.mirnaId}
                  onClick={() => setSelectedSite(i)}
                  className="text-xs px-3 py-1.5 rounded-md transition-colors cursor-pointer font-mono"
                  style={{
                    background: selectedSite === i ? "var(--primary-bg)" : "var(--bg-raised)",
                    color: selectedSite === i ? colourForMirna(i) : "var(--text-secondary)",
                    border: `1px solid ${selectedSite === i ? colourForMirna(i) : "var(--border)"}`,
                  }}
                >
                  {s.mirnaId}
                </button>
              ))}
            </div>
          )}

          <div className="flex flex-wrap items-center gap-2 sm:gap-3 text-xs mb-3" style={{ color: "var(--text-secondary)" }}>
            <span className="font-semibold font-mono" style={{ color: colourForMirna(selectedSite) }}>
              {site.mirnaId}
            </span>
            <span>{site.mirnaSeq.length} nt mature miRNA</span>
            <span>{site.siteSeq.length} nt binding site</span>
          </div>

          <DuplexView site={site} index={selectedSite} />

          {/* Annotation legend */}
          <div className="flex flex-wrap items-center gap-2 sm:gap-3 text-[10px] mt-3">
            <div className="flex items-center gap-1">
              <span style={{ color: "var(--text-tertiary)" }}>| W-C pair</span>
            </div>
            <div className="flex items-center gap-1">
              <span style={{ color: "var(--text-tertiary)" }}>: G-U wobble</span>
            </div>
            <div className="flex items-center gap-1">
              <span className="inline-block w-3 h-2 rounded-sm" style={{ background: colourForMirna(selectedSite) }} />
              <span style={{ color: "var(--text-tertiary)" }}>Binding site</span>
            </div>
          </div>
        </section>
      )}
    </div>
  );
}
