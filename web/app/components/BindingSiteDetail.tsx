"use client";

import { useState } from "react";
import type { SpongeSite } from "@/lib/sponge";

interface BindingSiteDetailProps {
  sites: SpongeSite[];
}

const SITE_COLOURS = [
  "#5B8FA8", "#E8766A", "#34D399", "#FBBF24",
  "#A78BFA", "#F472B6", "#67E8F9", "#FCA5A1",
];

/**
 * Build the alignment between miRNA (5'→3') and the binding site (3'→5').
 * Shows Watson-Crick pairs as |, G-U wobble as :, mismatches as space.
 */
function buildAlignment(site: SpongeSite) {
  const mirna = site.mirnaSeq;
  const siteSeq = site.siteSeq;

  // The binding site is the reverse complement of the miRNA,
  // read 3'→5', so we reverse it for alignment against the 5'→3' miRNA.
  const siteRev = siteSeq.split("").reverse().join("");

  const WC = new Set(["AU", "UA", "GC", "CG"]);
  const WOBBLE = new Set(["GU", "UG"]);

  // miRNA positions 2-8 (1-indexed) are the seed
  // Bulge is at positions 9-12 from 5' end of miRNA
  const lines = { mirna: "", bonds: "", site: "" };

  const len = Math.max(mirna.length, siteRev.length);
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

  return lines;
}

function DuplexView({ site, index }: { site: SpongeSite; index: number }) {
  const alignment = buildAlignment(site);
  const colour = SITE_COLOURS[index % SITE_COLOURS.length];
  const seedStart = 1; // 0-indexed position 1 = miRNA pos 2
  const seedEnd = 7;   // through position 7 = miRNA pos 8
  const bulgeStart = 8;
  const bulgeEnd = 11;

  return (
    <div className="flex flex-col gap-2">
      {/* Duplex alignment */}
      <div
        className="rounded-lg p-3 font-mono text-xs leading-relaxed overflow-x-auto"
        style={{ background: "var(--bg-inset)" }}
      >
        {/* miRNA strand */}
        <div className="flex items-center gap-2 whitespace-nowrap">
          <span className="text-[10px] w-8 shrink-0 text-right" style={{ color: "var(--text-tertiary)" }}>
            5&apos;
          </span>
          <span>
            {alignment.mirna.split("").map((ch, i) => {
              let bg = "transparent";
              if (i >= seedStart && i <= seedEnd) bg = "rgba(91, 143, 168, 0.15)";
              else if (i >= bulgeStart && i <= bulgeEnd) bg = "rgba(248, 113, 113, 0.12)";
              return (
                <span key={i} style={{ background: bg, color: "var(--text-primary)" }}>{ch}</span>
              );
            })}
          </span>
          <span className="text-[10px] shrink-0" style={{ color: "var(--text-tertiary)" }}>
            3&apos; miRNA
          </span>
        </div>

        {/* Bond line */}
        <div className="flex items-center gap-2 whitespace-nowrap">
          <span className="w-8 shrink-0" />
          <span style={{ color: "var(--text-tertiary)" }}>
            {alignment.bonds.split("").map((ch, i) => {
              let bg = "transparent";
              if (i >= seedStart && i <= seedEnd) bg = "rgba(91, 143, 168, 0.15)";
              else if (i >= bulgeStart && i <= bulgeEnd) bg = "rgba(248, 113, 113, 0.12)";
              return <span key={i} style={{ background: bg }}>{ch}</span>;
            })}
          </span>
        </div>

        {/* Binding site strand */}
        <div className="flex items-center gap-2 whitespace-nowrap">
          <span className="text-[10px] w-8 shrink-0 text-right" style={{ color: "var(--text-tertiary)" }}>
            3&apos;
          </span>
          <span>
            {alignment.site.split("").map((ch, i) => {
              let bg = "transparent";
              if (i >= seedStart && i <= seedEnd) bg = "rgba(91, 143, 168, 0.15)";
              else if (i >= bulgeStart && i <= bulgeEnd) bg = "rgba(248, 113, 113, 0.12)";
              return (
                <span key={i} style={{ background: bg, color: colour }}>{ch}</span>
              );
            })}
          </span>
          <span className="text-[10px] shrink-0" style={{ color: "var(--text-tertiary)" }}>
            5&apos; site
          </span>
        </div>
      </div>

      {/* Region annotations */}
      <div className="flex flex-wrap items-center gap-2 sm:gap-3 text-[10px]">
        <div className="flex items-center gap-1">
          <span className="inline-block w-3 h-2 rounded-sm" style={{ background: "rgba(91, 143, 168, 0.3)" }} />
          <span style={{ color: "var(--text-tertiary)" }}>Seed (pos 2-8)</span>
        </div>
        <div className="flex items-center gap-1">
          <span className="inline-block w-3 h-2 rounded-sm" style={{ background: "rgba(248, 113, 113, 0.25)" }} />
          <span style={{ color: "var(--text-tertiary)" }}>Bulge mismatch</span>
        </div>
        <div className="flex items-center gap-1">
          <span style={{ color: "var(--text-tertiary)" }}>| W-C pair</span>
        </div>
        <div className="flex items-center gap-1">
          <span style={{ color: "var(--text-tertiary)" }}>: wobble</span>
        </div>
      </div>
    </div>
  );
}

export default function BindingSiteDetail({ sites }: BindingSiteDetailProps) {
  const [selected, setSelected] = useState(0);
  const site = sites[selected];

  if (!site) return null;

  return (
    <div className="flex flex-col gap-4">
      {/* Tabs if multiple sites */}
      {sites.length > 1 && (
        <div className="flex flex-wrap gap-1.5">
          {sites.map((s, i) => (
            <button
              key={s.mirnaId}
              onClick={() => setSelected(i)}
              className="text-xs px-3 py-1.5 rounded-md transition-colors cursor-pointer font-mono"
              style={{
                background: selected === i ? "var(--primary-bg)" : "var(--bg-raised)",
                color: selected === i ? SITE_COLOURS[i % SITE_COLOURS.length] : "var(--text-secondary)",
                border: `1px solid ${selected === i ? SITE_COLOURS[i % SITE_COLOURS.length] : "var(--border)"}`,
              }}
            >
              {s.mirnaId}
            </button>
          ))}
        </div>
      )}

      {/* Selected site info */}
      <div className="flex flex-wrap items-center gap-2 sm:gap-3 text-xs" style={{ color: "var(--text-secondary)" }}>
        <span className="font-semibold font-mono" style={{ color: SITE_COLOURS[selected % SITE_COLOURS.length] }}>
          {site.mirnaId}
        </span>
        <span>{site.mirnaSeq.length} nt mature miRNA</span>
        <span>{site.siteSeq.length} nt binding site</span>
      </div>

      <DuplexView site={site} index={selected} />

      {/* Site sequence breakdown */}
      <div
        className="rounded-lg overflow-hidden border text-xs"
        style={{ borderColor: "var(--border)" }}
      >
        <table className="w-full">
          <thead style={{ background: "var(--bg-inset)" }}>
            <tr style={{ borderBottom: "1px solid var(--border)" }}>
              <th className="px-3 py-1.5 text-left font-medium" style={{ color: "var(--text-secondary)" }}>Domain</th>
              <th className="px-3 py-1.5 text-left font-medium" style={{ color: "var(--text-secondary)" }}>Sequence</th>
              <th className="px-3 py-1.5 text-right font-medium" style={{ color: "var(--text-secondary)" }}>Length</th>
            </tr>
          </thead>
          <tbody>
            <tr style={{ borderBottom: "1px solid var(--border)" }}>
              <td className="px-3 py-1.5" style={{ color: "var(--text-primary)" }}>3&apos; supplementary</td>
              <td className="px-3 py-1.5 font-mono" style={{ color: "var(--text-secondary)" }}>{site.threePrimeMatch}</td>
              <td className="px-3 py-1.5 text-right font-mono" style={{ color: "var(--text-tertiary)" }}>{site.threePrimeMatch.length}</td>
            </tr>
            <tr style={{ borderBottom: "1px solid var(--border)" }}>
              <td className="px-3 py-1.5" style={{ color: "var(--red)" }}>Bulge (mismatch)</td>
              <td className="px-3 py-1.5 font-mono" style={{ color: "var(--red)" }}>{site.bulgeMismatch}</td>
              <td className="px-3 py-1.5 text-right font-mono" style={{ color: "var(--text-tertiary)" }}>{site.bulgeMismatch.length}</td>
            </tr>
            <tr style={{ borderBottom: "1px solid var(--border)" }}>
              <td className="px-3 py-1.5" style={{ color: "var(--primary)" }}>Seed match (pos 2-8)</td>
              <td className="px-3 py-1.5 font-mono" style={{ color: "var(--primary)" }}>{site.seedMatch}</td>
              <td className="px-3 py-1.5 text-right font-mono" style={{ color: "var(--text-tertiary)" }}>{site.seedMatch.length}</td>
            </tr>
          </tbody>
        </table>
      </div>
    </div>
  );
}
