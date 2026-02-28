"use client";

import type { SpongeResult, SpongeRegion } from "@/lib/sponge";

interface SpongeArchitectureProps {
  sponge: SpongeResult;
}

/** Assign a stable colour to each miRNA by index */
const SITE_COLOURS = [
  "#5B8FA8", "#E8766A", "#34D399", "#FBBF24",
  "#A78BFA", "#F472B6", "#67E8F9", "#FCA5A1",
];

const REGION_META: Record<string, { label: string; colour: string }> = {
  stop: { label: "Stop", colour: "#F87171" },
  lead_in: { label: "Lead-in", colour: "#556378" },
  spacer: { label: "Spacer", colour: "#2A3A50" },
  lead_out: { label: "Lead-out", colour: "#556378" },
  poly_a: { label: "Poly-A signal", colour: "#8896AB" },
};

export default function SpongeArchitecture({ sponge }: SpongeArchitectureProps) {
  const totalLen = sponge.fullUtr3.length;

  // Unique miRNAs for legend
  const uniqueMirnas = sponge.sites.map((s, i) => ({
    id: s.mirnaId,
    colour: SITE_COLOURS[i % SITE_COLOURS.length],
  }));

  return (
    <div className="flex flex-col gap-4">
      {/* Architecture bar */}
      <div className="flex flex-col gap-1.5">
        <div className="flex items-center justify-between text-[10px] font-mono" style={{ color: "var(--text-tertiary)" }}>
          <span>5&apos;</span>
          <span>{totalLen} nt</span>
          <span>3&apos;</span>
        </div>

        {/* Segmented bar */}
        <div
          className="flex rounded-md overflow-hidden"
          style={{ height: 28 }}
        >
          {sponge.regions.map((region, i) => {
            const widthPct = ((region.end - region.start) / totalLen) * 100;
            const colour = region.type === "site"
              ? SITE_COLOURS[(region.mirnaIndex ?? 0) % SITE_COLOURS.length]
              : REGION_META[region.type]?.colour ?? "#2A3A50";

            return (
              <div
                key={i}
                className="relative group"
                style={{
                  width: `${widthPct}%`,
                  minWidth: region.type === "spacer" ? 2 : 1,
                  background: colour,
                  opacity: region.type === "spacer" ? 0.4 : 0.85,
                }}
                title={tooltipFor(region)}
              />
            );
          })}
        </div>

        {/* Position scale */}
        <div className="relative" style={{ height: 12 }}>
          {[0, 0.25, 0.5, 0.75, 1].map((frac) => (
            <span
              key={frac}
              className="absolute text-[8px] font-mono"
              style={{
                left: `${frac * 100}%`,
                transform: "translateX(-50%)",
                color: "var(--text-tertiary)",
              }}
            >
              {Math.round(frac * totalLen)}
            </span>
          ))}
        </div>
      </div>

      {/* Region breakdown */}
      <div
        className="rounded-lg overflow-x-auto border text-xs"
        style={{ borderColor: "var(--border)" }}
      >
        <table className="w-full">
          <thead style={{ background: "var(--bg-inset)" }}>
            <tr style={{ borderBottom: "1px solid var(--border)" }}>
              <th className="px-3 py-1.5 text-left font-medium" style={{ color: "var(--text-secondary)" }}>Region</th>
              <th className="px-3 py-1.5 text-right font-medium" style={{ color: "var(--text-secondary)" }}>Length</th>
              <th className="px-3 py-1.5 text-right font-medium" style={{ color: "var(--text-secondary)" }}>Position</th>
            </tr>
          </thead>
          <tbody>
            <RegionRow label="Stop codon" colour={REGION_META.stop.colour} length={3} start={0} />
            <RegionRow label="Lead-in" colour={REGION_META.lead_in.colour}
              length={sponge.regions.find((r) => r.type === "lead_in")?.seq.length ?? 0}
              start={sponge.regions.find((r) => r.type === "lead_in")?.start ?? 0} />

            {/* Binding sites grouped by miRNA */}
            {uniqueMirnas.map((m) => {
              const siteRegions = sponge.regions.filter((r) => r.type === "site" && r.mirnaId === m.id);
              const siteLen = siteRegions[0]?.seq.length ?? 0;
              return (
                <RegionRow key={m.id}
                  label={`${m.id} sites`}
                  colour={m.colour}
                  length={siteLen}
                  count={siteRegions.length}
                  start={siteRegions[0]?.start ?? 0}
                />
              );
            })}

            <RegionRow label="Spacers" colour={REGION_META.spacer.colour}
              length={4} count={sponge.numSites - 1}
              start={sponge.regions.find((r) => r.type === "spacer")?.start ?? 0} />
            <RegionRow label="Lead-out" colour={REGION_META.lead_out.colour}
              length={sponge.regions.find((r) => r.type === "lead_out")?.seq.length ?? 0}
              start={sponge.regions.find((r) => r.type === "lead_out")?.start ?? 0} />
            <RegionRow label="Poly-A signal" colour={REGION_META.poly_a.colour}
              length={sponge.regions.find((r) => r.type === "poly_a")?.seq.length ?? 0}
              start={sponge.regions.find((r) => r.type === "poly_a")?.start ?? 0} />
          </tbody>
        </table>
      </div>

      {/* Legend */}
      <div className="flex flex-wrap items-center gap-3">
        {uniqueMirnas.map((m) => (
          <div key={m.id} className="flex items-center gap-1.5">
            <span className="inline-block w-2.5 h-2.5 rounded-sm" style={{ background: m.colour }} />
            <span className="text-[10px] font-mono" style={{ color: "var(--text-secondary)" }}>{m.id}</span>
          </div>
        ))}
        <div className="flex items-center gap-1.5">
          <span className="inline-block w-2.5 h-2.5 rounded-sm" style={{ background: "#2A3A50" }} />
          <span className="text-[10px] font-mono" style={{ color: "var(--text-secondary)" }}>Spacer</span>
        </div>
        <div className="flex items-center gap-1.5">
          <span className="inline-block w-2.5 h-2.5 rounded-sm" style={{ background: "#8896AB" }} />
          <span className="text-[10px] font-mono" style={{ color: "var(--text-secondary)" }}>Poly-A</span>
        </div>
      </div>
    </div>
  );
}

function RegionRow({ label, colour, length, start, count }: {
  label: string; colour: string; length: number; start: number; count?: number;
}) {
  return (
    <tr
      className="transition-colors"
      style={{ borderBottom: "1px solid var(--border)" }}
      onMouseEnter={(e) => (e.currentTarget.style.background = "var(--bg-hover)")}
      onMouseLeave={(e) => (e.currentTarget.style.background = "transparent")}
    >
      <td className="px-3 py-1.5 flex items-center gap-2">
        <span className="inline-block w-2 h-2 rounded-sm shrink-0" style={{ background: colour }} />
        <span style={{ color: "var(--text-primary)" }}>
          {label}
          {count && count > 1 && (
            <span style={{ color: "var(--text-tertiary)" }}> ({count}x)</span>
          )}
        </span>
      </td>
      <td className="px-3 py-1.5 text-right font-mono" style={{ color: "var(--text-secondary)" }}>
        {count && count > 1 ? `${count} x ${length}` : length} nt
      </td>
      <td className="px-3 py-1.5 text-right font-mono" style={{ color: "var(--text-tertiary)" }}>
        {start + 1}
      </td>
    </tr>
  );
}

function tooltipFor(region: SpongeRegion): string {
  const label = region.type === "site"
    ? `Binding site: ${region.mirnaId}`
    : REGION_META[region.type]?.label ?? region.type;
  return `${label} (${region.start + 1}..${region.end}, ${region.end - region.start} nt)`;
}
