"use client";

import { useEffect, useMemo, useState } from "react";
import type { UtrDesignResult, BindingSite } from "@/lib/utr-design";

interface StructurePlotsProps {
  design: UtrDesignResult;
  mirnaNames: string[];
}

const SITE_COLOURS = [
  "#5B8FA8", "#E8766A", "#34D399", "#FBBF24",
  "#A78BFA", "#F472B6", "#67E8F9", "#FCA5A1",
];

function colourForMirna(index: number): string {
  return SITE_COLOURS[index % SITE_COLOURS.length];
}

const STRUCTURE_API = process.env.NEXT_PUBLIC_STRUCTURE_API ?? "http://localhost:8000";

interface FoldResponse {
  plot_1d: string | null;
  plot_2d: string | null;
  dot_bracket: string | null;
  mfe: number | null;
}

// ── API-based structure plots ─────────────────────────────────

function useStructureApi(design: UtrDesignResult, mirnaNames: string[]) {
  const [data, setData] = useState<FoldResponse | null>(null);
  const [loading, setLoading] = useState(false);
  const [error, setError] = useState<string | null>(null);

  const sequence = design.fullUtr3;
  const namesKey = mirnaNames.join(",");

  useEffect(() => {
    if (!sequence) return;

    let cancelled = false;
    setLoading(true);
    setError(null);

    fetch(`${STRUCTURE_API}/api/fold`, {
      method: "POST",
      headers: { "Content-Type": "application/json" },
      body: JSON.stringify({ sequence, mirna_names: mirnaNames }),
    })
      .then((res) => {
        if (!res.ok) throw new Error(`API returned ${res.status}`);
        return res.json();
      })
      .then((json: FoldResponse) => {
        if (!cancelled) setData(json);
      })
      .catch((err) => {
        if (!cancelled) setError(err.message);
      })
      .finally(() => {
        if (!cancelled) setLoading(false);
      });

    return () => { cancelled = true; };
  }, [sequence, namesKey]);

  return { data, loading, error };
}


// ── Duplex Detail (miRNA-site pairing) ───────────────────────

const WC = new Set(["AU", "UA", "GC", "CG"]);
const WOBBLE = new Set(["GU", "UG"]);

function DuplexView({ site, index }: { site: BindingSite; index: number }) {
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


// ── SVG Display ──────────────────────────────────────────────

function SvgPanel({ svg, label }: { svg: string; label: string }) {
  return (
    <div className="overflow-x-auto -mx-4 px-4 sm:mx-0 sm:px-0">
      <div
        className="w-full [&_svg]:w-full [&_svg]:h-auto [&_svg]:max-w-full"
        dangerouslySetInnerHTML={{ __html: svg }}
      />
    </div>
  );
}

function LoadingSpinner({ label }: { label: string }) {
  return (
    <div className="flex items-center justify-center gap-3 py-12 text-sm" style={{ color: "var(--text-secondary)" }}>
      <svg className="animate-spin w-5 h-5" style={{ color: "var(--primary)" }} viewBox="0 0 24 24" fill="none">
        <circle className="opacity-25" cx="12" cy="12" r="10" stroke="currentColor" strokeWidth="4" />
        <path className="opacity-75" fill="currentColor" d="M4 12a8 8 0 018-8v4a4 4 0 00-4 4H4z" />
      </svg>
      {label}
    </div>
  );
}

function ApiError({ message }: { message: string }) {
  return (
    <div
      className="rounded-lg px-4 py-3 text-xs"
      style={{ background: "var(--bg-inset)", color: "var(--text-tertiary)" }}
    >
      <p className="font-medium mb-1" style={{ color: "var(--text-secondary)" }}>Structure API unavailable</p>
      <p>Start the API server to see structure plots:</p>
      <code className="block mt-1.5 px-2 py-1 rounded" style={{ background: "var(--bg-raised)" }}>
        uv run uvicorn dashboard.api:app --port 8000
      </code>
    </div>
  );
}


// ── Main Export ──────────────────────────────────────────────

export default function StructurePlots({ design, mirnaNames }: StructurePlotsProps) {
  const [selectedSite, setSelectedSite] = useState(0);
  const site = design.sites[selectedSite];
  const { data, loading, error } = useStructureApi(design, mirnaNames);

  return (
    <div className="flex flex-col gap-4 sm:gap-6">
      {/* 1D Linear binding site map (from API) */}
      <section
        className="rounded-xl border p-4 sm:p-5"
        style={{ background: "var(--bg-surface)", borderColor: "var(--border)" }}
      >
        <h2 className="text-sm font-semibold mb-3 sm:mb-4" style={{ color: "var(--text-primary)" }}>
          3&apos;UTR cassette map
        </h2>
        {loading && <LoadingSpinner label="Generating cassette map..." />}
        {error && <ApiError message={error} />}
        {data?.plot_1d && <SvgPanel svg={data.plot_1d} label="1D cassette map" />}
      </section>

      {/* 2D Secondary structure (from API) */}
      <section
        className="rounded-xl border p-4 sm:p-5"
        style={{ background: "var(--bg-surface)", borderColor: "var(--border)" }}
      >
        <h2 className="text-sm font-semibold mb-3 sm:mb-4" style={{ color: "var(--text-primary)" }}>
          Predicted secondary structure
        </h2>
        {loading && <LoadingSpinner label="Folding with ViennaRNA..." />}
        {error && <ApiError message={error} />}
        {data?.plot_2d && (
          <>
            <SvgPanel svg={data.plot_2d} label="2D structure" />
            {data.mfe !== null && (
              <div className="flex flex-wrap items-center gap-3 mt-3 text-xs" style={{ color: "var(--text-secondary)" }}>
                <span>MFE: {data.mfe.toFixed(2)} kcal/mol</span>
                {data.dot_bracket && <span>{data.dot_bracket.length} nt</span>}
              </div>
            )}
          </>
        )}
      </section>

      {/* Binding site duplex detail (client-side, no API needed) */}
      {site && (
        <section
          className="rounded-xl border p-4 sm:p-5"
          style={{ background: "var(--bg-surface)", borderColor: "var(--border)" }}
        >
          <h2 className="text-sm font-semibold mb-3 sm:mb-4" style={{ color: "var(--text-primary)" }}>
            Binding site duplex detail
          </h2>

          {/* Site tabs */}
          {design.sites.length > 1 && (
            <div className="flex flex-wrap gap-1.5 mb-4">
              {design.sites.map((s, i) => (
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
