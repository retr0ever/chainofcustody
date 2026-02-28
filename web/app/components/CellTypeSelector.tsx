"use client";

import { useState, useMemo } from "react";

interface CellTypeSelectorProps {
  allCellTypes: string[];
  targets: string[];
  offTargets: string[];
  onTargetsChange: (targets: string[]) => void;
  onOffTargetsChange: (offTargets: string[]) => void;
}

type Panel = "targets" | "offTargets";

export default function CellTypeSelector({
  allCellTypes,
  targets,
  offTargets,
  onTargetsChange,
  onOffTargetsChange,
}: CellTypeSelectorProps) {
  const [search, setSearch] = useState("");
  const [activePanel, setActivePanel] = useState<Panel>("targets");

  const targetSet = useMemo(() => new Set(targets), [targets]);
  const offTargetSet = useMemo(() => new Set(offTargets), [offTargets]);

  const filtered = useMemo(() => {
    const q = search.toLowerCase();
    return allCellTypes.filter((ct) => ct.toLowerCase().includes(q));
  }, [allCellTypes, search]);

  function toggle(cellType: string) {
    if (activePanel === "targets") {
      const next = targetSet.has(cellType)
        ? targets.filter((t) => t !== cellType)
        : [...targets, cellType].filter((t) => !offTargetSet.has(t));
      if (offTargetSet.has(cellType)) {
        onOffTargetsChange(offTargets.filter((t) => t !== cellType));
      }
      onTargetsChange(next);
    } else {
      const next = offTargetSet.has(cellType)
        ? offTargets.filter((t) => t !== cellType)
        : [...offTargets, cellType].filter((t) => !targetSet.has(t));
      if (targetSet.has(cellType)) {
        onTargetsChange(targets.filter((t) => t !== cellType));
      }
      onOffTargetsChange(next);
    }
  }

  function clearAll() {
    if (activePanel === "targets") onTargetsChange([]);
    else onOffTargetsChange([]);
  }

  function selectAllNonTargetsAsOffTargets() {
    onOffTargetsChange(allCellTypes.filter((ct) => !targetSet.has(ct)));
  }

  function getState(ct: string): "target" | "offtarget" | "none" {
    if (targetSet.has(ct)) return "target";
    if (offTargetSet.has(ct)) return "offtarget";
    return "none";
  }

  return (
    <div className="flex flex-col gap-3">
      {/* Panel toggle */}
      <div
        className="flex rounded-lg border overflow-hidden text-sm font-medium"
        style={{ borderColor: "var(--border)" }}
      >
        <button
          onClick={() => setActivePanel("targets")}
          className="flex-1 px-3 py-2 flex items-center justify-center gap-2 transition-colors cursor-pointer"
          style={{
            background: activePanel === "targets" ? "var(--red-bg)" : "transparent",
            color: activePanel === "targets" ? "var(--red)" : "var(--text-secondary)",
            borderRight: "1px solid var(--border)",
          }}
        >
          <span className="w-2.5 h-2.5 rounded-full inline-block" style={{ background: "var(--red)" }} />
          Targets
          {targets.length > 0 && (
            <span
              className="ml-1 rounded-full text-xs px-1.5 py-0.5 leading-none"
              style={{ background: "var(--red-bg)", color: "var(--red)" }}
            >
              {targets.length}
            </span>
          )}
        </button>
        <button
          onClick={() => setActivePanel("offTargets")}
          className="flex-1 px-3 py-2 flex items-center justify-center gap-2 transition-colors cursor-pointer"
          style={{
            background: activePanel === "offTargets" ? "var(--primary-bg)" : "transparent",
            color: activePanel === "offTargets" ? "var(--primary)" : "var(--text-secondary)",
          }}
        >
          <span className="w-2.5 h-2.5 rounded-full inline-block" style={{ background: "var(--primary)" }} />
          Off-targets
          {offTargets.length > 0 && (
            <span
              className="ml-1 rounded-full text-xs px-1.5 py-0.5 leading-none"
              style={{ background: "var(--primary-bg)", color: "var(--primary)" }}
            >
              {offTargets.length}
            </span>
          )}
        </button>
      </div>

      {/* Quick-select */}
      {targets.length > 0 && (
        <button
          onClick={selectAllNonTargetsAsOffTargets}
          className="w-full flex items-center justify-center gap-1.5 rounded-lg border px-3 py-2 text-xs font-medium transition-colors cursor-pointer"
          style={{
            borderColor: "var(--border)",
            background: "var(--primary-bg)",
            color: "var(--primary)",
          }}
        >
          <span className="w-2 h-2 rounded-full flex-shrink-0" style={{ background: "var(--primary)" }} />
          Set all non-target cell types as off-targets
          <span
            className="ml-auto rounded-full px-1.5 py-0.5 leading-none"
            style={{ background: "var(--bg-raised)", color: "var(--text-secondary)" }}
          >
            {allCellTypes.length - targets.length}
          </span>
        </button>
      )}

      {/* Search */}
      <div className="relative">
        <svg
          className="absolute left-3 top-1/2 -translate-y-1/2 w-4 h-4 pointer-events-none"
          style={{ color: "var(--text-tertiary)" }}
          fill="none"
          viewBox="0 0 24 24"
          stroke="currentColor"
        >
          <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M21 21l-4.35-4.35M17 11A6 6 0 1 1 5 11a6 6 0 0 1 12 0z" />
        </svg>
        <input
          type="text"
          placeholder="Search cell types..."
          value={search}
          onChange={(e) => setSearch(e.target.value)}
          className="w-full pl-9 pr-3 py-2 text-sm rounded-lg border focus:outline-none"
          style={{
            background: "var(--bg-inset)",
            borderColor: "var(--border)",
            color: "var(--text-primary)",
          }}
        />
      </div>

      {/* List */}
      <div className="rounded-lg border overflow-hidden" style={{ borderColor: "var(--border)" }}>
        <div
          className="flex items-center justify-between px-3 py-1.5 border-b text-xs"
          style={{ background: "var(--bg-inset)", borderColor: "var(--border)", color: "var(--text-secondary)" }}
        >
          <span>{filtered.length} cell type{filtered.length !== 1 ? "s" : ""}</span>
          <button
            onClick={clearAll}
            className="transition-colors cursor-pointer"
            style={{ color: "var(--text-tertiary)" }}
          >
            Clear {activePanel === "targets" ? "targets" : "off-targets"}
          </button>
        </div>
        <ul className="max-h-72 overflow-y-auto">
          {filtered.map((ct) => {
            const state = getState(ct);
            const isActive =
              (activePanel === "targets" && state === "target") ||
              (activePanel === "offTargets" && state === "offtarget");
            const isOther =
              (activePanel === "targets" && state === "offtarget") ||
              (activePanel === "offTargets" && state === "target");

            return (
              <li key={ct} style={{ borderBottom: "1px solid var(--border)" }}>
                <button
                  onClick={() => toggle(ct)}
                  className="w-full text-left px-3 py-2.5 sm:py-2 text-sm flex items-center gap-2 transition-colors cursor-pointer"
                  style={{
                    background: isActive
                      ? activePanel === "targets" ? "var(--red-bg)" : "var(--primary-bg)"
                      : "transparent",
                    color: isActive
                      ? activePanel === "targets" ? "var(--red)" : "var(--primary)"
                      : isOther
                      ? "var(--text-tertiary)"
                      : "var(--text-primary)",
                    opacity: isOther ? 0.4 : 1,
                  }}
                >
                  <span
                    className="w-2 h-2 rounded-full flex-shrink-0"
                    style={{
                      background: state === "target" ? "var(--red)" : state === "offtarget" ? "var(--primary)" : "var(--border-strong)",
                    }}
                  />
                  <span className="truncate">{ct.replace(/_/g, " ")}</span>
                </button>
              </li>
            );
          })}
          {filtered.length === 0 && (
            <li className="px-3 py-4 text-sm text-center" style={{ color: "var(--text-tertiary)" }}>
              No cell types match your search.
            </li>
          )}
        </ul>
      </div>

      {/* Selected chips */}
      {(targets.length > 0 || offTargets.length > 0) && (
        <div className="flex flex-wrap gap-1.5">
          {targets.map((ct) => (
            <span
              key={`t-${ct}`}
              className="inline-flex items-center gap-1.5 px-2.5 py-1 rounded-full text-xs"
              style={{ background: "var(--red-bg)", color: "var(--red)" }}
            >
              {ct.replace(/_/g, " ")}
              <button
                onClick={() => onTargetsChange(targets.filter((t) => t !== ct))}
                className="leading-none cursor-pointer w-4 h-4 flex items-center justify-center rounded-full"
                style={{ background: "rgba(220, 38, 38, 0.15)" }}
              >
                <svg viewBox="0 0 12 12" fill="currentColor" className="w-2.5 h-2.5">
                  <path d="M3.05 3.05a.5.5 0 01.7 0L6 5.29l2.25-2.24a.5.5 0 01.7.7L6.71 6l2.24 2.25a.5.5 0 01-.7.7L6 6.71 3.75 8.95a.5.5 0 01-.7-.7L5.29 6 3.05 3.75a.5.5 0 010-.7z" />
                </svg>
              </button>
            </span>
          ))}
          {offTargets.map((ct) => (
            <span
              key={`o-${ct}`}
              className="inline-flex items-center gap-1.5 px-2.5 py-1 rounded-full text-xs"
              style={{ background: "var(--primary-bg)", color: "var(--primary)" }}
            >
              {ct.replace(/_/g, " ")}
              <button
                onClick={() => onOffTargetsChange(offTargets.filter((t) => t !== ct))}
                className="leading-none cursor-pointer w-4 h-4 flex items-center justify-center rounded-full"
                style={{ background: "rgba(74, 122, 148, 0.15)" }}
              >
                <svg viewBox="0 0 12 12" fill="currentColor" className="w-2.5 h-2.5">
                  <path d="M3.05 3.05a.5.5 0 01.7 0L6 5.29l2.25-2.24a.5.5 0 01.7.7L6.71 6l2.24 2.25a.5.5 0 01-.7.7L6 6.71 3.75 8.95a.5.5 0 01-.7-.7L5.29 6 3.05 3.75a.5.5 0 010-.7z" />
                </svg>
              </button>
            </span>
          ))}
        </div>
      )}
    </div>
  );
}
