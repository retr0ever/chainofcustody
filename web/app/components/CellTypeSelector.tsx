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
      // Remove from off-targets if now selected as target
      if (offTargetSet.has(cellType)) {
        onOffTargetsChange(offTargets.filter((t) => t !== cellType));
      }
      onTargetsChange(next);
    } else {
      const next = offTargetSet.has(cellType)
        ? offTargets.filter((t) => t !== cellType)
        : [...offTargets, cellType].filter((t) => !targetSet.has(t));
      // Remove from targets if now selected as off-target
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
      <div className="flex rounded-lg border border-zinc-200 overflow-hidden text-sm font-medium">
        <button
          onClick={() => setActivePanel("targets")}
          className={`flex-1 px-3 py-2 flex items-center justify-center gap-2 transition-colors ${
            activePanel === "targets"
              ? "bg-rose-50 text-rose-700 border-r border-zinc-200"
              : "bg-white text-zinc-500 border-r border-zinc-200 hover:bg-zinc-50"
          }`}
        >
          <span className="w-2.5 h-2.5 rounded-full bg-rose-400 inline-block" />
          Targets
          {targets.length > 0 && (
            <span className="ml-1 rounded-full bg-rose-200 text-rose-800 text-xs px-1.5 py-0.5 leading-none">
              {targets.length}
            </span>
          )}
        </button>
        <button
          onClick={() => setActivePanel("offTargets")}
          className={`flex-1 px-3 py-2 flex items-center justify-center gap-2 transition-colors ${
            activePanel === "offTargets"
              ? "bg-sky-50 text-sky-700"
              : "bg-white text-zinc-500 hover:bg-zinc-50"
          }`}
        >
          <span className="w-2.5 h-2.5 rounded-full bg-sky-400 inline-block" />
          Off-targets
          {offTargets.length > 0 && (
            <span className="ml-1 rounded-full bg-sky-200 text-sky-800 text-xs px-1.5 py-0.5 leading-none">
              {offTargets.length}
            </span>
          )}
        </button>
      </div>

      {/* Quick-select: all non-targets as off-targets */}
      {targets.length > 0 && (
        <button
          onClick={selectAllNonTargetsAsOffTargets}
          className="w-full flex items-center justify-center gap-1.5 rounded-lg border border-sky-200 bg-sky-50 px-3 py-2 text-xs font-medium text-sky-700 hover:bg-sky-100 transition-colors"
        >
          <span className="w-2 h-2 rounded-full bg-sky-400 flex-shrink-0" />
          Set all non-target cell types as off-targets
          <span className="ml-auto rounded-full bg-sky-200 text-sky-800 px-1.5 py-0.5 leading-none">
            {allCellTypes.length - targets.length}
          </span>
        </button>
      )}

      {/* Search */}
      <div className="relative">
        <svg
          className="absolute left-3 top-1/2 -translate-y-1/2 text-zinc-400 w-4 h-4 pointer-events-none"
          fill="none"
          viewBox="0 0 24 24"
          stroke="currentColor"
        >
          <path
            strokeLinecap="round"
            strokeLinejoin="round"
            strokeWidth={2}
            d="M21 21l-4.35-4.35M17 11A6 6 0 1 1 5 11a6 6 0 0 1 12 0z"
          />
        </svg>
        <input
          type="text"
          placeholder="Search cell types…"
          value={search}
          onChange={(e) => setSearch(e.target.value)}
          className="w-full pl-9 pr-3 py-2 text-sm rounded-lg border border-zinc-200 bg-white focus:outline-none focus:ring-2 focus:ring-sky-300"
        />
      </div>

      {/* List */}
      <div className="border border-zinc-200 rounded-lg overflow-hidden">
        <div className="flex items-center justify-between px-3 py-1.5 bg-zinc-50 border-b border-zinc-200 text-xs text-zinc-500">
          <span>
            {filtered.length} cell type{filtered.length !== 1 ? "s" : ""}
          </span>
          <button
            onClick={clearAll}
            className="text-zinc-400 hover:text-zinc-700 transition-colors"
          >
            Clear {activePanel === "targets" ? "targets" : "off-targets"}
          </button>
        </div>
        <ul className="max-h-72 overflow-y-auto divide-y divide-zinc-100">
          {filtered.map((ct) => {
            const state = getState(ct);
            const isActive =
              (activePanel === "targets" && state === "target") ||
              (activePanel === "offTargets" && state === "offtarget");
            const isOther =
              (activePanel === "targets" && state === "offtarget") ||
              (activePanel === "offTargets" && state === "target");

            return (
              <li key={ct}>
                <button
                  onClick={() => toggle(ct)}
                  className={`w-full text-left px-3 py-2 text-sm flex items-center gap-2 transition-colors ${
                    isActive
                      ? activePanel === "targets"
                        ? "bg-rose-50 text-rose-800"
                        : "bg-sky-50 text-sky-800"
                      : isOther
                      ? "opacity-40 cursor-default"
                      : "hover:bg-zinc-50 text-zinc-700"
                  }`}
                >
                  <span
                    className={`w-2 h-2 rounded-full flex-shrink-0 ${
                      state === "target"
                        ? "bg-rose-400"
                        : state === "offtarget"
                        ? "bg-sky-400"
                        : "bg-zinc-200"
                    }`}
                  />
                  <span className="truncate">{ct.replace(/_/g, " ")}</span>
                </button>
              </li>
            );
          })}
          {filtered.length === 0 && (
            <li className="px-3 py-4 text-sm text-zinc-400 text-center">
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
              className="inline-flex items-center gap-1 px-2 py-0.5 rounded-full bg-rose-100 text-rose-700 text-xs"
            >
              {ct.replace(/_/g, " ")}
              <button
                onClick={() => onTargetsChange(targets.filter((t) => t !== ct))}
                className="hover:text-rose-900 leading-none"
              >
                ×
              </button>
            </span>
          ))}
          {offTargets.map((ct) => (
            <span
              key={`o-${ct}`}
              className="inline-flex items-center gap-1 px-2 py-0.5 rounded-full bg-sky-100 text-sky-700 text-xs"
            >
              {ct.replace(/_/g, " ")}
              <button
                onClick={() =>
                  onOffTargetsChange(offTargets.filter((t) => t !== ct))
                }
                className="hover:text-sky-900 leading-none"
              >
                ×
              </button>
            </span>
          ))}
        </div>
      )}
    </div>
  );
}
