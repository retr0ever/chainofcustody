"use client";

import { createContext, useContext, useState, useCallback, type ReactNode } from "react";
import type { CandidateEvaluation } from "./types";
import { useCandidateData } from "./evaluation-data";

interface CandidatesContextValue {
  candidates: CandidateEvaluation[];
  loading: boolean;
  error: string | null;
  compareSet: Set<string>;
  toggleCompare: (id: string) => void;
  clearCompare: () => void;
  getCandidate: (id: string) => CandidateEvaluation | undefined;
}

const CandidatesContext = createContext<CandidatesContextValue | null>(null);

export function CandidatesProvider({ children }: { children: ReactNode }) {
  const { candidates, loading, error } = useCandidateData();
  const [compareSet, setCompareSet] = useState<Set<string>>(new Set());

  const toggleCompare = useCallback((id: string) => {
    setCompareSet((prev) => {
      const next = new Set(prev);
      if (next.has(id)) {
        next.delete(id);
      } else if (next.size < 4) {
        next.add(id);
      }
      return next;
    });
  }, []);

  const clearCompare = useCallback(() => {
    setCompareSet(new Set());
  }, []);

  const getCandidate = useCallback(
    (id: string) => candidates.find((c) => c.id === id),
    [candidates]
  );

  return (
    <CandidatesContext value={{
      candidates,
      loading,
      error,
      compareSet,
      toggleCompare,
      clearCompare,
      getCandidate,
    }}>
      {children}
    </CandidatesContext>
  );
}

export function useCandidates(): CandidatesContextValue {
  const ctx = useContext(CandidatesContext);
  if (!ctx) throw new Error("useCandidates must be used within CandidatesProvider");
  return ctx;
}
