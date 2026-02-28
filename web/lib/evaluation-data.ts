"use client";

import { useState, useEffect } from "react";
import type { CandidateEvaluation } from "./types";

let cached: CandidateEvaluation[] | null = null;
let pending: Promise<CandidateEvaluation[]> | null = null;

async function fetchCandidates(): Promise<CandidateEvaluation[]> {
  if (cached) return cached;
  if (pending) return pending;

  pending = fetch("/mock/candidates.json")
    .then((res) => {
      if (!res.ok) throw new Error(`Failed to load candidates: ${res.status}`);
      return res.json();
    })
    .then((data: CandidateEvaluation[]) => {
      cached = data;
      pending = null;
      return data;
    })
    .catch((err) => {
      pending = null;
      throw err;
    });

  return pending;
}

export function useCandidateData() {
  const [candidates, setCandidates] = useState<CandidateEvaluation[]>(cached ?? []);
  const [loading, setLoading] = useState(!cached);
  const [error, setError] = useState<string | null>(null);

  useEffect(() => {
    if (cached) {
      setCandidates(cached);
      setLoading(false);
      return;
    }

    fetchCandidates()
      .then((data) => {
        setCandidates(data);
        setLoading(false);
      })
      .catch((err) => {
        setError(err.message);
        setLoading(false);
      });
  }, []);

  return { candidates, loading, error };
}
