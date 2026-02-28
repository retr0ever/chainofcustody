"use client";

import { useState, useEffect } from "react";
import type { MirnaData } from "./types";

let cachedData: MirnaData | null = null;
let pendingFetch: Promise<MirnaData> | null = null;

async function fetchMirnaData(): Promise<MirnaData> {
  if (cachedData) return cachedData;
  if (pendingFetch) return pendingFetch;

  pendingFetch = fetch("/mirna_data.json")
    .then((res) => {
      if (!res.ok) throw new Error(`Failed to load mirna_data.json: ${res.status}`);
      return res.json() as Promise<MirnaData>;
    })
    .then((data) => {
      cachedData = data;
      pendingFetch = null;
      return data;
    });

  return pendingFetch;
}

export interface UseMirnaDataResult {
  data: MirnaData | null;
  loading: boolean;
  error: string | null;
}

export function useMirnaData(): UseMirnaDataResult {
  const [data, setData] = useState<MirnaData | null>(null);
  const [loading, setLoading] = useState(true);
  const [error, setError] = useState<string | null>(null);

  useEffect(() => {
    let cancelled = false;

    fetchMirnaData()
      .then((d) => {
        if (!cancelled) {
          setData(d);
          setLoading(false);
        }
      })
      .catch((err) => {
        if (!cancelled) {
          setError(String(err));
          setLoading(false);
        }
      });

    return () => {
      cancelled = true;
    };
  }, []);

  return { data, loading, error };
}
