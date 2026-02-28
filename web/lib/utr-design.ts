/**
 * TypeScript port of chainofcustody/three_prime/generate_utr3.py
 * Generates a 3'UTR from selected miRNA mature sequences.
 */

const COMPLEMENT: Record<string, string> = {
  A: "U", U: "A", G: "C", C: "G",
};

const MISMATCH: Record<string, string> = {
  A: "C", U: "G", G: "U", C: "A",
};

function reverseComplement(seq: string): string {
  return seq.split("").reverse().map((b) => COMPLEMENT[b] ?? b).join("");
}

function createMismatch(bulge: string): string {
  return bulge.split("").map((b) => MISMATCH[b] ?? b).join("");
}

const SPACERS = [
  "aauu", "ucga", "caag", "auac", "gaau",
  "cuua", "uuca", "agcu", "uacg", "gaua",
  "cuac", "acuc", "uguu", "caua", "ucuu", "agau",
];

const STOP_CODON = "UAA";
const LEAD_IN = "gcauac";
const LEAD_OUT = "gauc";
const POLY_A_SIGNAL =
  "CUCAGGUGCAGGCUGCCUAUCAGAAGGUGGUGGCUGGUGUGGCCAAUGCCCUGGCUCACAAAUACCACUGAGAUCUUUUUCCCUCUGCCAAAAAUUAUGGGGACAUCAUGAAGCCCCUUGAGCAUCUGACUUCUGGCUAAUAAAGGAAAUUUAUUUUCAUUGCAAUAGUGUGUUGGAAUUUUUUGUGUCUCUCACUCGGAAGGACAUAUGGGAGGGCAAAUCAUUUAAAACAUCAGAAUGAGUAUUUGGUUUAGAGUUUGGCA";

export type RegionType = "utr5" | "cds" | "stop" | "lead_in" | "site" | "spacer" | "lead_out" | "poly_a";

export interface UtrRegion {
  type: RegionType;
  start: number;
  end: number;
  seq: string;
  /** For "site" regions: which miRNA this site targets */
  mirnaId?: string;
  /** For "site" regions: index into the mirna list */
  mirnaIndex?: number;
}

export interface BindingSite {
  mirnaId: string;
  mirnaSeq: string;
  siteSeq: string;
  seedMatch: string;
  bulgeMismatch: string;
  threePrimeMatch: string;
}

export interface UtrDesignResult {
  fullUtr3: string;
  fullTranscript: string;
  utr5: string;
  cds: string;
  cassette: string;
  sites: BindingSite[];
  regions: UtrRegion[];
  numSites: number;
}

export function generateUtr(
  mirnaSequences: string[],
  mirnaIds: string[],
  numSites = 16,
  utr5Seq = "",
  cdsSeq = "",
): UtrDesignResult {
  if (mirnaSequences.length === 0) {
    return { 
      fullUtr3: "", 
      fullTranscript: "", 
      utr5: "", 
      cds: "", 
      cassette: "", 
      sites: [], 
      regions: [], 
      numSites: 0 
    };
  }

  const bindingSites: BindingSite[] = [];
  for (let idx = 0; idx < mirnaSequences.length; idx++) {
    const mirna = mirnaSequences[idx].toUpperCase().replace(/T/g, "U");
    const rc = reverseComplement(mirna);

    const seedMatch = rc.slice(-8);
    const bulgeRc = rc.slice(-12, -8);
    const threePrimeMatch = rc.slice(0, -12);
    const bulgeMismatch = createMismatch(bulgeRc);

    bindingSites.push({
      mirnaId: mirnaIds[idx] ?? `miRNA-${idx + 1}`,
      mirnaSeq: mirna,
      siteSeq: threePrimeMatch + bulgeMismatch + seedMatch,
      seedMatch,
      bulgeMismatch,
      threePrimeMatch,
    });
  }

  // Build full transcript and track regions
  const regions: UtrRegion[] = [];
  let pos = 0;

  // 5' UTR
  const cleanUtr5 = utr5Seq.toUpperCase().replace(/T/g, "U");
  if (cleanUtr5) {
    regions.push({ type: "utr5", start: pos, end: pos + cleanUtr5.length, seq: cleanUtr5 });
    pos += cleanUtr5.length;
  }

  // CDS
  const cleanCds = cdsSeq.toUpperCase().replace(/T/g, "U");
  if (cleanCds) {
    regions.push({ type: "cds", start: pos, end: pos + cleanCds.length, seq: cleanCds });
    pos += cleanCds.length;
  }

  // Stop codon (only if not already at end of CDS)
  regions.push({ type: "stop", start: pos, end: pos + STOP_CODON.length, seq: STOP_CODON });
  pos += STOP_CODON.length;

  // Lead-in
  regions.push({ type: "lead_in", start: pos, end: pos + LEAD_IN.length, seq: LEAD_IN });
  pos += LEAD_IN.length;

  // Cassette: sites + spacers
  let cassette = "";
  for (let i = 0; i < numSites; i++) {
    const site = bindingSites[i % bindingSites.length];
    regions.push({
      type: "site",
      start: pos,
      end: pos + site.siteSeq.length,
      seq: site.siteSeq,
      mirnaId: site.mirnaId,
      mirnaIndex: i % bindingSites.length,
    });
    cassette += site.siteSeq;
    pos += site.siteSeq.length;

    if (i < numSites - 1) {
      const spacer = SPACERS[i % SPACERS.length];
      regions.push({ type: "spacer", start: pos, end: pos + spacer.length, seq: spacer });
      cassette += spacer;
      pos += spacer.length;
    }
  }

  // Lead-out
  regions.push({ type: "lead_out", start: pos, end: pos + LEAD_OUT.length, seq: LEAD_OUT });
  pos += LEAD_OUT.length;

  // Poly-A signal
  regions.push({ type: "poly_a", start: pos, end: pos + POLY_A_SIGNAL.length, seq: POLY_A_SIGNAL });
  pos += POLY_A_SIGNAL.length;

  const fullUtr3 = `${STOP_CODON}${LEAD_IN}${cassette}${LEAD_OUT}${POLY_A_SIGNAL}`;
  // Use lowercase for utr5 and cds in the full transcript so the backend 
  // structure API doesn't misidentify them as binding sites (which are uppercase)
  const fullTranscript = `${cleanUtr5.toLowerCase()}${cleanCds.toLowerCase()}${fullUtr3}`;

  return { 
    fullUtr3, 
    fullTranscript,
    utr5: cleanUtr5,
    cds: cleanCds,
    cassette, 
    sites: bindingSites, 
    regions, 
    numSites 
  };
}
