"use client";

interface SliderProps {
  label: string;
  hint: string;
  value: number;
  min: number;
  max: number;
  step: number;
  format?: (v: number) => string;
  onChange: (v: number) => void;
}

function Slider({ label, hint, value, min, max, step, format, onChange }: SliderProps) {
  const display = format ? format(value) : String(value);
  return (
    <div className="flex flex-col gap-1">
      <div className="flex justify-between items-baseline">
        <label className="text-xs font-medium" style={{ color: "var(--text-secondary)" }}>{label}</label>
        <span
          className="text-xs font-mono font-semibold rounded px-1.5 py-0.5"
          style={{ background: "var(--bg-inset)", color: "var(--text-primary)" }}
        >
          {display}
        </span>
      </div>
      <input
        type="range"
        min={min}
        max={max}
        step={step}
        value={value}
        onChange={(e) => onChange(Number(e.target.value))}
        className="w-full"
        style={{ accentColor: "var(--primary)" }}
      />
      <p className="text-xs" style={{ color: "var(--text-tertiary)" }}>{hint}</p>
    </div>
  );
}

interface TextInputProps {
  label: string;
  hint: string;
  value: string;
  placeholder: string;
  onChange: (v: string) => void;
  onAction?: () => void;
  actionLabel?: string;
  loading?: boolean;
}

function TextInput({ label, hint, value, placeholder, onChange, onAction, actionLabel, loading }: TextInputProps) {
  return (
    <div className="flex flex-col gap-1.5">
      <label className="text-xs font-medium" style={{ color: "var(--text-secondary)" }}>{label}</label>
      <div className="flex gap-2">
        <input
          type="text"
          value={value}
          placeholder={placeholder}
          onChange={(e) => onChange(e.target.value.toUpperCase().replace(/[^A-Z0-9]/g, ""))}
          className="flex-1 px-3 py-2 text-xs font-mono rounded-md border transition-colors outline-none"
          style={{ 
            background: "var(--bg-inset)", 
            borderColor: "var(--border)",
            color: "var(--text-primary)"
          }}
        />
        {onAction && (
          <button
            onClick={onAction}
            disabled={loading || !value}
            className="px-3 py-2 text-xs font-medium rounded-md transition-colors cursor-pointer disabled:opacity-50 disabled:cursor-not-allowed shrink-0"
            style={{ 
              background: "var(--primary-bg)", 
              color: "var(--primary)",
              border: "1px solid var(--primary)"
            }}
          >
            {loading ? "..." : actionLabel}
          </button>
        )}
      </div>
      <p className="text-xs" style={{ color: "var(--text-tertiary)" }}>{hint}</p>
    </div>
  );
}

interface ParameterControlsProps {
  targetThreshold: number;
  coverThreshold: number;
  maxMirnas: number;
  utr5Seq: string;
  cdsSeq: string;
  geneSymbol: string;
  fetchingCds: boolean;
  optimizationMode: "manual" | "optimize";
  optimizing: boolean;
  onTargetThresholdChange: (v: number) => void;
  onCoverThresholdChange: (v: number) => void;
  onMaxMirnasChange: (v: number) => void;
  onUtr5SeqChange: (v: string) => void;
  onCdsSeqChange: (v: string) => void;
  onGeneSymbolChange: (v: string) => void;
  onFetchCds: () => void;
  onOptimizationModeChange: (v: "manual" | "optimize") => void;
  onRunOptimization: () => void;
}

export default function ParameterControls({
  targetThreshold,
  coverThreshold,
  maxMirnas,
  utr5Seq,
  cdsSeq,
  geneSymbol,
  fetchingCds,
  optimizationMode,
  optimizing,
  onTargetThresholdChange,
  onCoverThresholdChange,
  onMaxMirnasChange,
  onUtr5SeqChange,
  onCdsSeqChange,
  onGeneSymbolChange,
  onFetchCds,
  onOptimizationModeChange,
  onRunOptimization,
}: ParameterControlsProps) {
  return (
    <div className="flex flex-col gap-5">
      <div className="space-y-4">
        <div className="flex items-center justify-between">
          <h3 className="text-[10px] font-bold uppercase tracking-wider text-[var(--text-tertiary)]">5&apos; UTR Mode</h3>
          <div className="flex bg-[var(--bg-inset)] p-0.5 rounded-md border border-[var(--border)]">
            <button
              onClick={() => onOptimizationModeChange("manual")}
              className={`px-2 py-1 text-[10px] rounded-sm transition-colors ${optimizationMode === "manual" ? "bg-[var(--bg-surface)] shadow-sm text-[var(--primary)] font-semibold" : "text-[var(--text-tertiary)]"}`}
            >
              Manual
            </button>
            <button
              onClick={() => onOptimizationModeChange("optimize")}
              className={`px-2 py-1 text-[10px] rounded-sm transition-colors ${optimizationMode === "optimize" ? "bg-[var(--bg-surface)] shadow-sm text-[var(--primary)] font-semibold" : "text-[var(--text-tertiary)]"}`}
            >
              Optimise
            </button>
          </div>
        </div>

        {optimizationMode === "manual" ? (
          <TextInput
            label="Custom 5' UTR"
            hint="Sequence upstream of the CDS"
            value={utr5Seq}
            placeholder="e.g. GCAUAC..."
            onChange={onUtr5SeqChange}
          />
        ) : (
          <div className="space-y-3">
            <button
              onClick={onRunOptimization}
              disabled={optimizing || !geneSymbol}
              className="w-full py-2 text-xs font-bold rounded-md transition-all flex items-center justify-center gap-2 disabled:opacity-50"
              style={{ 
                background: "var(--primary)", 
                color: "white",
                boxShadow: "0 2px 4px rgba(0,0,0,0.1)"
              }}
            >
              {optimizing ? (
                <>
                  <svg className="animate-spin h-3 w-3 text-white" fill="none" viewBox="0 0 24 24">
                    <circle className="opacity-25" cx="12" cy="12" r="10" stroke="currentColor" strokeWidth="4" />
                    <path className="opacity-75" fill="currentColor" d="M4 12a8 8 0 018-8v4a4 4 0 00-4 4H4z" />
                  </svg>
                  Running Genetic Algorithm...
                </>
              ) : "Run 5' UTR Optimisation"}
            </button>
            <p className="text-[10px] text-[var(--text-tertiary)] text-center leading-relaxed">
              Uses a Genetic Algorithm (NSGA-III) to evolve a 5&apos;UTR that maximizes Translation Efficiency while maintaining stability.
            </p>
          </div>
        )}
        
        <div className="pt-2 border-t" style={{ borderColor: "var(--border)" }} />
        
        <h3 className="text-[10px] font-bold uppercase tracking-wider text-[var(--text-tertiary)]">Coding Sequence</h3>
        <TextInput
          label="Gene Symbol"
          hint="Fetch canonical CDS from Ensembl"
          value={geneSymbol}
          placeholder="e.g. POU5F1"
          onChange={onGeneSymbolChange}
          onAction={onFetchCds}
          actionLabel="Fetch CDS"
          loading={fetchingCds}
        />

        <TextInput
          label="CDS (Coding Sequence)"
          hint="The protein-coding region"
          value={cdsSeq}
          placeholder="e.g. AUG..."
          onChange={onCdsSeqChange}
        />
      </div>

      <div className="space-y-4 pt-4 border-t" style={{ borderColor: "var(--border)" }}>
        <h3 className="text-[10px] font-bold uppercase tracking-wider text-[var(--text-tertiary)]">3&apos; UTR Algorithm</h3>
        <Slider
          label="Target silence threshold"
          hint="Max mean RPM in target cells for a miRNA to be a candidate"
          value={targetThreshold}
          min={1}
          max={200}
          step={1}
          format={(v) => `< ${v} RPM`}
          onChange={onTargetThresholdChange}
        />
        <Slider
          label="Off-target coverage threshold"
          hint="Min mean RPM in an off-target cell to count it as covered"
          value={coverThreshold}
          min={50}
          max={5000}
          step={50}
          format={(v) => `\u2265 ${v} RPM`}
          onChange={onCoverThresholdChange}
        />
        <Slider
          label="Max miRNAs to select"
          hint="Algorithm stops after this many miRNAs regardless of coverage"
          value={maxMirnas}
          min={1}
          max={40}
          step={1}
          onChange={onMaxMirnasChange}
        />
      </div>
    </div>
  );
}
