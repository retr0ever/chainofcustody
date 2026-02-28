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
        <label className="text-xs font-medium text-zinc-600">{label}</label>
        <span className="text-xs font-mono font-semibold text-zinc-800 bg-zinc-100 rounded px-1.5 py-0.5">
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
        className="w-full accent-sky-600"
      />
      <p className="text-xs text-zinc-400">{hint}</p>
    </div>
  );
}

interface ParameterControlsProps {
  targetThreshold: number;
  coverThreshold: number;
  maxMirnas: number;
  onTargetThresholdChange: (v: number) => void;
  onCoverThresholdChange: (v: number) => void;
  onMaxMirnasChange: (v: number) => void;
}

export default function ParameterControls({
  targetThreshold,
  coverThreshold,
  maxMirnas,
  onTargetThresholdChange,
  onCoverThresholdChange,
  onMaxMirnasChange,
}: ParameterControlsProps) {
  return (
    <div className="flex flex-col gap-4">
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
        format={(v) => `≥ ${v} RPM`}
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
  );
}
