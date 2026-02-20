"use client";

interface AlignmentStatsProps {
  totalReads: number;
  mappedReads: number;
  mappingRate: number;
  properlyPaired?: number;
  meanCoverage?: number;
}

export function AlignmentStats({
  totalReads,
  mappedReads,
  mappingRate,
  properlyPaired,
  meanCoverage,
}: AlignmentStatsProps) {
  const stats = [
    { label: "Total Reads", value: totalReads.toLocaleString() },
    { label: "Mapped Reads", value: mappedReads.toLocaleString() },
    { label: "Mapping Rate", value: `${(mappingRate * 100).toFixed(1)}%` },
    ...(properlyPaired !== undefined
      ? [{ label: "Properly Paired", value: properlyPaired.toLocaleString() }]
      : []),
    ...(meanCoverage !== undefined
      ? [{ label: "Mean Coverage", value: `${meanCoverage.toFixed(1)}x` }]
      : []),
  ];

  return (
    <div className="rounded-xl border border-slate-700/50 bg-slate-800/40 p-4">
      <h3 className="font-medium text-slate-200 mb-3">Alignment Statistics</h3>
      <div className="grid grid-cols-2 gap-3">
        {stats.map(({ label, value }) => (
          <div key={label} className="rounded-lg bg-slate-900/50 p-3 text-center">
            <p className="text-lg font-semibold text-cyan-400">{value}</p>
            <p className="text-xs text-slate-500 mt-1">{label}</p>
          </div>
        ))}
      </div>
    </div>
  );
}
