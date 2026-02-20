"use client";

import { useEffect, useRef } from "react";

interface Gene {
  gene: string;
  log2_fold_change: number;
  neg_log10_padj: number;
  significant: boolean;
}

interface VolcanoPlotProps {
  genes: Gene[];
  fcThreshold?: number;
  pThreshold?: number;
}

export function VolcanoPlot({ genes, fcThreshold = 1, pThreshold = 1.301 }: VolcanoPlotProps) {
  const containerRef = useRef<HTMLDivElement>(null);

  useEffect(() => {
    let mounted = true;

    async function render() {
      if (!containerRef.current || genes.length === 0) return;

      try {
        const Plotly = await import("plotly.js-dist-min");
        if (!mounted || !containerRef.current) return;

        const sig = genes.filter((g) => g.significant);
        const nonSig = genes.filter((g) => !g.significant);

        const traces = [
          {
            x: nonSig.map((g) => g.log2_fold_change),
            y: nonSig.map((g) => g.neg_log10_padj),
            text: nonSig.map((g) => g.gene),
            mode: "markers" as const,
            type: "scatter" as const,
            name: "Not significant",
            marker: { color: "#64748b", size: 4, opacity: 0.5 },
          },
          {
            x: sig.map((g) => g.log2_fold_change),
            y: sig.map((g) => g.neg_log10_padj),
            text: sig.map((g) => g.gene),
            mode: "markers" as const,
            type: "scatter" as const,
            name: "Significant",
            marker: { color: "#22d3ee", size: 6 },
          },
        ];

        const layout = {
          paper_bgcolor: "rgba(0,0,0,0)",
          plot_bgcolor: "#0f172a",
          font: { color: "#94a3b8", size: 11 },
          xaxis: { title: "log2 Fold Change", gridcolor: "#1e293b", zerolinecolor: "#334155" },
          yaxis: { title: "-log10(p-adj)", gridcolor: "#1e293b" },
          margin: { t: 20, r: 20, b: 50, l: 60 },
          showlegend: true,
          legend: { x: 0.01, y: 0.99 },
          shapes: [
            { type: "line" as const, x0: -fcThreshold, x1: -fcThreshold, y0: 0, y1: 1, yref: "paper" as const, line: { dash: "dash" as const, color: "#475569" } },
            { type: "line" as const, x0: fcThreshold, x1: fcThreshold, y0: 0, y1: 1, yref: "paper" as const, line: { dash: "dash" as const, color: "#475569" } },
            { type: "line" as const, x0: 0, x1: 1, xref: "paper" as const, y0: pThreshold, y1: pThreshold, line: { dash: "dash" as const, color: "#475569" } },
          ],
        };

        Plotly.newPlot(containerRef.current, traces, layout, { responsive: true });
      } catch {
        if (containerRef.current) {
          containerRef.current.innerHTML = '<p class="text-slate-500 text-sm p-4">Plotly not available</p>';
        }
      }
    }

    render();
    return () => { mounted = false; };
  }, [genes, fcThreshold, pThreshold]);

  return (
    <div ref={containerRef} className="rounded-xl border border-slate-700/50 overflow-hidden h-80" />
  );
}
