"use client";

import { useEffect, useRef } from "react";

interface MolViewer3DProps {
  smiles: string;
  width?: number;
  height?: number;
}

export function MolViewer3D({ smiles, width = 400, height = 300 }: MolViewer3DProps) {
  const containerRef = useRef<HTMLDivElement>(null);

  useEffect(() => {
    let mounted = true;

    async function render() {
      if (!containerRef.current || !smiles) return;

      try {
        const $3Dmol = await import("3dmol");
        if (!mounted || !containerRef.current) return;

        containerRef.current.innerHTML = "";
        const viewer = $3Dmol.createViewer(containerRef.current, {
          backgroundColor: "0x1e293b",
        });

        // Use SDF from a conversion (or inline SMILES if 3Dmol supports it)
        // 3Dmol can load from SMILES via addModel with "smi" format in some builds
        // Fallback: show a placeholder
        try {
          viewer.addModel(smiles, "smi");
          viewer.setStyle({}, { stick: { colorscheme: "cyanCarbon" } });
          viewer.zoomTo();
          viewer.render();
        } catch {
          containerRef.current.innerHTML =
            '<p class="text-slate-500 text-sm p-4">3D rendering requires SDF/PDB format. Use the mol_visualizer tool to generate 3D HTML.</p>';
        }
      } catch {
        if (containerRef.current) {
          containerRef.current.innerHTML =
            '<p class="text-slate-500 text-sm p-4">3Dmol.js not available</p>';
        }
      }
    }

    render();
    return () => { mounted = false; };
  }, [smiles]);

  return (
    <div
      ref={containerRef}
      style={{ width, height }}
      className="rounded-lg border border-slate-700 overflow-hidden"
    />
  );
}
