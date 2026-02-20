"use client";

import { CheckCircle, XCircle, AlertTriangle } from "lucide-react";

interface MolecularProperties {
  formula?: string;
  molecular_weight?: number;
  logp?: number;
  hbd?: number;
  hba?: number;
  tpsa?: number;
  rotatable_bonds?: number;
  qed?: number;
  lipinski?: { pass: boolean; violations: number };
  pains_alerts?: string[];
}

export function PropertyTable({ properties }: { properties: MolecularProperties }) {
  const rows = [
    { label: "Formula", value: properties.formula },
    { label: "MW", value: properties.molecular_weight?.toFixed(2) },
    { label: "LogP", value: properties.logp?.toFixed(2) },
    { label: "HBD", value: properties.hbd },
    { label: "HBA", value: properties.hba },
    { label: "TPSA", value: properties.tpsa?.toFixed(2) },
    { label: "Rotatable Bonds", value: properties.rotatable_bonds },
    { label: "QED", value: properties.qed?.toFixed(3) },
  ].filter((r) => r.value !== undefined);

  return (
    <div className="rounded-xl border border-slate-700/50 bg-slate-800/40 overflow-hidden">
      <table className="w-full text-sm">
        <tbody>
          {rows.map(({ label, value }) => (
            <tr key={label} className="border-b border-slate-700/30 last:border-0">
              <td className="px-4 py-2 font-medium text-slate-400">{label}</td>
              <td className="px-4 py-2 text-slate-200">{value}</td>
            </tr>
          ))}
          {properties.lipinski && (
            <tr className="border-b border-slate-700/30">
              <td className="px-4 py-2 font-medium text-slate-400">Lipinski Ro5</td>
              <td className="px-4 py-2 flex items-center gap-2">
                {properties.lipinski.pass ? (
                  <>
                    <CheckCircle className="h-4 w-4 text-emerald-400" />
                    <span className="text-emerald-400">Pass</span>
                  </>
                ) : (
                  <>
                    <XCircle className="h-4 w-4 text-red-400" />
                    <span className="text-red-400">{properties.lipinski.violations} violations</span>
                  </>
                )}
              </td>
            </tr>
          )}
          {properties.pains_alerts && properties.pains_alerts.length > 0 && (
            <tr>
              <td className="px-4 py-2 font-medium text-slate-400">PAINS</td>
              <td className="px-4 py-2 flex items-center gap-2">
                <AlertTriangle className="h-4 w-4 text-amber-400" />
                <span className="text-amber-400">{properties.pains_alerts.join(", ")}</span>
              </td>
            </tr>
          )}
        </tbody>
      </table>
    </div>
  );
}
