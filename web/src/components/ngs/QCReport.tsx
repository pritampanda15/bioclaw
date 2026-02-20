"use client";

import { CheckCircle, AlertTriangle, XCircle } from "lucide-react";

interface QCModule {
  module: string;
  status: "pass" | "warn" | "fail";
}

interface QCReportProps {
  modules: QCModule[];
  basicStats?: Record<string, string | number>;
}

const statusIcon = {
  pass: <CheckCircle className="h-4 w-4 text-emerald-400" />,
  warn: <AlertTriangle className="h-4 w-4 text-amber-400" />,
  fail: <XCircle className="h-4 w-4 text-red-400" />,
};

export function QCReport({ modules, basicStats }: QCReportProps) {
  return (
    <div className="rounded-xl border border-slate-700/50 bg-slate-800/40 p-4 space-y-4">
      <h3 className="font-medium text-slate-200">Quality Control Report</h3>

      {basicStats && (
        <div className="grid grid-cols-2 gap-2 text-sm">
          {Object.entries(basicStats).map(([key, value]) => (
            <div key={key} className="flex justify-between rounded-lg bg-slate-900/50 px-3 py-1.5">
              <span className="text-slate-400">{key}</span>
              <span className="text-slate-200">{value}</span>
            </div>
          ))}
        </div>
      )}

      <div className="space-y-1">
        {modules.map(({ module, status }) => (
          <div key={module} className="flex items-center gap-2 text-sm">
            {statusIcon[status]}
            <span className="text-slate-300">{module}</span>
          </div>
        ))}
      </div>
    </div>
  );
}
