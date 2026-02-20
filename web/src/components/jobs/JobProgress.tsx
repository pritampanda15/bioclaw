"use client";

import { Clock, CheckCircle, XCircle, Loader2, Ban } from "lucide-react";

export interface Job {
  id: string;
  tool_name: string;
  status: "pending" | "running" | "completed" | "failed" | "cancelled";
  progress: number | null;
  created_at: string;
  started_at?: string;
  completed_at?: string;
  error?: string;
  result?: {
    success: boolean;
    summary: string;
  };
}

const statusConfig = {
  pending: { icon: Clock, color: "text-blue-400", bg: "bg-blue-400/10", label: "Pending" },
  running: { icon: Loader2, color: "text-amber-400", bg: "bg-amber-400/10", label: "Running" },
  completed: { icon: CheckCircle, color: "text-emerald-400", bg: "bg-emerald-400/10", label: "Completed" },
  failed: { icon: XCircle, color: "text-red-400", bg: "bg-red-400/10", label: "Failed" },
  cancelled: { icon: Ban, color: "text-slate-400", bg: "bg-slate-400/10", label: "Cancelled" },
};

export function JobProgress({ job }: { job: Job }) {
  const config = statusConfig[job.status];
  const Icon = config.icon;

  return (
    <div className="rounded-xl border border-slate-700/50 bg-slate-800/40 p-4 space-y-2">
      <div className="flex items-center justify-between">
        <div className="flex items-center gap-3">
          <Icon
            className={`h-5 w-5 ${config.color} ${job.status === "running" ? "animate-spin" : ""}`}
          />
          <div>
            <p className="font-medium text-slate-200">{job.tool_name}</p>
            <p className="text-xs text-slate-500">{job.id.slice(0, 8)}</p>
          </div>
        </div>
        <span className={`rounded-full px-2.5 py-0.5 text-xs font-medium ${config.color} ${config.bg}`}>
          {config.label}
        </span>
      </div>

      {job.status === "running" && job.progress !== null && (
        <div className="h-1.5 rounded-full bg-slate-700 overflow-hidden">
          <div
            className="h-full rounded-full bg-cyan-500 transition-all duration-300"
            style={{ width: `${job.progress}%` }}
          />
        </div>
      )}

      {job.result?.summary && (
        <p className="text-xs text-slate-400 truncate">{job.result.summary}</p>
      )}
      {job.error && <p className="text-xs text-red-400 truncate">{job.error}</p>}

      <p className="text-xs text-slate-600">
        {new Date(job.created_at).toLocaleString()}
      </p>
    </div>
  );
}
