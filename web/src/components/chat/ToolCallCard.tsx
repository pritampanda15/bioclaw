"use client";

import { useState } from "react";
import { ChevronDown, ChevronRight, Wrench, CheckCircle, XCircle } from "lucide-react";

export interface ToolCall {
  name: string;
  input: Record<string, unknown>;
  id: string;
  result?: {
    success: boolean;
    summary: string;
    data?: Record<string, unknown>;
    error?: string;
  };
}

export function ToolCallCard({ toolCall }: { toolCall: ToolCall }) {
  const [expanded, setExpanded] = useState(false);
  const hasResult = !!toolCall.result;
  const success = toolCall.result?.success;

  return (
    <div className="ml-11 rounded-xl border border-slate-700/50 bg-slate-800/40 overflow-hidden">
      <button
        onClick={() => setExpanded(!expanded)}
        className="flex w-full items-center gap-3 px-4 py-2.5 text-sm hover:bg-slate-800/60 transition-colors"
      >
        <Wrench className="h-4 w-4 text-amber-400 shrink-0" />
        <span className="font-medium text-slate-200">{toolCall.name}</span>
        <span className="text-xs text-slate-500 truncate flex-1 text-left">
          {Object.entries(toolCall.input)
            .map(([k, v]) => `${k}=${typeof v === "string" ? v : JSON.stringify(v)}`)
            .join(", ")
            .slice(0, 80)}
        </span>
        {hasResult &&
          (success ? (
            <CheckCircle className="h-4 w-4 text-emerald-400 shrink-0" />
          ) : (
            <XCircle className="h-4 w-4 text-red-400 shrink-0" />
          ))}
        {expanded ? (
          <ChevronDown className="h-4 w-4 text-slate-500 shrink-0" />
        ) : (
          <ChevronRight className="h-4 w-4 text-slate-500 shrink-0" />
        )}
      </button>

      {expanded && (
        <div className="border-t border-slate-700/50 px-4 py-3 space-y-3 text-xs">
          <div>
            <p className="font-medium text-slate-400 mb-1">Input</p>
            <pre className="rounded-lg bg-slate-900 p-3 overflow-x-auto text-slate-300">
              {JSON.stringify(toolCall.input, null, 2)}
            </pre>
          </div>
          {toolCall.result && (
            <div>
              <p className="font-medium text-slate-400 mb-1">Result</p>
              {toolCall.result.summary && (
                <p className="text-slate-300 mb-2">{toolCall.result.summary}</p>
              )}
              {toolCall.result.error && (
                <p className="text-red-400 mb-2">{toolCall.result.error}</p>
              )}
              {toolCall.result.data && Object.keys(toolCall.result.data).length > 0 && (
                <pre className="rounded-lg bg-slate-900 p-3 overflow-x-auto text-slate-300 max-h-64 overflow-y-auto">
                  {JSON.stringify(toolCall.result.data, null, 2)}
                </pre>
              )}
            </div>
          )}
        </div>
      )}
    </div>
  );
}
