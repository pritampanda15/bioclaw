"use client";

import { useState, useEffect } from "react";
import { Settings, Key, Cpu, Wrench, Save, Check } from "lucide-react";

interface ToolInfo {
  name: string;
  description: string;
  enabled: boolean;
}

interface SettingsState {
  apiKey: string;
  model: string;
  tools: ToolInfo[];
}

const MODELS = [
  { id: "claude-sonnet-4-20250514", label: "Claude Sonnet 4" },
  { id: "claude-opus-4-20250514", label: "Claude Opus 4" },
  { id: "gpt-4o", label: "GPT-4o" },
  { id: "gpt-4o-mini", label: "GPT-4o Mini" },
];

export default function SettingsPage() {
  const [settings, setSettings] = useState<SettingsState>({
    apiKey: "",
    model: "claude-sonnet-4-20250514",
    tools: [],
  });
  const [saved, setSaved] = useState(false);
  const [loadingTools, setLoadingTools] = useState(true);

  useEffect(() => {
    // Load saved settings from localStorage
    const stored = localStorage.getItem("bioclaw_settings");
    if (stored) {
      try {
        const parsed = JSON.parse(stored);
        setSettings((prev) => ({ ...prev, ...parsed }));
      } catch {
        // ignore parse errors
      }
    }

    // Fetch available tools
    fetch("http://localhost:8000/api/tools")
      .then((res) => res.json())
      .then((data) => {
        const tools: ToolInfo[] = (Array.isArray(data) ? data : data.tools ?? []).map(
          (t: { name: string; description?: string }) => ({
            name: t.name,
            description: t.description ?? "",
            enabled: true,
          })
        );
        setSettings((prev) => ({
          ...prev,
          tools: tools.map((t) => {
            const existing = prev.tools.find((pt) => pt.name === t.name);
            return existing ? { ...t, enabled: existing.enabled } : t;
          }),
        }));
      })
      .catch(() => {
        // Tools endpoint not available
      })
      .finally(() => setLoadingTools(false));
  }, []);

  const handleSave = () => {
    localStorage.setItem("bioclaw_settings", JSON.stringify(settings));
    setSaved(true);
    setTimeout(() => setSaved(false), 2000);
  };

  const toggleTool = (name: string) => {
    setSettings((prev) => ({
      ...prev,
      tools: prev.tools.map((t) =>
        t.name === name ? { ...t, enabled: !t.enabled } : t
      ),
    }));
  };

  return (
    <div className="flex h-full flex-col">
      {/* Header */}
      <header className="flex items-center justify-between border-b border-slate-800 bg-slate-900/50 px-6 py-4">
        <div className="flex items-center gap-3">
          <Settings className="h-5 w-5 text-cyan-400" />
          <h1 className="text-lg font-semibold text-white">Settings</h1>
        </div>
        <button
          onClick={handleSave}
          className={`flex items-center gap-2 rounded-lg px-4 py-2 text-sm font-medium transition-colors ${
            saved
              ? "bg-green-600 text-white"
              : "bg-cyan-600 text-white hover:bg-cyan-500"
          }`}
        >
          {saved ? (
            <>
              <Check className="h-4 w-4" /> Saved
            </>
          ) : (
            <>
              <Save className="h-4 w-4" /> Save Settings
            </>
          )}
        </button>
      </header>

      <div className="flex-1 overflow-y-auto p-6">
        <div className="mx-auto max-w-2xl space-y-8">
          {/* API Key */}
          <section className="rounded-xl border border-slate-800 bg-slate-900/50 p-6">
            <div className="mb-4 flex items-center gap-2">
              <Key className="h-5 w-5 text-cyan-400" />
              <h2 className="text-base font-semibold text-white">API Key</h2>
            </div>
            <input
              type="password"
              value={settings.apiKey}
              onChange={(e) =>
                setSettings((prev) => ({ ...prev, apiKey: e.target.value }))
              }
              placeholder="Enter your API key..."
              className="w-full rounded-lg border border-slate-700 bg-slate-800 px-4 py-2.5 text-sm text-slate-200 placeholder-slate-500 focus:border-cyan-500 focus:outline-none focus:ring-1 focus:ring-cyan-500"
            />
            <p className="mt-2 text-xs text-slate-500">
              Your API key is stored locally and never sent to third parties.
            </p>
          </section>

          {/* Model Selection */}
          <section className="rounded-xl border border-slate-800 bg-slate-900/50 p-6">
            <div className="mb-4 flex items-center gap-2">
              <Cpu className="h-5 w-5 text-cyan-400" />
              <h2 className="text-base font-semibold text-white">Model</h2>
            </div>
            <div className="grid grid-cols-2 gap-3">
              {MODELS.map((m) => (
                <button
                  key={m.id}
                  onClick={() =>
                    setSettings((prev) => ({ ...prev, model: m.id }))
                  }
                  className={`rounded-lg border px-4 py-3 text-left text-sm transition-colors ${
                    settings.model === m.id
                      ? "border-cyan-500 bg-cyan-500/10 text-cyan-400"
                      : "border-slate-700 bg-slate-800 text-slate-400 hover:border-slate-600"
                  }`}
                >
                  {m.label}
                </button>
              ))}
            </div>
          </section>

          {/* Tool Toggles */}
          <section className="rounded-xl border border-slate-800 bg-slate-900/50 p-6">
            <div className="mb-4 flex items-center gap-2">
              <Wrench className="h-5 w-5 text-cyan-400" />
              <h2 className="text-base font-semibold text-white">Tools</h2>
            </div>

            {loadingTools ? (
              <p className="text-sm text-slate-500">Loading tools...</p>
            ) : settings.tools.length === 0 ? (
              <p className="text-sm text-slate-500">
                No tools available. Make sure the backend is running.
              </p>
            ) : (
              <div className="space-y-3">
                {settings.tools.map((tool) => (
                  <div
                    key={tool.name}
                    className="flex items-center justify-between rounded-lg border border-slate-700 bg-slate-800 px-4 py-3"
                  >
                    <div className="min-w-0 flex-1">
                      <p className="text-sm font-medium text-slate-200">
                        {tool.name}
                      </p>
                      {tool.description && (
                        <p className="mt-0.5 truncate text-xs text-slate-500">
                          {tool.description}
                        </p>
                      )}
                    </div>
                    <button
                      onClick={() => toggleTool(tool.name)}
                      className={`relative ml-4 h-6 w-11 shrink-0 rounded-full transition-colors ${
                        tool.enabled ? "bg-cyan-600" : "bg-slate-600"
                      }`}
                    >
                      <span
                        className={`absolute top-0.5 h-5 w-5 rounded-full bg-white shadow transition-transform ${
                          tool.enabled ? "left-[22px]" : "left-0.5"
                        }`}
                      />
                    </button>
                  </div>
                ))}
              </div>
            )}
          </section>
        </div>
      </div>
    </div>
  );
}
