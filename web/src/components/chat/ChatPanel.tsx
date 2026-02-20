"use client";

import { useCallback, useEffect, useRef, useState } from "react";
import { Send, RotateCcw, Loader2 } from "lucide-react";
import { MessageBubble, Message } from "./MessageBubble";
import { ToolCallCard, ToolCall } from "./ToolCallCard";

const WS_URL = "ws://localhost:8000/ws/chat";

export function ChatPanel() {
  const [messages, setMessages] = useState<Message[]>([]);
  const [input, setInput] = useState("");
  const [isConnected, setIsConnected] = useState(false);
  const [isLoading, setIsLoading] = useState(false);
  const wsRef = useRef<WebSocket | null>(null);
  const scrollRef = useRef<HTMLDivElement>(null);
  const pendingTextRef = useRef("");
  const pendingToolsRef = useRef<ToolCall[]>([]);

  const scrollToBottom = () => {
    scrollRef.current?.scrollTo({ top: scrollRef.current.scrollHeight, behavior: "smooth" });
  };

  const connect = useCallback(() => {
    const ws = new WebSocket(WS_URL);

    ws.onopen = () => setIsConnected(true);
    ws.onclose = () => {
      setIsConnected(false);
      setTimeout(connect, 3000);
    };
    ws.onerror = () => ws.close();

    ws.onmessage = (evt) => {
      const { type, data } = JSON.parse(evt.data);

      if (type === "text_delta") {
        pendingTextRef.current += data;
        setMessages((prev) => {
          const last = prev[prev.length - 1];
          if (last?.role === "assistant" && !last.toolCall) {
            return [...prev.slice(0, -1), { ...last, content: pendingTextRef.current }];
          }
          return [...prev, { role: "assistant", content: pendingTextRef.current }];
        });
      } else if (type === "tool_call") {
        pendingToolsRef.current.push({ name: data.name, input: data.input, id: data.id });
        setMessages((prev) => [
          ...prev,
          { role: "assistant", content: "", toolCall: { name: data.name, input: data.input, id: data.id } },
        ]);
      } else if (type === "tool_result") {
        setMessages((prev) => {
          const idx = prev.findLastIndex((m) => m.toolCall?.name === data.name && !m.toolCall.result);
          if (idx >= 0) {
            const updated = [...prev];
            updated[idx] = { ...updated[idx], toolCall: { ...updated[idx].toolCall!, result: data.result } };
            return updated;
          }
          return prev;
        });
      } else if (type === "done") {
        setIsLoading(false);
        pendingTextRef.current = "";
        pendingToolsRef.current = [];
      } else if (type === "error") {
        setIsLoading(false);
        setMessages((prev) => [...prev, { role: "assistant", content: `Error: ${data}` }]);
      }
    };

    wsRef.current = ws;
  }, []);

  useEffect(() => {
    connect();
    return () => wsRef.current?.close();
  }, [connect]);

  useEffect(scrollToBottom, [messages]);

  const send = () => {
    const text = input.trim();
    if (!text || !wsRef.current || wsRef.current.readyState !== WebSocket.OPEN) return;
    setMessages((prev) => [...prev, { role: "user", content: text }]);
    setInput("");
    setIsLoading(true);
    pendingTextRef.current = "";
    wsRef.current.send(JSON.stringify({ message: text }));
  };

  const clear = () => {
    wsRef.current?.send(JSON.stringify({ action: "clear" }));
    setMessages([]);
  };

  return (
    <div className="flex h-full flex-col">
      {/* Header */}
      <div className="flex h-14 items-center justify-between border-b border-slate-800 px-6">
        <h1 className="text-lg font-semibold">Chat</h1>
        <div className="flex items-center gap-3">
          <span className={`h-2 w-2 rounded-full ${isConnected ? "bg-emerald-400" : "bg-red-400"}`} />
          <button onClick={clear} className="rounded-lg p-2 text-slate-400 hover:bg-slate-800 hover:text-slate-200">
            <RotateCcw className="h-4 w-4" />
          </button>
        </div>
      </div>

      {/* Messages */}
      <div ref={scrollRef} className="flex-1 overflow-y-auto px-6 py-4 space-y-4">
        {messages.length === 0 && (
          <div className="flex h-full items-center justify-center text-slate-500">
            <div className="text-center">
              <p className="text-xl font-medium mb-2">Welcome to BioClaw</p>
              <p className="text-sm">Ask about molecular properties, drug discovery, or NGS analysis.</p>
            </div>
          </div>
        )}
        {messages.map((msg, i) =>
          msg.toolCall ? (
            <ToolCallCard key={i} toolCall={msg.toolCall} />
          ) : (
            <MessageBubble key={i} message={msg} />
          )
        )}
        {isLoading && (
          <div className="flex items-center gap-2 text-cyan-400">
            <Loader2 className="h-4 w-4 animate-spin" />
            <span className="text-sm">Thinking...</span>
          </div>
        )}
      </div>

      {/* Input */}
      <div className="border-t border-slate-800 p-4">
        <div className="flex items-center gap-3 rounded-xl bg-slate-900 border border-slate-700 px-4 py-2 focus-within:border-cyan-500/50">
          <input
            value={input}
            onChange={(e) => setInput(e.target.value)}
            onKeyDown={(e) => e.key === "Enter" && !e.shiftKey && (e.preventDefault(), send())}
            placeholder="Ask BioClaw..."
            className="flex-1 bg-transparent text-sm outline-none placeholder:text-slate-500"
          />
          <button
            onClick={send}
            disabled={!input.trim() || isLoading}
            className="rounded-lg bg-cyan-600 p-2 text-white transition-colors hover:bg-cyan-500 disabled:opacity-40"
          >
            <Send className="h-4 w-4" />
          </button>
        </div>
      </div>
    </div>
  );
}
