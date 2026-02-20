"use client";

import "./globals.css";
import { useState } from "react";
import Link from "next/link";
import { usePathname } from "next/navigation";
import {
  MessageSquare,
  ListTodo,
  Settings,
  ChevronLeft,
  ChevronRight,
  Dna,
} from "lucide-react";

const navItems = [
  { href: "/", label: "Chat", icon: MessageSquare },
  { href: "/jobs", label: "Jobs", icon: ListTodo },
  { href: "/settings", label: "Settings", icon: Settings },
];

export default function RootLayout({
  children,
}: {
  children: React.ReactNode;
}) {
  const pathname = usePathname();
  const [collapsed, setCollapsed] = useState(false);

  return (
    <html lang="en" className="dark">
      <head>
        <title>BioClaw</title>
        <meta name="description" content="Bioinformatics AI Agent Platform" />
      </head>
      <body className="bg-slate-950 text-slate-200 antialiased">
        <div className="flex h-screen overflow-hidden">
          {/* Sidebar */}
          <aside
            className={`relative flex flex-col border-r border-slate-800 bg-slate-900 transition-all duration-200 ${
              collapsed ? "w-16" : "w-56"
            }`}
          >
            {/* Logo */}
            <div className="flex h-14 items-center gap-2 border-b border-slate-800 px-4">
              <Dna className="h-6 w-6 shrink-0 text-cyan-400" />
              {!collapsed && (
                <span className="text-lg font-bold tracking-tight text-white">
                  BioClaw
                </span>
              )}
            </div>

            {/* Nav links */}
            <nav className="flex-1 space-y-1 p-2">
              {navItems.map(({ href, label, icon: Icon }) => {
                const active =
                  href === "/" ? pathname === "/" : pathname.startsWith(href);
                return (
                  <Link
                    key={href}
                    href={href}
                    className={`flex items-center gap-3 rounded-lg px-3 py-2.5 text-sm font-medium transition-colors ${
                      active
                        ? "bg-cyan-500/10 text-cyan-400"
                        : "text-slate-400 hover:bg-slate-800 hover:text-slate-200"
                    }`}
                  >
                    <Icon className="h-5 w-5 shrink-0" />
                    {!collapsed && <span>{label}</span>}
                  </Link>
                );
              })}
            </nav>

            {/* Collapse toggle */}
            <button
              onClick={() => setCollapsed(!collapsed)}
              className="flex h-10 items-center justify-center border-t border-slate-800 text-slate-500 hover:text-slate-300 transition-colors"
            >
              {collapsed ? (
                <ChevronRight className="h-4 w-4" />
              ) : (
                <ChevronLeft className="h-4 w-4" />
              )}
            </button>
          </aside>

          {/* Main content */}
          <main className="flex-1 overflow-hidden">{children}</main>
        </div>
      </body>
    </html>
  );
}
