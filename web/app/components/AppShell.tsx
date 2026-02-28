"use client";

import { useState, useEffect } from "react";
import Image from "next/image";
import { useTheme } from "@/lib/theme";

function useIsMobile(breakpoint = 1024) {
  const [isMobile, setIsMobile] = useState(false);
  useEffect(() => {
    const mq = window.matchMedia(`(max-width: ${breakpoint - 1}px)`);
    setIsMobile(mq.matches);
    const handler = (e: MediaQueryListEvent) => setIsMobile(e.matches);
    mq.addEventListener("change", handler);
    return () => mq.removeEventListener("change", handler);
  }, [breakpoint]);
  return isMobile;
}

export default function AppShell({ children }: { children: React.ReactNode }) {
  const [collapsed, setCollapsed] = useState(false);
  const [mobileOpen, setMobileOpen] = useState(false);
  const { theme, toggle: toggleTheme } = useTheme();
  const isMobile = useIsMobile();

  // Close mobile drawer on route-like changes or resize to desktop
  useEffect(() => {
    if (!isMobile) setMobileOpen(false);
  }, [isMobile]);

  const sidebarVisible = isMobile ? mobileOpen : true;
  const sidebarWidth = isMobile ? 260 : collapsed ? 60 : 220;

  return (
    <div className="flex min-h-screen">
      {/* Mobile backdrop */}
      {isMobile && mobileOpen && (
        <div
          className="fixed inset-0 z-30 bg-black/50"
          onClick={() => setMobileOpen(false)}
        />
      )}

      {/* Mobile hamburger */}
      {isMobile && !mobileOpen && (
        <button
          onClick={() => setMobileOpen(true)}
          className="fixed top-3 left-3 z-20 flex items-center justify-center w-10 h-10 rounded-lg cursor-pointer"
          style={{
            background: "var(--bg-surface)",
            border: "1px solid var(--border)",
            color: "var(--text-secondary)",
          }}
        >
          <svg viewBox="0 0 20 20" fill="currentColor" className="w-5 h-5">
            <path fillRule="evenodd" d="M3 5a1 1 0 011-1h12a1 1 0 110 2H4a1 1 0 01-1-1zM3 10a1 1 0 011-1h12a1 1 0 110 2H4a1 1 0 01-1-1zM3 15a1 1 0 011-1h12a1 1 0 110 2H4a1 1 0 01-1-1z" clipRule="evenodd" />
          </svg>
        </button>
      )}

      {/* Sidebar */}
      <aside
        className={`
          fixed top-0 left-0 h-screen flex flex-col
          border-r transition-all duration-200 ease-out
          ${isMobile ? "z-40" : "z-20"}
          ${isMobile && !mobileOpen ? "-translate-x-full" : "translate-x-0"}
        `}
        style={{
          width: sidebarWidth,
          background: "var(--bg-surface)",
          borderColor: "var(--border)",
        }}
      >
        {/* Brand */}
        <div
          className="flex flex-col items-center justify-center gap-1.5 px-4 py-4 border-b shrink-0"
          style={{ borderColor: "var(--border)" }}
        >
          <Image
            src="/chain-of-custody.png"
            alt="Chain of Custody"
            width={(!isMobile && collapsed) ? 36 : 80}
            height={(!isMobile && collapsed) ? 20 : 44}
            className="shrink-0"
            style={{ objectFit: "contain" }}
          />
          {(isMobile || !collapsed) && (
            <span
              className="text-sm font-semibold tracking-tight"
              style={{ color: "var(--text-primary)" }}
            >
              Chain of Custody
            </span>
          )}
        </div>

        {/* Nav - single item */}
        <nav className="flex-1 py-3 px-2 flex flex-col gap-0.5">
          <div
            className={`flex items-center gap-2.5 rounded-lg px-3 py-2 text-sm font-medium transition-colors ${
              !isMobile && collapsed ? "justify-center" : ""
            }`}
            style={{
              background: "var(--bg-raised)",
              color: "var(--text-primary)",
            }}
          >
            <svg viewBox="0 0 20 20" fill="currentColor" className="w-[18px] h-[18px] shrink-0">
              <path d="M7 3a1 1 0 000 2h6a1 1 0 100-2H7zM4 7a1 1 0 011-1h10a1 1 0 110 2H5a1 1 0 01-1-1zM2 11a2 2 0 012-2h12a2 2 0 012 2v4a2 2 0 01-2 2H4a2 2 0 01-2-2v-4z" />
            </svg>
            {(isMobile || !collapsed) && <span>The Optimizer</span>}
          </div>
        </nav>

        {/* Theme toggle */}
        <button
          onClick={toggleTheme}
          className="flex items-center justify-center gap-2 h-10 border-t shrink-0 transition-colors cursor-pointer"
          style={{
            borderColor: "var(--border)",
            color: "var(--text-tertiary)",
          }}
          onMouseEnter={(e) => (e.currentTarget.style.color = "var(--text-secondary)")}
          onMouseLeave={(e) => (e.currentTarget.style.color = "var(--text-tertiary)")}
        >
          {theme === "dark" ? (
            <svg viewBox="0 0 20 20" fill="currentColor" className="w-4 h-4">
              <path fillRule="evenodd" d="M10 2a1 1 0 011 1v1a1 1 0 11-2 0V3a1 1 0 011-1zm4 8a4 4 0 11-8 0 4 4 0 018 0zm-.464 4.95l.707.707a1 1 0 001.414-1.414l-.707-.707a1 1 0 00-1.414 1.414zm2.12-10.607a1 1 0 010 1.414l-.706.707a1 1 0 11-1.414-1.414l.707-.707a1 1 0 011.414 0zM17 11a1 1 0 100-2h-1a1 1 0 100 2h1zm-7 4a1 1 0 011 1v1a1 1 0 11-2 0v-1a1 1 0 011-1zM5.05 6.464A1 1 0 106.465 5.05l-.708-.707a1 1 0 00-1.414 1.414l.707.707zm1.414 8.486l-.707.707a1 1 0 01-1.414-1.414l.707-.707a1 1 0 011.414 1.414zM4 11a1 1 0 100-2H3a1 1 0 000 2h1z" clipRule="evenodd" />
            </svg>
          ) : (
            <svg viewBox="0 0 20 20" fill="currentColor" className="w-4 h-4">
              <path d="M17.293 13.293A8 8 0 016.707 2.707a8.001 8.001 0 1010.586 10.586z" />
            </svg>
          )}
          {(isMobile || !collapsed) && <span className="text-xs">{theme === "dark" ? "Light mode" : "Dark mode"}</span>}
        </button>

        {/* Collapse toggle (desktop only) / Close button (mobile) */}
        <button
          onClick={() => {
            if (isMobile) setMobileOpen(false);
            else setCollapsed(!collapsed);
          }}
          className="flex items-center justify-center h-10 border-t shrink-0 transition-colors cursor-pointer"
          style={{
            borderColor: "var(--border)",
            color: "var(--text-tertiary)",
          }}
          onMouseEnter={(e) => (e.currentTarget.style.color = "var(--text-secondary)")}
          onMouseLeave={(e) => (e.currentTarget.style.color = "var(--text-tertiary)")}
        >
          {isMobile ? (
            <svg viewBox="0 0 20 20" fill="currentColor" className="w-4 h-4">
              <path fillRule="evenodd" d="M4.293 4.293a1 1 0 011.414 0L10 8.586l4.293-4.293a1 1 0 111.414 1.414L11.414 10l4.293 4.293a1 1 0 01-1.414 1.414L10 11.414l-4.293 4.293a1 1 0 01-1.414-1.414L8.586 10 4.293 5.707a1 1 0 010-1.414z" clipRule="evenodd" />
            </svg>
          ) : (
            <svg
              viewBox="0 0 20 20"
              fill="currentColor"
              className={`w-4 h-4 transition-transform ${collapsed ? "rotate-180" : ""}`}
            >
              <path
                fillRule="evenodd"
                d="M12.707 5.293a1 1 0 010 1.414L9.414 10l3.293 3.293a1 1 0 01-1.414 1.414l-4-4a1 1 0 010-1.414l4-4a1 1 0 011.414 0z"
                clipRule="evenodd"
              />
            </svg>
          )}
        </button>
      </aside>

      {/* Main content */}
      <main
        className={`flex-1 min-w-0 transition-all duration-200 ease-out ${
          isMobile ? "ml-0" : collapsed ? "ml-[60px]" : "ml-[220px]"
        }`}
      >
        {/* Spacer for mobile hamburger */}
        {isMobile && <div className="h-14" />}
        {children}
      </main>
    </div>
  );
}
