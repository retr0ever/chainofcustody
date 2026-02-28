"use client";

import Link from "next/link";

interface NavLinkProps {
  href: string;
  label: string;
  icon: React.ReactNode;
  active?: boolean;
  collapsed?: boolean;
}

export default function NavLink({ href, label, icon, active, collapsed }: NavLinkProps) {
  return (
    <Link
      href={href}
      className={`
        flex items-center gap-2.5 rounded-lg text-sm font-medium
        transition-colors duration-150 relative
        ${collapsed ? "justify-center px-0 py-2.5" : "px-3 py-2"}
      `}
      style={{
        color: active ? "var(--text-primary)" : "var(--text-secondary)",
        background: active ? "var(--primary-bg)" : "transparent",
      }}
      onMouseEnter={(e) => {
        if (!active) {
          e.currentTarget.style.background = "var(--bg-hover)";
          e.currentTarget.style.color = "var(--text-primary)";
        }
      }}
      onMouseLeave={(e) => {
        if (!active) {
          e.currentTarget.style.background = "transparent";
          e.currentTarget.style.color = "var(--text-secondary)";
        }
      }}
    >
      {active && (
        <span
          className="absolute left-0 top-1/2 -translate-y-1/2 w-[3px] h-5 rounded-r-full"
          style={{ background: "var(--primary)" }}
        />
      )}
      <span className="shrink-0" style={{ color: active ? "var(--primary)" : undefined }}>
        {icon}
      </span>
      {!collapsed && <span className="truncate">{label}</span>}
    </Link>
  );
}
