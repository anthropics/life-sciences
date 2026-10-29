#!/usr/bin/env python3
# Copyright 2026 Anthropic PBC
# SPDX-License-Identifier: Apache-2.0
"""Draws the figure specs a report writer returned (the spec format is FIG_SCHEMA below; render and render_all draw
them). Needs matplotlib; without it no figure is drawn and the report keeps its text.

    python3 ft_figures.py --specs <file holding a JSON list of figure specs> --out <folder>

Prints one JSON object as its last line: {"ok": true, "drawn": {filename: caption}, "problem": ""}."""

import argparse
import json
import math
import re
import sys
from pathlib import Path

FIG_SCHEMA = {
    "type": "array",
    "items": {
        "type": "object",
        "properties": {
            "filename": {"type": "string"},
            "kind": {"type": "string", "enum": ["forest", "bar", "scatter", "waterfall", "curve", "heatmap"]},
            "title": {"type": "string"},
            "caption": {"type": "string"},
            "x_label": {"type": "string"},
            "y_label": {"type": "string"},
            "series": {
                "type": "array",
                "items": {
                    "type": "object",
                    "properties": {
                        "label": {"type": "string"},
                        "type": {"type": "string", "enum": ["points", "line", "bars", "band"]},
                        "labels": {"type": "array", "items": {"type": "string"}},
                        "x": {"type": "array", "items": {"type": "number"}},
                        "y": {"type": "array", "items": {"type": "number"}},
                        "lo": {"type": "array", "items": {"type": "number"}},
                        "hi": {"type": "array", "items": {"type": "number"}},
                        "matrix": {"type": "array", "items": {"type": "array", "items": {"type": "number"}}},
                    },
                    "required": ["label"],
                },
            },
            "reference_lines": {
                "type": "array",
                "items": {
                    "type": "object",
                    "properties": {
                        "axis": {"type": "string", "enum": ["x", "y"]},
                        "value": {"type": "number"},
                        "label": {"type": "string"},
                    },
                    "required": ["axis", "value", "label"],
                },
            },
        },
        "required": ["filename", "kind", "title", "caption", "series"],
    },
}

PALETTE = ["#1F7A8C", "#C08A2D", "#5B6776", "#8C3B1F", "#3B7A57", "#6B4C9A"]


def _safe_name(fn: str) -> str:
    fn = re.sub(r"[^A-Za-z0-9_.-]", "_", Path(fn).name)
    return fn if fn.lower().endswith(".png") else fn + ".png"


def _num(v, default=float("nan")):
    try:
        f = float(v)
        return f if math.isfinite(f) else default
    except Exception:
        return default


def _plain(spec: dict) -> dict:
    """A copy of the spec with "$" escaped in every text that is drawn: matplotlib reads text between two "$" as math."""
    def t(v):
        if isinstance(v, list):
            return [t(x) for x in v]
        return v.replace(r"\$", "$").replace("$", r"\$") if isinstance(v, str) else v  # a "$" already escaped stays escaped once

    def texts(d, keys):
        return dict(d, **{k: t(d[k]) for k in keys if k in d}) if isinstance(d, dict) else d
    out = texts(spec, ("title", "x_label", "y_label"))
    for key, keys in (("series", ("label", "labels")), ("reference_lines", ("label",))):
        if isinstance(spec.get(key), list):
            out[key] = [texts(x, keys) for x in spec[key]]
    return out


def render(spec: dict, out_dir: Path) -> Path | None:
    """Draw one spec; returns the file path or None (a malformed spec is skipped, never fatal)."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    spec = _plain(spec)
    out_dir.mkdir(parents=True, exist_ok=True)
    fn = out_dir / _safe_name(spec["filename"])
    kind = spec["kind"]
    series = [s for s in spec.get("series", []) if isinstance(s, dict)]
    if not series:
        return None
    fig, ax = plt.subplots(figsize=(10, 5.4) if spec.get("kind") == "waterfall" else (8.5, 5.0), dpi=150)
    try:
        if kind == "forest":
            s0 = series[0]
            labels = list(s0.get("labels") or [])
            y = [_num(v) for v in (s0.get("y") or [])]
            n = min(len(labels), len(y))
            if n == 0:
                raise ValueError("forest needs labels and y")
            lo = [_num(v) for v in (s0.get("lo") or [])][:n]
            hi = [_num(v) for v in (s0.get("hi") or [])][:n]
            ys = list(range(n))[::-1]
            if len(lo) == n and len(hi) == n:
                ax.hlines(ys, lo, hi, color=PALETTE[0], lw=2)
            ax.plot(y[:n], ys, "o", color=PALETTE[0], ms=7, label=s0.get("label", ""))
            for extra, col in zip(series[1:], PALETTE[1:]):
                ey = [_num(v) for v in (extra.get("y") or [])][:n]
                if len(ey) == n:
                    elo = [_num(v) for v in (extra.get("lo") or [])][:n]
                    ehi = [_num(v) for v in (extra.get("hi") or [])][:n]
                    if len(elo) == n and len(ehi) == n:
                        ax.hlines([v - 0.18 for v in ys], elo, ehi, color=col, lw=2, alpha=0.8)
                    ax.plot(ey, [v - 0.18 for v in ys], "D", color=col, ms=6, label=extra.get("label", ""))
            ax.set_yticks(ys)
            ax.set_yticklabels(labels[:n])
            ax.set_xlabel(spec.get("x_label", ""))
        elif kind == "bar":
            s0 = series[0]
            labels = list(s0.get("labels") or [])
            k = len([s for s in series if s.get("y")])
            width = 0.8 / max(k, 1)
            xs = list(range(len(labels)))
            j = 0
            for s, col in zip(series, PALETTE):
                y = [_num(v) for v in (s.get("y") or [])][: len(labels)]
                if not y:
                    continue
                off = [x + (j - (k - 1) / 2) * width for x in xs[: len(y)]]
                lo = [_num(v) for v in (s.get("lo") or [])][: len(y)]
                hi = [_num(v) for v in (s.get("hi") or [])][: len(y)]
                err = None
                if len(lo) == len(y) and len(hi) == len(y):
                    err = [[max(a - b, 0) for a, b in zip(y, lo)], [max(b - a, 0) for a, b in zip(y, hi)]]
                ax.bar(off, y, width=width * 0.95, color=col, label=s.get("label", ""), yerr=err, capsize=3)
                j += 1
            ax.set_xticks(xs)
            ax.set_xticklabels(labels, rotation=20, ha="right")
            ax.set_ylabel(spec.get("y_label", ""))
        elif kind == "scatter":
            for s, col in zip(series, PALETTE):
                x = [_num(v) for v in (s.get("x") or [])]
                y = [_num(v) for v in (s.get("y") or [])]
                n = min(len(x), len(y))
                if n == 0:
                    continue
                if s.get("type") == "line":
                    ax.plot(x[:n], y[:n], "-", color=col, lw=2, label=s.get("label", ""))
                else:
                    ax.plot(x[:n], y[:n], "o", color=col, ms=7, label=s.get("label", ""))
                    lo = [_num(v) for v in (s.get("lo") or [])][:n]
                    hi = [_num(v) for v in (s.get("hi") or [])][:n]
                    if len(lo) == n and len(hi) == n:
                        ax.vlines(x[:n], lo, hi, color=col, lw=1.2, alpha=0.7)
                    for lab, xi, yi in zip(s.get("labels") or [], x, y):
                        ax.annotate(str(lab), (xi, yi), textcoords="offset points", xytext=(5, 4), fontsize=8)
            ax.set_xlabel(spec.get("x_label", ""))
            ax.set_ylabel(spec.get("y_label", ""))
        elif kind == "waterfall":
            s0 = series[0]
            labels = list(s0.get("labels") or [])
            deltas = [_num(v) for v in (s0.get("y") or [])]
            n = min(len(labels), len(deltas))
            if n < 2:
                raise ValueError("waterfall needs a start value and at least one step")
            level = deltas[0]
            ax.bar(0, level, color=PALETTE[2])
            ax.text(0, level, f"{level:.3g}", ha="center", va="bottom", fontsize=8)
            for i in range(1, n):
                d = deltas[i]
                col = PALETTE[4] if d >= 0 else PALETTE[3]
                ax.bar(i, d, bottom=level, color=col)
                ax.text(i, level + d, f"{d:+.3g}", ha="center", va="bottom", fontsize=8)
                level += d
            ax.bar(n, level, color=PALETTE[0])
            ax.text(n, level, f"{level:.3g}", ha="center", va="bottom", fontsize=8)
            import textwrap as _tw

            ax.set_xticks(range(n + 1))
            ax.set_xticklabels(
                ["\n".join(_tw.wrap(str(l), 16)) for l in labels[:n]] + ["result"], rotation=0, fontsize=8
            )
            ax.set_ylabel(spec.get("y_label", ""))
        elif kind == "curve":
            for s, col in zip(series, PALETTE):
                x = [_num(v) for v in (s.get("x") or [])]
                y = [_num(v) for v in (s.get("y") or [])]
                n = min(len(x), len(y))
                if n == 0:
                    continue
                lo = [_num(v) for v in (s.get("lo") or [])][:n]
                hi = [_num(v) for v in (s.get("hi") or [])][:n]
                if s.get("type") == "band" and len(lo) == n and len(hi) == n:
                    ax.fill_between(x[:n], lo, hi, color=col, alpha=0.2, label=s.get("label", ""))
                else:
                    ax.plot(x[:n], y[:n], "-", color=col, lw=2, label=s.get("label", ""))
                    if len(lo) == n and len(hi) == n:
                        ax.fill_between(x[:n], lo, hi, color=col, alpha=0.15)
            ax.set_xlabel(spec.get("x_label", ""))
            ax.set_ylabel(spec.get("y_label", ""))
        elif kind == "heatmap":
            s0 = series[0]
            M = [[_num(v) for v in row] for row in (s0.get("matrix") or []) if isinstance(row, list)]
            if not M:
                raise ValueError("heatmap needs a matrix")
            im = ax.imshow(M, cmap="viridis", aspect="auto")
            fig.colorbar(im, ax=ax, label=spec.get("y_label", ""))
            xl = s0.get("labels") or []
            yl = series[1].get("labels") if len(series) > 1 else []
            if xl:
                ax.set_xticks(range(len(xl)))
                ax.set_xticklabels(xl, rotation=20, ha="right")
            if yl:
                ax.set_yticks(range(len(yl)))
                ax.set_yticklabels(yl)
            for i, row in enumerate(M):
                for j, v in enumerate(row):
                    if math.isfinite(v):
                        ax.text(j, i, f"{v:.2g}", ha="center", va="center", fontsize=7, color="white")
            ax.set_xlabel(spec.get("x_label", ""))
        import textwrap

        for j, rl in enumerate(spec.get("reference_lines") or []):
            v = _num(rl.get("value"))
            if not math.isfinite(v):
                continue
            lab = textwrap.shorten(str(rl.get("label", "")), 48)
            if rl.get("axis") == "x":
                ax.axvline(v, ls="--", color="#8C3B1F", lw=1.2)
                # stagger successive labels down the axis so nearby lines do not overprint each other
                ax.annotate(
                    f" {lab}",
                    xy=(v, 1.0 - 0.05 * j),
                    xycoords=("data", "axes fraction"),
                    va="top",
                    fontsize=7.5,
                    color="#8C3B1F",
                    bbox=dict(boxstyle="round,pad=0.15", fc="white", ec="none", alpha=0.8),
                )
            else:
                ax.axhline(v, ls="--", color="#8C3B1F", lw=1.2)
                ax.annotate(
                    f" {lab}",
                    xy=(0.02 + 0.3 * (j % 3), v),
                    xycoords=("axes fraction", "data"),
                    va="bottom",
                    fontsize=7.5,
                    color="#8C3B1F",
                    bbox=dict(boxstyle="round,pad=0.15", fc="white", ec="none", alpha=0.8),
                )
        ax.set_title("\n".join(textwrap.wrap(str(spec.get("title", "")), 80)), fontsize=10.5, loc="left")
        if any(s.get("label") for s in series) and kind != "heatmap":
            ax.legend(fontsize=8, frameon=False)
        ax.grid(alpha=0.25)
        fig.tight_layout()
        fig.savefig(fn)
        return fn
    except Exception:
        return None
    finally:
        plt.close(fig)


def render_all(specs: list, out_dir: Path) -> dict:
    """Render every spec; returns {filename: caption} for the ones that drew."""
    done = {}
    for sp in specs or []:
        try:
            f = render(sp, out_dir)
        except Exception:
            f = None
        if f:
            done[f.name] = sp.get("caption", "")
    return done


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--specs", required=True)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    try:
        specs = json.loads(Path(a.specs).read_text(encoding="utf-8"))
    except (OSError, ValueError) as e:
        print(json.dumps({"ok": False, "drawn": {}, "problem": f"cannot read {a.specs}: {e}"}))
        return 1
    try:
        import matplotlib  # noqa: F401
    except ImportError:
        print(json.dumps({"ok": True, "drawn": {}, "problem": "matplotlib is not installed, so no figure was drawn"}))
        return 0
    drawn = render_all(specs if isinstance(specs, list) else [], Path(a.out))
    print(json.dumps({"ok": True, "drawn": drawn, "problem": ""}))
    return 0


if __name__ == "__main__":
    sys.exit(main())
