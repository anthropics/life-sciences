# Copyright 2026 Anthropic PBC
# SPDX-License-Identifier: Apache-2.0
"""The report, the claim ledger and the deliverables of the clinical-trial-forecast tools."""

import collections
import hashlib
import json
import re
import shutil
import subprocess
import sys
from pathlib import Path

from ft_base import *  # noqa: F401,F403
from ft_research import *  # noqa: F401,F403


# ---- the report: narrative layout, condense and copy-edit passes -----------------------------------------------------------
CEILINGS = {"summary": 900, "narrative": 3500, "appendix": 2500, "simulation": 5000, "analysis": 1800}  # words incl. tables
FIG_REF = r"figures/[A-Za-z0-9_.-]+\.png"


def ceiling(name):
    return 3000 if name.startswith("chapter_") else CEILINGS.get(name, 4000)


def _scalar(v, prefer=("headline", "pts", "overall")):
    if isinstance(v, bool):
        return None
    if isinstance(v, (int, float)):
        return float(v)
    if isinstance(v, dict):
        for k, x in v.items():
            if any(w in str(k).lower() for w in prefer) and isinstance(x, (int, float)) and not isinstance(x, bool):
                return float(x)
        for x in v.values():
            if isinstance(x, (int, float)) and not isinstance(x, bool):
                return float(x)
            if isinstance(x, dict):
                y = _scalar(x, prefer)
                if y is not None:
                    return y
    return None


def load_bearing_text(P):
    """A deterministic ranking of what carries the forecast's uncertainty (tornado, conditional power,
    Monte Carlo error), given to every writer so words follow load, not mechanics."""
    r = read_results(P / "sim")
    meta = read_json(P / "sim" / "parameter_metadata.json", {}) or {}
    lines = ["LOAD-BEARING RANKING (from sim/results.json; PTS swing = PTS at the parameter's 90th minus 10th percentile, all else sampled):"]
    tor = (r.get("sensitivity") or {}).get("tornado") or [] if isinstance(r.get("sensitivity"), dict) else []
    if isinstance(tor, dict):
        tor = [{"param": k, **v} for k, v in tor.items() if isinstance(v, dict)]
    rows = []
    for t in tor if isinstance(tor, list) else []:
        if not isinstance(t, dict):
            continue
        pid = t.get("param") or t.get("parameter") or t.get("id")
        lo = _scalar(next((t[k] for k in ("pts_at_p10", "at_q10", "p10", "low") if k in t), None))
        hi = _scalar(next((t[k] for k in ("pts_at_p90", "at_q90", "p90", "high") if k in t), None))
        if lo is None or hi is None:
            continue
        rows.append((abs(hi - lo), str(pid), lo, hi))
    rows.sort(reverse=True)
    for span, pid, lo, hi in rows[:12]:
        name = str((meta.get(pid) or {}).get("display_name") or pid) if isinstance(meta, dict) else pid
        lines.append(f"- {pid} ({name}): PTS {float(lo):.3f} → {float(hi):.3f}, swing {100 * span:.1f} pp"
                     + ("  ← load-bearing" if span >= 0.03 else ("  (table row, not a section)" if span < 0.02 else "")))
    cp = _scalar(r.get("conditional_power_at_point_estimates"))
    pts = _scalar(r.get("pts_overall"))
    mcse = _scalar(r.get("mc_se"))
    if cp is not None and pts is not None:
        if cp >= pts:
            lines.append(f"- Aleatoric part: conditional power at point estimates {cp:.3f} vs headline {pts:.3f} — trial luck accounts for ~{100 * (1 - cp):.1f} pp of the shortfall; the rest is epistemic (the effect-size derivation) and operational. Powering at the expected effect is therefore not a driver: one sentence.")
        else:
            lines.append(f"- Conditional power at the point estimates is only {cp:.3f} against a headline of {pts:.3f}: at the central assumptions the test does NOT clear its bar, and the headline probability comes from the epistemic upside (worlds where the true effect is larger than the point estimate). Say this plainly; it is the single most important fact about this forecast.")
    if mcse is not None:
        lines.append(f"- Monte Carlo SE {mcse:.4f}: negligible — one sentence, never a subsection.")
    return "\n".join(lines) if len(lines) > 1 else "(no tornado available — allocate words by the parameter sheet's stated sensitivities)"


def stages_text(P):
    f = read_json(P / "sim" / "flow.json", {}) or {}
    out = []
    for i, st in enumerate((f.get("stages") if isinstance(f, dict) else None) or [], 1):
        if isinstance(st, dict):
            out.append(f"{i}. {st.get('label')} — {st.get('sublabel', '')} · params: {', '.join(map(str, st.get('params_used') or []))}")
    return "\n".join(out) or "(no flow.json — follow the order of sim/sim_report.md's model description)"


def tidy_markdown(md):
    """Clean-up of writer output: inline ' - a - b - c' runs become lists, lists get a blank line before."""
    md = md.replace("\u2010", "-").replace("\u2011", "-").replace("\u2012", "-").replace("\u00ad", "")
    out = []
    for line in md.split("\n"):
        raw = line.rstrip()
        if (raw and not re.match(r"\s*([-*+]\s|#|\||>|!\[|\$\$|\d+\.\s)", raw)
                and len(re.findall(r"(?<=[.:)!?\*])\s+-\s+(?=\S)", raw)) >= 2):
            parts = re.split(r"(?<=[.:)!?\*])\s+-\s+(?=\S)", raw)
            head = parts[0].strip()
            if head:
                out.append(head)
            out.append("")
            out += [f"- {x.strip()}" for x in parts[1:] if x.strip()]
            continue
        out.append(raw)
    fixed = []
    for line in out:
        if line.startswith("- ") and fixed and fixed[-1].strip() and not re.match(r"\s*[-*+]\s", fixed[-1]):
            fixed.append("")
        fixed.append(line)
    return "\n".join(fixed)


def parameter_table_md(meta):
    """The full parameter table from sim/parameter_metadata.json, grouped by category."""
    if not isinstance(meta, dict) or not meta:
        return ""
    groups = {}
    for pid, m in meta.items():
        if isinstance(m, dict):
            groups.setdefault(str(m.get("category") or "Other"), []).append((pid, m))
    out = ["### A.8 Full parameter table", "", "Every parameter of the accepted model in compact form; the evidence anchors behind each are in `model/parameter_metadata.json`, delivered with this report.", ""]

    def c(v):
        return str(v if v is not None else "").replace("|", "\\|").replace("\n", " ")

    for g in sorted(groups):
        out += [f"#### {g}", "", "| id | what it is | distribution | base | 90% range | key | anchor |", "|---|---|---|---|---|---|---|"]
        for pid, m in sorted(groups[g], key=lambda x: (not bool(x[1].get("is_key")), x[0])):
            ev = (m.get("evidence") or [{}])[0] if isinstance(m.get("evidence"), list) and m.get("evidence") else {}
            anchor = str((ev.get("source", "") or ev.get("citation", "")) if isinstance(ev, dict) else "")[:60]
            out.append(f"| `{pid}` | {c(m.get('display_name') or m.get('description'))[:90]} | {c(m.get('distribution'))[:80]} | {c(m.get('base_value'))} | "
                       f"{c(m.get('low'))} – {c(m.get('high'))} | {'★' if m.get('is_key') else ''} | {c(anchor)} |")
        out.append("")
    return "\n".join(out)


def applied(P, f):
    """True when this answer file was already applied to the report parts (a relaunch must not apply it again over
    later edits). Marks it applied otherwise; callers mark only answers they used or deliberately rejected."""
    led = P / "report" / "applied.json"
    seen = read_json(led, {}) or {}
    key = rel(P, f)
    return key in seen


def mark_applied(P, f, what):
    led = P / "report" / "applied.json"
    seen = read_json(led, {}) or {}
    seen[rel(P, f)] = what
    write_json(led, seen)


def parts_dir(P):
    d = P / "report" / "parts"
    d.mkdir(parents=True, exist_ok=True)
    return d


def assemble(P):
    """Title → Summary → narrative → appendix → parameter table → Methods, written
    to report/report.md and report/report.html; the model's figures copied next to it."""
    pj = project(P)
    pd_ = parts_dir(P)

    def rd_(n):
        return tidy_markdown((pd_ / n).read_text(encoding="utf-8").strip()) if (pd_ / n).exists() else ""

    parts = [f"# {pj.get('title', 'Trial forecast')}", f"*Forecast of {(pj.get('focal_trial') or {}).get('nct_id', '')} · generated {now_iso()}*", ""]
    for n in ("summary.md", "narrative.md", "appendix.md"):
        t = rd_(n)
        if t:
            parts += [t, "", "---", ""]
    pt = parameter_table_md(read_json(P / "sim" / "parameter_metadata.json", {}) or {})
    if pt:
        parts += [pt, "", "---", ""]
    t = rd_("methods.md")
    if t:
        parts += [t, ""]
    md = "\n".join(parts)
    if (P / "sim" / "figures").exists():
        copy_pngs(P / "sim" / "figures", P / "report" / "figures")  # a figure the writers drew under the same name stays
    old = (P / "report" / "report.md").read_text(encoding="utf-8") if (P / "report" / "report.md").is_file() else ""
    if old.split("\n")[:1] + old.split("\n")[2:] == md.split("\n")[:1] + md.split("\n")[2:]:
        return old  # same text: keep the file (and its sha256, which the claim ledger is bound to)
    write_text(P / "report" / "report.md", md)
    write_text(P / "report" / "report.html", md_to_html(md, pj.get("title", "Trial forecast")))
    return md


def _escape_raw_html(md):
    """The markdown with every "<" written as &lt;, so HTML that a model or a web page put into the report is shown as
    text in report.html, never rendered or run. Inside code the converter escapes the "&" again; _restore_code puts the
    single-escaped form back there."""
    return md.replace("<", "&lt;")


def _restore_code(body):
    """Inside <code> elements, &amp;lt; (from _escape_raw_html followed by the converter's own escaping) back to &lt;, so
    code shows "<" as written."""
    return re.sub(r"(<code[^>]*>)(.*?)(</code>)", lambda m: m.group(1) + m.group(2).replace("&amp;lt;", "&lt;") + m.group(3),
                  body, flags=re.S)


def _safe_address(value, image=False):
    """True for an address report.html may keep: an http(s) link, an in-page #anchor, or a path with no scheme; for an
    image, only a path under figures/. The address is read as a browser reads it: entities decoded, and whitespace and
    control characters removed."""
    import html as _html
    v = re.sub(r"[\x00-\x20\x7f]+", "", _html.unescape(value))
    if image:
        return v.startswith("figures/") and ".." not in v
    if v.startswith("#") or re.match(r"(?i)https?://", v):
        return True
    return not re.match(r"[^/?#]*:", v)  # a scheme comes before any "/", "?" or "#"


def _drop_unsafe_links(body):
    """Links whose address _safe_address rejects keep their text and lose the link; images it rejects are dropped."""
    body = re.sub(r'<a\b([^>]*)>(.*?)</a>',
                  lambda m: m.group(0) if all(_safe_address(h) for h in re.findall(r'\bhref="([^"]*)"', m.group(1))) else m.group(2),
                  body, flags=re.S | re.I)
    return re.sub(r'<img\b[^>]*>',
                  lambda m: m.group(0) if all(_safe_address(s, image=True) for s in re.findall(r'\bsrc="([^"]*)"', m.group(0))) else "",
                  body, flags=re.I)


def md_to_html(md, title):
    import html as _html
    try:
        import markdown  # optional: without it the page shows the markdown as text
        body = _drop_unsafe_links(_restore_code(markdown.markdown(_escape_raw_html(md), extensions=["tables", "fenced_code", "toc", "footnotes"])))
    except Exception:
        body = "<pre>" + _html.escape(md) + "</pre>"
    css = ("img{max-width:100%;height:auto}body{font:15px/1.5 -apple-system,Segoe UI,Helvetica,Arial,sans-serif;max-width:980px;margin:2rem auto;"
           "padding:0 1rem;color:#1a1a1a}table{border-collapse:collapse;font-size:13px;margin:1rem 0;display:block;overflow-x:auto}th,td{border:1px solid #ccc;"
           "padding:4px 8px;vertical-align:top}th{background:#f2f2f2}h1,h2,h3{margin-top:2rem}code{background:#f4f4f4;padding:1px 4px}")
    return f"<!doctype html><html><head><meta charset='utf-8'><title>{_html.escape(title)}</title><style>{css}</style></head><body>{body}</body></html>"


def ledger_note(P, with_ledger):
    if with_ledger and (P / "report" / "claim_ledger.md").exists():
        return (" Also read report/claim_ledger.md — an independent verification of the previous draft's claims: apply every corrected "
                "wording, add the omitted conditions, and reword any proposal-as-fact clause it flags; do not cite the ledger itself.")
    return ""


def cmd_report_prep(a):
    """The narrative and appendix writers' prompts. Needs the finalized model (sim/flow.json, parameter_metadata.json)."""
    P = Path(a.project)
    if not (P / "sim" / "results.json").exists():
        return fail("sim/results.json is missing: the model phase must publish a model first")
    # a report written on an earlier published model (a relaunch that then finished a later round) is set aside whole,
    # with its claim checks, so the writers and checkers run again on the model now published
    pub = (read_json(P / "sim" / "published_round.json", {}) or {}).get("round")
    made = (read_json(P / "report" / "model_round.json", {}) or {}).get("round")
    superseded = ""
    if (P / "report").exists() and made != pub:
        dst = P / f"report.superseded-{now_iso().replace(':', '').replace('-', '')}"
        (P / "report").rename(dst)
        if (P / "checks" / "claims").exists():
            (P / "checks" / "claims").rename(dst / "claim_checks")
        superseded = rel(P, dst)
    write_json(P / "report" / "model_round.json", {"round": pub})
    out = P / "report" / "prompts"
    lb, stg = load_bearing_text(P), stages_text(P)
    l1 = render_prompt(P, "narrative_writer", out / "narrative.md", stages=stg, load_bearing=lb, ledger_note="")
    l2 = render_prompt(P, "appendix_writer", out / "appendix.md", stages=stg, ledger_note="")
    return done({"ok": True, "prompts": {"narrative": str(out / "narrative.md"), "appendix": str(out / "appendix.md")},
                 "results": {"narrative": str(P / "report" / "out" / "narrative.json"), "appendix": str(P / "report" / "out" / "appendix.json")},
                 "has_flow": (P / "sim" / "flow.json").exists(), "superseded": superseded, "unfilled": sorted(set(l1 + l2))})


def cmd_report_body_apply(a):
    """The writers' answers -> report/parts/narrative.md (its figure specs drawn; references to figures that did not
    draw are removed) and appendix.md; then the summary editor's prompt."""
    P = Path(a.project)
    pd_ = parts_dir(P)
    fig_dir = P / "report" / "figures"
    got, problems = {}, []
    for part, schema in (("narrative", "chapter"), ("appendix", "section")):
        src_json = P / "report" / "out" / f"{part}.json"
        if applied(P, src_json) and (pd_ / f"{part}.md").exists():
            got[part] = {"figures": None, "figure_note": "applied on an earlier launch"}
            continue
        o, errs = load_checked(src_json, schema)
        if errs:
            problems.append(f"{part}: " + "; ".join(errs[:3]))
            continue
        body = o["section_markdown"].strip()
        n_figs, note = 0, ""
        if part == "narrative" and o.get("figure_specs"):
            spec_file = P / "report" / "out" / "figure_specs.json"
            write_json(spec_file, o["figure_specs"])
            # at most 200 s, so the whole command stays under the helper's bash limit (about 295 s for one command)
            try:
                r = subprocess.run([sys.executable, str(SKILL / "scripts" / "ft_figures.py"), "--specs", str(spec_file), "--out", str(fig_dir)],
                                   capture_output=True, text=True, timeout=200)
                res = json.loads((r.stdout.strip().splitlines() or ["{}"])[-1])
            except subprocess.TimeoutExpired:
                res = {"drawn": {}, "problem": "drawing the figures took over 200 seconds, so they were left out"}
            except ValueError:
                res = {"drawn": {}, "problem": (r.stderr or "")[-300:]}
            drawn = res.get("drawn") or {}
            n_figs, note = len(drawn), res.get("problem", "")
            for m in re.findall(r"!\[[^\]]*\]\(figures/([^)]+)\)", body):
                if m not in drawn and not (P / "sim" / "figures" / m).exists():
                    body = re.sub(r"!\[[^\]]*\]\(figures/" + re.escape(m) + r"\)\n?", "", body)
        write_text(pd_ / f"{part}.md", body + "\n")
        mark_applied(P, src_json, "written")
        got[part] = {"figures": n_figs, "figure_note": note}
    if "narrative" not in got:
        return fail("the narrative could not be used: " + "; ".join(problems))
    f = P / "report" / "prompts" / "summary_editor.md"
    left = render_prompt(P, "narrative_editor", f, run_label=f"Forecast run of {now_iso()[:10]}", ledger_note="", load_bearing=load_bearing_text(P))
    return done({"ok": True, "parts": got, "problems": problems, "step_notes": problems, "prompt": str(f),
                 "result_file": str(P / "report" / "out" / "summary_editor.json"), "unfilled": left})


def polish_targets(P):
    return sorted(f for f in parts_dir(P).glob("*.md") if f.name != "methods.md" and not f.name.endswith((".draft.md", ".long.md", ".preledger.md")))


def condense_prompts(P, rnd):
    out = []
    lb = load_bearing_text(P)
    for f in polish_targets(P):
        txt = f.read_text(encoding="utf-8")
        n, cap = len(txt.split()), ceiling(f.stem)
        if n <= int(cap * (1.3 if f.stem == "summary" else 1.1)):
            continue
        title = next((l.lstrip("# ").strip() for l in txt.splitlines() if l.startswith("#")), f.stem)
        pf = P / "report" / "prompts" / f"condense-{rnd}-{f.stem}.md"
        render_prompt(P, "condense_part", pf, title=title, section=txt, max_words=cap, n_words=n, must_cut=max(0, n - cap), load_bearing=lb)
        out.append({"part": f.name, "prompt": str(pf), "result_file": str(P / "report" / "out" / f"condense-{rnd}-{f.stem}.json")})
    return out


def cmd_report_editor_apply(a):
    """The summary editor's answer -> parts/summary.md and parts/methods.md; then the first round of condense prompts
    (for every part over its word ceiling)."""
    P = Path(a.project)
    src_json = P / "report" / "out" / "summary_editor.json"
    o, errs = load_checked(src_json, "summary_editor")
    pd_ = parts_dir(P)
    if applied(P, src_json):
        problem = ""
    elif errs:
        problem = "the summary editor's answer could not be used: " + "; ".join(errs[:3])
    else:
        problem = ""
        write_text(pd_ / "summary.md", o["summary_markdown"].strip() + "\n")
        write_text(pd_ / "methods.md", o["methods_markdown"].strip() + "\n")
        mark_applied(P, src_json, "written")
    return done({"ok": True, "editor_problem": problem, "step_notes": [problem] if problem else [], "condense": condense_prompts(P, 1)})


def cmd_condense_apply(a):
    """Round rnd's condensed parts, each kept only when it passes these checks: not under 30% of the ceiling or 20% of the
    part, keeps at least half the figures, starts with a heading, shorter than before. Then the next round's prompts
    while parts are still over length, for at most 3 rounds."""
    P = Path(a.project)
    pd_ = parts_dir(P)
    log = []
    for f in sorted(x for x in (P / "report" / "out").glob(f"condense-{a.round}-*.json") if not x.name.endswith(".rejected.json")):
        part = f.stem.split("-", 2)[2]
        src_f = pd_ / f"{part}.md"
        if applied(P, f):
            continue
        o, errs = load_checked(f, "polish")
        if errs or not src_f.exists():
            log.append(f"{part}: answer unusable")
            continue
        mark_applied(P, f, "seen")
        src, out = src_f.read_text(encoding="utf-8"), o["section_markdown"].strip()
        cap = ceiling(part)
        n1, n2 = len(src.split()), len(out.split())
        fs, fo = set(re.findall(FIG_REF, src)), set(re.findall(FIG_REF, out))
        if (n2 < max(0.3 * cap, 0.2 * n1) or (fs and len(fo & fs) < 0.5 * len(fs)) or not out.lstrip().startswith("#") or n2 >= n1):
            log.append(f"{part}: rejected (words {n1} -> {n2}); kept the original")
            continue
        if not (pd_ / f"{part}.long.md").exists():
            write_text(pd_ / f"{part}.long.md", src)
        write_text(src_f, out + "\n")
        log.append(f"{part}: {n1} -> {n2} words (ceiling {cap})")
    nxt = condense_prompts(P, a.round + 1) if a.round < 3 else []
    return done({"ok": True, "log": log, "condense": nxt})


def cmd_copyedit_prep(a):
    P = Path(a.project)
    out = []
    for f in polish_targets(P):
        txt = f.read_text(encoding="utf-8")
        if len(txt) < 400:
            continue
        title = next((l.lstrip("# ").strip() for l in txt.splitlines() if l.startswith("#")), f.stem)
        pf = P / "report" / "prompts" / f"copyedit-{f.stem}.md"
        render_prompt(P, "copyedit_section", pf, title=title, section=txt)
        out.append({"part": f.name, "prompt": str(pf), "result_file": str(P / "report" / "out" / f"copyedit-{f.stem}.json")})
    return done({"ok": True, "copyedit": out})


def cmd_copyedit_apply(a):
    """Each copy-edited part kept only when it passes these checks (more than half the original's length, every
    figure kept, not longer than 1.05 times the original or 1.1 times the ceiling); then report/report.md."""
    P = Path(a.project)
    pd_ = parts_dir(P)
    log = []
    for f in sorted(x for x in (P / "report" / "out").glob("copyedit-*.json") if not x.name.endswith(".rejected.json")):
        part = f.stem.split("-", 1)[1]
        src_f = pd_ / f"{part}.md"
        if applied(P, f):
            log.append(f"{part}: applied on an earlier launch")
            continue
        o, errs = load_checked(f, "polish")
        if errs or not src_f.exists():
            log.append(f"{part}: answer unusable")
            continue
        mark_applied(P, f, "seen")
        src, out = src_f.read_text(encoding="utf-8"), o["section_markdown"].strip()
        if len(out) > 0.5 * len(src) and all(m in out for m in re.findall(FIG_REF, src)) and \
                len(out.split()) <= max(int(len(src.split()) * 1.05), int(ceiling(part) * 1.1)):
            write_text(pd_ / f"{part}.draft.md", src)
            write_text(src_f, out + "\n")
            log.append(f"{part}: copy-edited")
        else:
            log.append(f"{part}: dropped a figure, grew past the ceiling or came back too short; kept the draft")
    md = assemble(P)
    return done({"ok": True, "log": log, "words": len(md.split())})


# ---- claim ledger ------------------------------------------------------------------------------------------------------
def check_file(P, i):
    """The i-th claim's (0-based) check file. Claim checks live outside report/, which the checkers may not open."""
    return P / "checks" / "claims" / f"c{i + 1:02d}.json"


def cmd_freeze(a):
    """Records the sha256 of report/report.md and writes report/claims_source.md, the text the extractor reads (the report up to
    its first Appendix heading, at most 90,000 characters)."""
    R = Path(a.project) / "report"
    if not (R / "report.md").is_file():
        return fail("report/report.md is missing")
    md = (R / "report.md").read_text(encoding="utf-8")
    write_text(R / "claims_source.md", re.split(r"\n#+ Appendix", md)[0][:90000])
    digest = hashlib.sha256((R / "report.md").read_bytes()).hexdigest()
    write_text(R / "report.sha256", digest + "\n")
    return done({"ok": True, "report_sha256": digest})


def cmd_claims(a):
    """After extraction: keeps at most --max claims, numbers them C01, C02, ... (check file cNN.json belongs to claim CNN),
    reads sponsor_only strictly, and prints each claim's public wording for the checkers."""
    R = Path(a.project) / "report"
    got = read_json(R / "claims.json")
    claims = [c for c in (got or {}).get("claims") or [] if isinstance(c, dict) and str(c.get("claim") or "").strip()] if isinstance(got, dict) else []
    if not claims:
        return fail("report/claims.json is missing or holds no claims: delete it, then launch the workflow again")
    claims = claims[: a.max]
    for i, c in enumerate(claims):
        c["id"] = f"C{i + 1:02d}"  # renumbered: model-written ids may repeat
        c["sponsor_only"] = truthy(c.get("sponsor_only")) or not str(c.get("public_claim") or "").strip()
    got["claims"] = claims
    if read_json(R / "claims.json") != got:
        write_json(R / "claims.json", got)
    # checks older than this claims.json belong to an earlier report: remove them, so every claim is checked afresh
    stale = [f for f in (Path(a.project) / "checks" / "claims").glob("c*.json") if f.stat().st_mtime < (R / "claims.json").stat().st_mtime - 1]
    for f in stale:
        f.unlink()
    out = []
    for i, c in enumerate(claims):
        if c["sponsor_only"]:
            continue
        pf = Path(a.project) / "checks" / "claims" / f"c{i + 1:02d}.prompt.md"
        render_prompt(Path(a.project), "verify_claim", pf, claim=str(c.get("public_claim") or ""), context=str(c.get("public_context") or ""))
        out.append({"id": c["id"], "prompt": str(pf), "result_file": str(check_file(Path(a.project), i))})
    return done({"ok": True, "claims": len(claims), "sponsor_only": sum(1 for c in claims if c["sponsor_only"]), "checks": out})


def cmd_ledger(a):
    """One ledger row per extracted claim, bound to the sha256 of report/report.md. A sponsor-only claim, or one whose
    check is missing, says so."""
    P = Path(a.project)
    R = P / "report"
    draft = R / "report.md"
    if (R / "ledger_applied.json").is_file() and (R / "claim_ledger.json").is_file():
        # built and applied on an earlier launch: the report has since taken the ledger's edits, so keep the ledger
        old = read_json(R / "claim_ledger.json", {}) or {}
        tally = collections.Counter(cell(c.get("verdict")) for c in old.get("claims") or [])
        return done({"ok": True, "claims": len(old.get("claims") or []), "verdicts": dict(tally), "note": "built on an earlier launch"})
    claims = read_json(R / "claims.json")
    if not draft.is_file() or not isinstance(claims, dict) or not claims.get("claims"):
        return fail("report/report.md is missing, or report/claims.json is missing or empty")
    frozen = (R / "report.sha256").read_text(encoding="utf-8").strip() if (R / "report.sha256").is_file() else ""
    if frozen != hashlib.sha256(draft.read_bytes()).hexdigest():
        return fail("report/report.md changed after its claims were extracted (or freeze was not run). To redo the checks, delete "
                    "report/claims.json and checks/claims/, then launch the workflow again")
    now = now_iso()
    rows = []
    for i, c in enumerate(claims.get("claims") or []):
        base = {"id": c.get("id") or f"c{i + 1:02d}", "section": c.get("section"), "kind": c.get("kind"), "claim": c.get("claim"), "checked": now}
        empty = {"assertions": [], "conditions": [], "examples": [], "precedent_count": "", "corrected_wording": "", "sources": []}
        chk = read_json(check_file(P, i))
        if truthy(c.get("sponsor_only")):
            rows.append({**empty, "verdict": "not checked (rests on the person's request or background material, or is the run's own judgment)",
                         "confidence": "L", "misread_risk": "", **base})
        elif not isinstance(chk, dict) or not chk.get("verdict") or wrong_shape(chk, "verify"):
            rows.append({**empty, "verdict": "could-not-verify (the checker did not finish)", "confidence": "L", "misread_risk": "", **base})
        else:
            rows.append({**{k: v for k, v in chk.items() if not str(k).startswith("_")}, **base})
    digest = hashlib.sha256(draft.read_bytes()).hexdigest()
    write_json(R / "claim_ledger.json", {"report_sha256": digest, "built": now, "claims": rows})
    for r in rows:
        r["verdict"] = cell(r.get("verdict"))
    tally = collections.Counter(r["verdict"] for r in rows)
    table = [{"id": x["id"], "verdict": x["verdict"], "confidence": x.get("confidence"), "claim": cut(x["claim"], 110),
              "conditions": cut("; ".join(map(str, x.get("conditions") or [])), 160),
              "examples": cut("; ".join(f"{e.get('name')} ({e.get('id')})" for e in (x.get("examples") or [])[:3] if isinstance(e, dict)), 160),
              "possible misreading": cut(x.get("misread_risk"), 140),
              "primary sources": cut("; ".join(cell(s.get("primary_source")) for s in x.get("assertions") or []
                                               if isinstance(s, dict) and s.get("source_type") == "primary"), 140)} for x in rows]
    md = (f"# Claim ledger · {now}\n\nIndependent checkers' findings on the forecast report's headline claims. Verdicts: "
          + ", ".join(f"{n} {v}" for v, n in sorted(tally.items())) + ".\n\n"
          + md_table(table, ["id", "verdict", "confidence", "claim", "conditions", "examples", "possible misreading", "primary sources"], 160)
          + "\n\n## Corrected wordings\n\n"
          + "\n".join(f"- **{x['id']}** ({x['verdict']}): {x.get('corrected_wording', '')}"
                      + (f" — proposal-as-fact: {x['proposal_as_fact']}" if x.get("proposal_as_fact") and str(x["proposal_as_fact"]).lower() != "none" else "")
                      for x in rows if x.get("corrected_wording") or (x.get("proposal_as_fact") and str(x["proposal_as_fact"]).lower() != "none")) + "\n")
    write_text(R / "claim_ledger.md", md)
    return done({"ok": True, "claims": len(rows), "verdicts": dict(tally), "report_sha256": digest})


def cmd_ledger_prep(a):
    """The ledger-editing prompt (prompts/apply_ledger.md) for each report part of at least 200 characters."""
    P = Path(a.project)
    R = P / "report"
    led = read_json(R / "claim_ledger.json")
    if not isinstance(led, dict) or not (R / "report.md").is_file():
        return fail("report/claim_ledger.json or report/report.md is missing")
    if hashlib.sha256((R / "report.md").read_bytes()).hexdigest() != led.get("report_sha256"):
        if not (R / "ledger_applied.json").exists():
            return fail("report/report.md changed after the claim ledger was built. To redo the checks, delete report/claims.json "
                        "and checks/claims/, then launch the workflow again")
    ledger = (R / "claim_ledger.md").read_text(encoding="utf-8")
    out = []
    for f in sorted(parts_dir(P).glob("*.md")):
        if f.name.endswith((".draft.md", ".long.md", ".preledger.md")) or len(f.read_text(encoding="utf-8")) < 200:
            continue
        txt = f.read_text(encoding="utf-8")
        title = next((l.lstrip("# ").strip() for l in txt.splitlines() if l.startswith("#")), f.stem)
        pf = R / "prompts" / f"ledger-{f.stem}.md"
        render_prompt(P, "apply_ledger", pf, title=title, section=txt, ledger=ledger)
        out.append({"part": f.name, "prompt": str(pf), "result_file": str(R / "out" / f"ledger-{f.stem}.json")})
    return done({"ok": True, "edits": out})


def cmd_ledger_apply(a):
    """Each ledger-edited part kept only when it passes these checks (every figure kept, within about 20% of
    the checked part's length); then the report is assembled again."""
    P = Path(a.project)
    R = P / "report"
    pd_ = parts_dir(P)
    log = []
    for f in sorted(x for x in (R / "out").glob("ledger-*.json") if not x.name.endswith(".rejected.json")):
        part = f.stem.split("-", 1)[1]
        src_f = pd_ / f"{part}.md"
        if applied(P, f):
            log.append(f"{part}: applied on an earlier launch")
            continue
        o, errs = load_checked(f, "polish")
        if errs or not src_f.exists():
            log.append(f"{part}: answer unusable")
            continue
        mark_applied(P, f, "seen")
        src, out = src_f.read_text(encoding="utf-8"), o["section_markdown"].strip()
        n_src, n_out = len(src.split()), len(out.split())
        figs_kept = all(m in out for m in re.findall(FIG_REF, src))
        if out and figs_kept and 0.8 * n_src <= n_out <= 1.2 * n_src + 40:
            write_text(pd_ / f"{part}.preledger.md", src)
            write_text(src_f, out + "\n")
            log.append(f"{part}: edited")
        else:
            log.append(f"{part}: kept the checked text (words {n_src} -> {n_out}, figures kept: {figs_kept})")
    md = assemble(P)
    write_json(R / "ledger_applied.json", {"applied": now_iso(), "log": log})
    return done({"ok": True, "log": log, "words": len(md.split())})


def cmd_deliver(a):
    """The deliverables: the forecast report (numbered footnotes) as markdown and HTML with its figures, the claim
    ledger, the model (parameter sheet, code, results, the builder's chapter, critiques) and a one-file summary."""
    P = Path(a.project)
    out = Path(a.out)
    R = P / "report"
    if not (R / "report.md").is_file():
        return fail("report/report.md is missing")
    out.mkdir(parents=True, exist_ok=True)
    pj = project(P)
    # trials whose registry facts came from web search (a run with no registry route): the report says so
    web = sorted(f.stem for f in (P / "raw" / "ctgov").glob("NCT*.json")
                 if str((read_json(f, {}) or {}).get("source") or "").startswith("web search"))
    body = (R / "report.md").read_text(encoding="utf-8")
    if web:
        body += ("\n\n## A note on the trial data\n\nThis run had no access to the trial registry, so the registry facts for "
                 f"{len(web)} trial{'s' if len(web) > 1 else ''} ({', '.join(web[:40])}{', …' if len(web) > 40 else ''}) were gathered "
                 "from web search results instead; a fact no result stated was left out.\n")
    md = footnotize(body)
    write_text(out / "forecast_report.md", md)
    write_text(out / "forecast_report.html", md_to_html(md, pj.get("title", "Trial forecast")))
    files = ["forecast_report.md", "forecast_report.html"]
    if (R / "figures").exists():
        copy_pngs(R / "figures", out / "figures", replace=True)
        files.append("figures/")
    if (R / "claim_ledger.md").is_file():
        shutil.copy2(R / "claim_ledger.md", out / "forecast_claim_ledger.md")
        files.append("forecast_claim_ledger.md")
    mdir = out / "model"
    mdir.mkdir(exist_ok=True)
    for f in ("params.json", "simulate.py", "results.json", "sim_report.md", "run_log.txt", "flow.json", "parameter_metadata.json",
              "summary_stats.json", "status.json"):
        if (P / "sim" / f).is_file():
            shutil.copy2(P / "sim" / f, mdir / f)
    for f in sorted((P / "sim").glob("critique_r*.md")):
        shutil.copy2(f, mdir / f.name)
    files.append("model/")
    res = read_results(P / "sim")
    st = read_json(P / "sim" / "status.json", {}) or {}
    led = read_json(R / "claim_ledger.json", {}) or {}
    summary = {"trial": (pj.get("focal_trial") or {}).get("nct_id"), "title": pj.get("title"),
               "success_definition": (pj.get("simulation") or {}).get("success_definition"),
               "pts_overall": res.get("pts_overall"), "mc_se": res.get("mc_se"), "n_sims": res.get("n_sims"),
               "pts_efficacy_only": res.get("pts_efficacy_only"), "conditional_power_at_point_estimates": res.get("conditional_power_at_point_estimates"),
               "reproducibility": res.get("reproducibility"), "model_rounds": st.get("final_round"), "critic_verdict": st.get("verdict"),
               "claim_verdicts": dict(collections.Counter(cell(c.get("verdict")) for c in led.get("claims") or [])),
               "trials_from_web_search": web, "built": now_iso()}
    write_json(out / "forecast_summary.json", summary)
    files.append("forecast_summary.json")
    return done({"ok": True, "deliverables": files, "pts_overall": res.get("pts_overall"), "mc_se": res.get("mc_se"),
                 "critic_verdict": st.get("verdict"), "claim_verdicts": summary["claim_verdicts"], "words": len(md.split())})


# ---- footnotes and deliverables ----------------------------------------------------------------------------------------
# The footnote renderer: bracket citations become numbered footnotes, listed under "Sources"; NCT ids, PMIDs and DOIs
# get links to their pages, and other sources are listed as written.
CIT = re.compile(
    r"(?<!!)\[((?:reported|derived|sponsor claim|secondary|estimate|label|guidance|registry|URL|PMID|DOI|CT\.gov|NCT)\b[^\[\]]*|NCT\d{8}[^\[\]]*|PMID[^\[\]]*|DOI[^\[\]]*|https?://[^\[\]]*)\](?!\()",
    re.I,
)
PLAIN = {"estimate", "derived", "reported", "inferred", "sponsor claim"}


def _key(txt: str) -> str:
    t = txt.strip()
    # strip provenance prefix for the key: "reported: DOI 10.xxxx/yyyy" -> "DOI 10.xxxx/yyyy"
    t2 = re.sub(
        r"^(reported|derived|sponsor claim|secondary|estimate)\s*:\s*",
        "",
        t,
        flags=re.I,
    )
    return t2.strip().rstrip(".;,")


def footnotize(md: str) -> str:
    notes: dict[str, int] = {}
    order: list[str] = []

    def repl(m):
        raw = m.group(1).strip()
        low = raw.lower().rstrip(".")
        if low in PLAIN:  # bare provenance word → plain words
            return {
                "estimate": "(our estimate)",
                "derived": "(derived)",
                "reported": "",
                "inferred": "(inferred)",
                "sponsor claim": "(sponsor claim)",
            }.get(low, "")
        key = _key(raw)
        if not key or len(key) < 4:
            return ""
        if key not in notes:
            notes[key] = len(notes) + 1
            order.append(key)
        n = notes[key]
        prefix = ""
        pm = re.match(
            r"^(derived|sponsor claim|secondary|estimate)\s*:\s*(.*)$",
            raw,
            flags=re.I | re.S,
        )
        if (
            pm
        ):  # "[derived: in-house census]" → "(derived from in-house census)" + footnote
            word = pm.group(1).lower()
            rest = pm.group(2).strip()
            lead = {
                "derived": "derived from",
                "sponsor claim": "sponsor claim,",
                "secondary": "secondary source:",
                "estimate": "estimate based on",
            }[word]
            prefix = (
                f"({lead} {rest[:60]}) "
                if rest and not re.match(r"^(NCT|PMID|DOI|https?)", rest, re.I)
                else f"({word}) "
            )
        return f"{prefix}[^{n}]"

    # "[PMID 1; NCT0…]" → "[PMID 1][NCT0…]" so each source gets its own footnote and link
    md = re.sub(
        r"\[([^\[\]]+;[^\[\]]+)\]",
        lambda m: "".join(f"[{x.strip()}]" for x in m.group(1).split(";"))
        if all(CIT.fullmatch(f"[{x.strip()}]") for x in m.group(1).split(";"))
        else m.group(0),
        md,
    )
    body = CIT.sub(repl, md)
    # tidy: collapse "[^1][^2]" spacing and space-before-punct
    body = re.sub(r"(?<=[\w)\]%.,;])[ \t]+(\[\^\d+\])", r"\1", body)
    body = re.sub(r"(\[\^\d+\])\s+([.,;:])", r"\1\2", body)
    if not order:
        return md

    def link(k: str) -> str:
        if re.match(r"^NCT\d{8}", k):
            return f"[{k}](https://clinicaltrials.gov/study/{k[:11]})"
        if re.match(r"^PMID\s*\d+", k, re.I):
            ids = re.findall(r"\d{6,9}", k)
            return (
                "; ".join(
                    f"[PMID {i}](https://pubmed.ncbi.nlm.nih.gov/{i}/)" for i in ids
                )
                or k
            )
        if re.match(r"^DOI\s*", k, re.I):
            d = re.sub(r"^DOI\s*", "", k, flags=re.I).split()[0]
            return f"[DOI {d}](https://doi.org/{d})"
        if k.lower().startswith("url "):
            return k[4:]
        if k.startswith("http"):
            return k
        return k

    fn = (
        "\n\n## Sources\n\n"
        + "\n".join(f"[^{i}]: {link(k)}" for i, k in enumerate(order, 1))
        + "\n"
    )
    return body + fn



