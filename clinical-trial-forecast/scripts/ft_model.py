# Copyright 2026 Anthropic PBC
# SPDX-License-Identifier: Apache-2.0
"""The model rounds of the clinical-trial-forecast tools: parameter sheet, build, critic, publish, finalize."""

import collections
import hashlib
import json
import re
import shutil
import subprocess
import sys
from pathlib import Path

from ft_base import *  # noqa: F401,F403


# ---- the model: parameter sheet -> build -> critic rounds, publish, finalize -------------------------------------------------
def sha20(obj):
    """A short content hash: sha256 of json.dumps(obj, sort_keys=True, ensure_ascii=False), first 20 hex digits."""
    s = obj if isinstance(obj, (str, bytes)) else json.dumps(obj, sort_keys=True, ensure_ascii=False)
    return hashlib.sha256(s.encode() if isinstance(s, str) else s).hexdigest()[:20]


def params_sha(params):
    """The parameter sheet's hash as prompts/sim_build.md describes it to the builder: sha256 of
    json.dumps(the sheet without its top-level _meta key, sort_keys=True), first 20 hex digits."""
    return hashlib.sha256(json.dumps({k: v for k, v in params.items() if k != "_meta"}, sort_keys=True).encode()).hexdigest()[:20]


def clean(s, n=600):
    s = re.sub(r"[\x00-\x08\x0b-\x1f\x7f]", "", str(s or ""))
    return s if len(s) <= n else s[: n - 1] + "…"


def requests_text(crit, target):
    """The critic's change requests as a numbered, length-capped list, the only critic text a builder sees."""
    if not crit:
        return "(none — first round)"
    rows = sorted((c for c in crit.get("change_requests", []) if isinstance(c, dict) and c.get("target") in target),
                  key=lambda c: c.get("priority", 3))
    if not rows:
        return "(none)"
    return "\n".join(f"{i}. [{c.get('target')} · priority {c.get('priority')}] {clean(c.get('what'))}\n   why: {clean(c.get('why'), 400)}\n"
                     f"   evidence: {clean(c.get('evidence'), 400)}" for i, c in enumerate(rows, 1))


def evidence_sha(P):
    parts = []
    for pat in ("knowledge/*.md", "sim/critique_*.md", "sim/params.json"):
        for f in sorted(P.glob(pat)):
            parts.append(f"{rel(P, f)}:{sha20(f.read_bytes())}")
    for f in sorted((P / "intake").glob("**/*")):
        if f.is_file() and f.suffix in (".txt", ".md"):
            parts.append(f"intake/{f.name}:{sha20(f.read_bytes())}")
    return sha20("\n".join(parts))


def validate_params(pr):
    """Structural checks of the parameter sheet beyond its schema."""
    errs = []
    ids = [x["id"] for x in pr.get("parameters", [])]
    if len(ids) != len(set(ids)):
        errs.append("duplicate parameter ids")
    if not pr["trial"]["arms"]:
        errs.append("no arms")
    if not any(e["role"] in ("primary", "co_primary") for e in pr["trial"]["endpoints"]):
        errs.append("no primary endpoint")
    for x in pr.get("parameters", []):
        if x["kind"] == "epistemic" and x["distribution"]["family"] == "fixed":
            errs.append(f"epistemic parameter {x['id']} has a fixed distribution")
        if not x.get("basis") or not x.get("sources"):
            errs.append(f"parameter {x['id']} lacks basis/sources")
    known = set(ids) | {e["id"] for e in pr["trial"]["endpoints"]} | {a["id"] for a in pr["trial"]["arms"]} | {pr["trial"]["name"]}
    for d in pr.get("dependencies", []):
        for b in d.get("between", []):
            if b not in known and not any(k.lower() in b.lower() for k in known):
                errs.append(f"dependency references unknown id {b}")
    return errs


def sim_dirs(P):
    sd = P / "sim"
    ws = sd / "workspace"
    for d in (sd, ws, sd / "history", ws / "figures"):
        d.mkdir(parents=True, exist_ok=True)
    return sd, ws, sd / "history"


def last_critique(P, rd):
    """The critique of the latest round before rd that has one (None on the first round)."""
    for r in range(rd - 1, 0, -1):
        c = read_json(P / "sim" / "history" / f"round{r:02d}" / "critique.json")
        if isinstance(c, dict):
            return c
    return None


def cmd_sim_prep(a):
    """Round rd: writes the parameter-sheet prompt when this round needs a sheet (round 1, or the last critic asked for
    parameter changes); otherwise says so, and the round goes straight to the build."""
    P = Path(a.project)
    sd, ws, hist = sim_dirs(P)
    rd = a.round
    rdir = hist / f"round{rd:02d}"
    rdir.mkdir(parents=True, exist_ok=True)
    own = read_json(rdir / "critique.json")
    if isinstance(own, dict) and own.get("verdict"):
        # this round finished on an earlier launch: the script moves on (or stops) without redoing it
        return done({"ok": True, "round_done": True, "verdict": own["verdict"], "another_round": own["verdict"] != "accept" and rd < a.rounds,
                     "need_spec": False, "accepted_already": False})
    crit = last_critique(P, rd)
    params = read_json(sd / "params.json")
    if crit and crit.get("verdict") == "accept":
        return done({"ok": True, "round_done": False, "need_spec": False, "accepted_already": True})
    need = params is None or rd == 1 or (crit and crit.get("verdict") in ("revise_params", "revise_both"))
    if not need:
        return done({"ok": True, "round_done": False, "need_spec": False, "accepted_already": False})
    prev = json.dumps(params, ensure_ascii=False)[:60000] if params else "(none yet)"
    f = rdir / "spec.prompt.md"
    left = render_prompt(P, "sim_spec", f, round=rd, previous_params=prev, critic_requests=requests_text(crit, {"params"}),
                         critic_assessment=clean(crit.get("overall_assessment"), 1500) if crit else "(none — first round)",
                         evidence_sha=evidence_sha(P))
    return done({"ok": True, "round_done": False, "need_spec": True, "accepted_already": False, "prompt": str(f), "result_file": str(rdir / "params.json"), "unfilled": left})


def cmd_sim_spec_apply(a):
    """The round's parameter sheet (when it wrote one) -> checks, sim/params.json; then the workspace: params.json,
    brief.md and the builder's prompt."""
    P = Path(a.project)
    sd, ws, hist = sim_dirs(P)
    rd = a.round
    rdir = hist / f"round{rd:02d}"
    S = size_of(a.size)
    if (rdir / "spec.prompt.md").is_file():
        pr, errs = load_checked(rdir / "params.json", "params")
        if errs:
            return fail("the parameter sheet does not match the schema: " + "; ".join(errs[:6]))
        warn = validate_params(pr)
        pr["_meta"] = {"round": rd, "saved": now_iso(), "warnings": warn}
        write_json(rdir / "params.json", pr)
        write_json(sd / "params.json", pr)
    else:
        warn = []
    params = read_json(sd / "params.json")
    if not isinstance(params, dict):
        return fail("sim/params.json is missing")
    shutil.copy2(sd / "params.json", ws / "params.json")
    psha = params_sha(params)
    shutil.copy2(SKILL / "scripts" / "bg.py", ws / "bg.py")  # runs the builder's long commands in the background
    render_prompt(P, "sim_brief", ws / "brief.md", n_sims=S["n_sims"])
    crit = last_critique(P, rd)
    mode = "revise" if (ws / "results.json").exists() and rd > 1 else "build"
    # a builder's result stamped with this round (by sim-build-apply, on an earlier launch) is kept; any other is removed,
    # so the builder runs again (a launch that ended between the build and sim-build-apply rebuilds)
    br = ws / "build_result.json"
    if br.is_file() and (read_json(br, {}) or {}).get("_round") != rd:
        br.unlink()
    f = ws / "build.prompt.md"
    left = render_prompt(P, "sim_build", f, mode=mode, round=rd, n_sims=S["n_sims"], params_sha=psha,
                         critic_requests=requests_text(crit, {"code", "report"}),
                         critic_assessment=clean(crit.get("overall_assessment"), 1500) if crit else "(none — first round)")
    return done({"ok": True, "mode": mode, "params_sha": psha, "parameters": len(params.get("parameters") or []),
                 "warnings": warn, "step_notes": warn[:3], "prompt": str(f), "result_file": str(br), "unfilled": left})


def cmd_sim_build_apply(a):
    """The builder's result -> the round's copies; then the critic's prompt."""
    P = Path(a.project)
    sd, ws, hist = sim_dirs(P)
    rd = a.round
    rdir = hist / f"round{rd:02d}"
    bo, errs = load_checked(ws / "build_result.json", "build")
    if not errs and bo.get("status") == "failed":
        errs = ["the builder reported status failed"]
    built_sha = (read_json(ws / "results.json", {}) or {}).get("params_sha")
    sheet = read_json(sd / "params.json")
    if not errs and built_sha and isinstance(sheet, dict) and built_sha != params_sha(sheet):
        errs = ["sim/workspace/results.json was built from another parameter sheet"]
    if errs or not (ws / "results.json").exists():
        return fail("the build did not finish: " + "; ".join(errs[:4] + ([] if (ws / "results.json").exists() else ["sim/workspace/results.json is missing"])))
    bo["_round"] = rd
    write_json(ws / "build_result.json", bo)
    for f in ("simulate.py", "results.json", "sim_report.md", "run_log.txt"):
        if (ws / f).exists():
            shutil.copy2(ws / f, rdir / f)
    if (ws / "figures").exists():
        shutil.copytree(ws / "figures", rdir / "figures", dirs_exist_ok=True)
    write_json(rdir / "build_result.json", bo)
    f = rdir / "critic.prompt.md"
    left = render_prompt(P, "sim_critic", f, round=rd, rounds=a.rounds, build_summary=clean(bo.get("summary"), 3000),
                         build_issues="\n".join(f"- {clean(x)}" for x in bo.get("issues", [])) or "- (none reported)",
                         evidence_sha=evidence_sha(P))
    return done({"ok": True, "status": bo["status"], "pts_overall": bo.get("pts_overall"), "mc_se": bo.get("mc_se"),
                 "n_sims": bo.get("n_sims"), "prompt": str(f), "result_file": str(rdir / "critique.json"), "unfilled": left})


def cmd_sim_critic_apply(a):
    P = Path(a.project)
    sd, ws, hist = sim_dirs(P)
    rd = a.round
    rdir = hist / f"round{rd:02d}"
    crit, errs = load_checked(rdir / "critique.json", "critic")
    if errs:
        return fail("the critique does not match the schema: " + "; ".join(errs[:5]))
    write_text(sd / f"critique_r{rd}.md",
               f"# Simulation critique — round {rd}\n\n*verdict: {crit['verdict']} · scores {json.dumps(crit['scores'])}*\n\n"
               + crit["critique_markdown"] + "\n\n## Change requests\n" + requests_text(crit, {"params", "code", "report"})
               + "\n\n## Open research questions\n" + "\n".join(f"- {x}" for x in crit.get("open_research_questions", [])))
    return done({"ok": True, "verdict": crit["verdict"], "scores": crit["scores"], "change_requests": len(crit["change_requests"]),
                 "another_round": crit["verdict"] != "accept" and rd < a.rounds})


def cmd_sim_publish(a):
    """The last built round -> sim/ (simulate.py, results.json, sim_report.md, run_log.txt, figures/) and sim/status.json;
    then the finalize prompt."""
    P = Path(a.project)
    sd, ws, hist = sim_dirs(P)
    built = sorted(int(d.name[5:]) for d in hist.glob("round[0-9][0-9]") if (d / "build_result.json").exists())
    if not built:
        return fail("no model round finished its build, so there is nothing to publish")
    final = built[-1]
    fdir = hist / f"round{final:02d}"
    # the sheet that round built: its own, else the latest earlier round's
    sheet = next((hist / f"round{r:02d}" / "params.json" for r in range(final, 0, -1) if (hist / f"round{r:02d}" / "params.json").exists()), None)
    if sheet is not None:
        shutil.copy2(sheet, sd / "params.json")
        shutil.copy2(sheet, ws / "params.json")
    for f in ("simulate.py", "results.json", "sim_report.md", "run_log.txt"):
        if (fdir / f).exists():
            shutil.copy2(fdir / f, sd / f)
            shutil.copy2(fdir / f, ws / f)
    shutil.copy2(SKILL / "scripts" / "bg.py", ws / "bg.py")  # for the finalize step's long commands
    if (fdir / "figures").exists():
        shutil.copytree(fdir / "figures", sd / "figures", dirs_exist_ok=True)
    crit = read_json(fdir / "critique.json", {}) or {}
    res = read_results(fdir)
    write_json(sd / "status.json", {"final_round": final, "verdict": crit.get("verdict"), "pts_overall": res.get("pts_overall"), "finished": now_iso()})
    params = read_json(sd / "params.json", {}) or {}
    S = size_of(a.size)
    f = ws / "finalize.prompt.md"
    fr = ws / "finalize_result.json"
    if fr.is_file() and (read_json(fr, {}) or {}).get("_published_round") != final:
        fr.unlink()
    write_json(sd / "published_round.json", {"round": final})
    left = render_prompt(P, "sim_finalize", f, n_sims=S["n_sims"], params_sha=params_sha(params))
    return done({"ok": True, "final_round": final, "verdict": crit.get("verdict"), "pts_overall": res.get("pts_overall"),
                 "mc_se": res.get("mc_se"), "prompt": str(f), "result_file": str(fr), "unfilled": left})


def cmd_sim_finalize_apply(a):
    P = Path(a.project)
    sd, ws, hist = sim_dirs(P)
    bo, errs = load_checked(ws / "finalize_result.json", "build")
    if errs:
        return fail("the finalize step did not finish: " + "; ".join(errs[:4]))
    bo["_published_round"] = (read_json(sd / "published_round.json", {}) or {}).get("round")
    write_json(ws / "finalize_result.json", bo)
    produced = []
    for f in ("flow.json", "parameter_metadata.json", "summary_stats.json", "finalize_log.txt", "finalize.py"):
        if (ws / f).exists():
            shutil.copy2(ws / f, sd / f)
            produced.append(f)
    if (ws / "figures").exists():
        shutil.copytree(ws / "figures", sd / "figures", dirs_exist_ok=True)
    res = read_results(sd)
    try:
        drift = abs(float(bo.get("pts_overall") or 0) - float(res.get("pts_overall") or 0))
    except (TypeError, ValueError):
        drift = None
    return done({"ok": True, "status": bo.get("status"), "files": produced, "drift": drift,
                 "figures": len(list((sd / "figures").glob("*.png"))), "issues": [clean(x, 300) for x in bo.get("issues") or []][:8]})


