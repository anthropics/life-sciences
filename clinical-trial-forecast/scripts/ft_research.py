# Copyright 2026 Anthropic PBC
# SPDX-License-Identifier: Apache-2.0
"""Trial cards and the research loop of the clinical-trial-forecast tools."""

import collections
import hashlib
import json
import re
import shutil
import subprocess
import sys
from pathlib import Path

from ft_base import *  # noqa: F401,F403


# ---- trial cards, on the trial-dossier checklist -------------------------------------------------------------------------------
def load_bank():
    """The checklist and its worker groups: reference/checklist.json holds prompts/trial_dossier_bank.yaml and
    prompts/trial_dossier_groups.yaml as JSON, because the standard library has no YAML reader. Edit both together."""
    c = read_json(SKILL / "reference" / "checklist.json")
    bank, groups = c["bank"], c["groups"]
    items = bank["items"] if isinstance(bank, dict) else bank
    if isinstance(groups, dict):
        groups = groups.get("groups") or []
    return items, groups


def applicable(item, t):
    """Which checklist items apply to a trial."""
    aw = str(item.get("applies_when", "always")).lower()
    status = t.get("status", "")
    phases = set(t.get("phase") or [])
    if aw in ("always", "has_dmc", ""):
        return True
    if aw == "completed_with_results":
        return status not in ("NOT_YET_RECRUITING", "RECRUITING", "WITHDRAWN")
    if aw == "ongoing":
        return status in ACTIVE
    if aw.startswith("phase2"):
        return bool(phases & {"PHASE2", "PHASE3", "PHASE4"})
    return True


def slim(tj):
    t = tj["trial"]
    keys = ("nct_id", "title", "brief_title", "acronym", "status", "phase", "study_type", "conditions", "interventions", "sponsor",
            "collaborators", "enrollment", "start_date", "primary_completion_date", "completion_date", "brief_summary",
            "detailed_description", "eligibility_criteria", "minimum_age", "maximum_age", "sex", "primary_outcomes",
            "secondary_outcomes", "other_outcomes", "url", "has_results")
    return {"trial": {k: t.get(k) for k in keys}, "n_locations": len(t.get("locations") or []),
            "countries": sorted({l.get("country") for l in (t.get("locations") or []) if isinstance(l, dict) and l.get("country")}),
            "classification": tj.get("classification"), "tiers": tj.get("tiers")}


def fmt_items(items):
    return "\n".join(f"- id: {i['id']}  [priority {i.get('priority')}] \n  question: {i['question']}\n  where to look: {i.get('inputs')}\n"
                     f"  output spec: {i.get('output')}\n  why: {i.get('why')}" for i in items)


def cmd_card_prompts(a):
    P = Path(a.project)
    n = a.nct.strip().upper()
    tj = read_json(P / "trials" / f"{n}.json")
    if not isinstance(tj, dict):
        return fail(f"trials/{n}.json is missing: the index step must run first")
    items, groups = load_bank()
    by_bank = {}
    for it in items:
        by_bank.setdefault(it["group"], []).append(it)
    tslim = json.dumps(slim(tj), ensure_ascii=False)[:50000]
    pj = project(P)
    out = []
    for g in groups:
        gitems = [it for bg in g["bank_groups"] for it in by_bank.get(bg, []) if applicable(it, tj["trial"])]
        gitems.sort(key=lambda x: x.get("priority", 3))
        if not gitems:
            continue
        f = P / "analyses" / n / "prompts" / f"group-{g['id']}.md"
        left = render_prompt(P, "trial_checklist_worker", f, group_title=g["title"], fetched_at=tj.get("fetched_at"), items=fmt_items(gitems),
                             trial_json=tslim, reference_population=(pj.get("indication") or {}).get("reference_population", ""))
        out.append({"group": g["id"], "title": g["title"], "prompt": str(f), "items": [i["id"] for i in gitems], "unfilled": left})
    return done({"ok": bool(out), "groups": out, "problem": "" if out else "no checklist item applies to this trial"})


def render_group_md(g, answers):
    L = [f"## {g['title']}", ""]
    for x in answers:
        L += [f"### {x['item']['id']} — {x['item']['question'][:160]}", f"*confidence: {x['confidence']}*", "", x["answer_markdown"], ""]
        ans = x.get("answer") or {}
        for k, v in ans.items():
            if isinstance(v, list) and v and isinstance(v[0], dict):
                cols = list({c for row in v if isinstance(row, dict) for c in row.keys()})[:10]
                L += [f"**{k}**", md_table(v, cols, max_len=60), ""]
        if x.get("gaps"):
            L.append("Gaps: " + "; ".join(map(str, x["gaps"])))
        L.append("")
    return "\n".join(L)


def cmd_card_apply(a):
    """analyses/<NCT>/groups.json ({"groups": [{"group": id, ...the trial_checklist schema}]}) -> items/<item>.json and groups/<group>.md,
    then the card writer's prompt."""
    P = Path(a.project)
    n = a.nct.strip().upper()
    d = P / "analyses" / n
    got = read_json(d / "groups.json")
    if not isinstance(got, dict) or not isinstance(got.get("groups"), list):
        return fail(f"analyses/{n}/groups.json is missing or has no groups list")
    items, groups = load_bank()
    by_item = {i["id"]: i for i in items}
    gmap = {g["id"]: g for g in groups}
    done_groups, bad, reqs = [], [], []
    for gr in got["groups"]:
        if not isinstance(gr, dict):
            continue
        gid = gr.get("group")
        errs = schema_errors(gr, schemas()["trial_checklist"])
        if gid not in gmap or errs:
            bad.append(f"{gid}: " + ("unknown group" if gid not in gmap else "; ".join(errs[:3])))
            continue
        answers = []
        for x in gr.get("items") or []:
            it = by_item.get(x.get("id"))
            if not it:
                continue
            rec = {"nct_id": n, "group": gid, "item": it, "answer": x.get("answer"), "answer_markdown": x.get("answer_markdown", ""),
                   "confidence": x.get("confidence"), "sources": x.get("sources", []), "gaps": x.get("gaps", []), "updated": now_iso()}
            write_json(d / "items" / f"{it['id']}.json", rec)
            answers.append(rec)
        adj = gr.get("adjacent_findings") or []
        write_text(d / "groups" / f"{gid}.md", render_group_md(gmap[gid], answers)
                   + ("\n\nAdjacent findings:\n" + "\n".join(f"- {x}" for x in adj) if adj else ""))
        reqs += gr.get("paper_requests") or []
        done_groups.append(gid)
    if reqs:
        write_json(d / "paper_requests.json", reqs)
    if not done_groups:
        return fail("no group answer could be used: " + "; ".join(bad[:5]))
    f = d / "prompts" / "condense.md"
    left = render_prompt(P, "trial_checklist_condense", f, nct=n, group_titles="; ".join(g["title"] for g in groups))
    return done({"ok": True, "groups": done_groups, "rejected": bad, "condense_prompt": str(f), "unfilled": left})


def cmd_card_finish(a):
    """analyses/<NCT>/condense.json -> analysis.md, card.md, analysis.json; then the number checker's prompt (the card's
    key facts that hold a number, the card, and every source its items cited)."""
    P = Path(a.project)
    n = a.nct.strip().upper()
    d = P / "analyses" / n
    o, errs = load_checked(d / "condense.json", "condense_trial")
    if errs:
        return fail("; ".join(errs[:5]))
    write_text(d / "analysis.md", o["analysis_markdown"] + "\n")
    write_text(d / "card.md", o["card_markdown"].strip() + f"\n\n<!-- {n} · card {now_iso()} -->\n")
    gids = sorted(f.stem for f in (d / "groups").glob("*.md"))
    write_json(d / "analysis.json", {"nct_id": n, "groups": gids, "key_facts": o["key_facts"], "gaps": o["gaps"],
                                     "items": sorted(f.stem for f in (d / "items").glob("*.json")), "source": "trial_checklist"})
    srcs, seen = [], set()
    for f in sorted((d / "items").glob("*.json")):
        for x in (read_json(f, {}) or {}).get("sources", []) or []:
            if x and x not in seen:
                seen.add(x)
                srcs.append(x)
    nums = [k for k in o["key_facts"] if any(ch.isdigit() for ch in str(k))][:30]
    if not nums:
        return done({"ok": True, "check_prompt": "", "note": "no key fact holds a number, so there is nothing to check"})
    card = (d / "card.md").read_text(encoding="utf-8")
    f = d / "prompts" / "check.md"
    left = render_prompt(P, "card_checker", f,
                         key_numbers="\n".join(f"- {k}" for k in nums) + "\n\nCARD (for context of arms/timepoints):\n" + card[:3000],
                         sources="\n".join(f"- {x}" for x in srcs[:40]) + f"\n- trials/{n}.json (registry record incl. posted results)"
                         f"\n- analyses/{n}/items/*.json (the item answers with their own citations)")
    return done({"ok": True, "check_prompt": str(f), "key_numbers": len(nums), "unfilled": left})


def check_lines(chk):
    lines = []
    for c in chk.get("checks") or []:
        if not isinstance(c, dict):
            continue
        if c.get("verdict") == "corrected":
            lines.append(f"- CORRECTED: {c.get('key_number', '')} → **{c.get('corrected_value', '')}** ({c.get('what_went_wrong', '')}; source: {cut(c.get('source_quote'), 160)})")
        elif c.get("verdict") == "unverifiable":
            lines.append(f"- unverifiable: {c.get('key_number', '')} ({cut(c.get('what_went_wrong'), 120)})")
        else:
            lines.append(f"- confirmed: {c.get('key_number', '')}")
    return lines


def cmd_card_check_apply(a):
    P = Path(a.project)
    n = a.nct.strip().upper()
    d = P / "analyses" / n
    if not (d / "check.json").exists():
        return done({"ok": True, "note": f"analyses/{n}/check.json does not exist (no number check for this card)"})
    chk, errs = load_checked(d / "check.json", "check")
    if errs:
        # not a failure of this command: the next trials' checks still get applied, and the check reruns next launch
        return done({"ok": True, "note": f"analyses/{n}/check.json was not usable: " + "; ".join(errs[:3])})
    aj = read_json(d / "analysis.json", {}) or {}
    if aj.get("check") == chk:
        return done({"ok": True, "note": "already applied"})
    aj["check"] = chk
    write_json(d / "analysis.json", aj)
    checks = [c for c in chk["checks"] if isinstance(c, dict)]
    lines = ["\n\n**Number check (independent re-derivation from sources)**"]
    lines += [x for x in check_lines(chk) if not x.startswith("- confirmed")]
    if len(lines) == 1:
        lines.append(f"- all {len(checks)} checked numbers confirmed")
    lines.append(f"- Checker summary: {chk['summary']}")
    txt = (d / "card.md").read_text(encoding="utf-8")
    marker = f"\n\n<!-- {n} · card"
    i = txt.rfind(marker)
    write_text(d / "card.md", (txt[:i] if i > 0 else txt) + "\n".join(lines) + (txt[i:] if i > 0 else "\n"))
    tally = collections.Counter(c.get("verdict") for c in checks)
    return done({"ok": True, "verdicts": dict(tally)})


def cmd_digest_prompt(a):
    P = Path(a.project)
    cards = list((P / "analyses").glob("*/card.md"))
    n_cards = len(cards)
    dg = P / "knowledge" / "trials_digest.md"
    if dg.exists() and any(c.stat().st_mtime > dg.stat().st_mtime for c in cards):
        dg.rename(dg.with_name("trials_digest.superseded.md"))  # a card written after the digest: write the digest again
    f = P / "knowledge" / "prompts" / "digest.md"
    left = render_prompt(P, "condense_trials", f, n_cards=n_cards)
    return done({"ok": n_cards > 0, "prompt": str(f), "cards": n_cards, "unfilled": left,
                 "problem": "" if n_cards else "no trial card exists"})



# ---- research loop, with a follow-up wave chosen by a triage step --------------------------------------------------------------
def iter_dir(P, it):
    return P / "research" / f"iter{it:02d}"


def iter_dirs(P):
    """(iteration number, folder) for every research/iterNN folder, in order; other folders are ignored."""
    out = []
    for d in (P / "research").glob("iter*"):
        m = re.fullmatch(r"iter(\d+)", d.name)
        if m and d.is_dir():
            out.append((int(m.group(1)), d))
    return sorted(out)


# knowledge/<section>.md shares its folder with these names
RESERVED = {"trials_digest", "learning_log", "open_questions", "blocked_items", "topics", "history", "prompts"}


def section_id(raw):
    t = slug(raw, "general")
    return f"section_{t}" if t in RESERVED else t


def wave_file(P, it):
    """research/iterNN/questions.json: every question this iteration runs, wave 1 (planner) then wave 2 (follow-ups)."""
    return iter_dir(P, it) / "questions.json"


def answer_file(P, it, i):
    """The i-th question's (0-based) answer file in iteration it: {"question", "answer", "check"}."""
    return iter_dir(P, it) / f"q{i + 1:02d}.json"


def read_answer(P, it, i):
    """The i-th answer file as a dict ({} if there is none). An answer in which a list or an object is of the wrong type
    counts as missing, so that one malformed answer cannot stop a command."""
    got = read_json(answer_file(P, it, i))
    if not isinstance(got, dict):
        return {}
    if isinstance(got.get("answer"), dict) and wrong_shape(got["answer"], "answer"):
        return {k: v for k, v in got.items() if k != "answer"}
    return got


def iteration_answers(P, it):
    out = []
    for i, _ in enumerate((read_json(wave_file(P, it), {}) or {}).get("questions") or []):
        got = read_answer(P, it, i)
        if isinstance(got.get("answer"), dict):
            out.append(got["answer"])
    return out


def cmd_planner_prompt(a):
    P = Path(a.project)
    f = iter_dir(P, a.iter) / "planner.prompt.md"
    left = render_prompt(P, "research_planner", f, iteration=a.iter, k=a.k)
    return done({"ok": True, "prompt": str(f), "unfilled": left})


def researcher_prompt(P, it, i, q):
    qid = f"q{it:02d}_{i + 1:02d}"
    f = iter_dir(P, it) / f"q{i + 1:02d}.prompt.md"
    render_prompt(P, "researcher", f, qid=qid, section=q["section"], question=q["question"], why=q["why_it_matters"],
                  what_would_answer=q["what_would_answer_it"])
    return str(f)


def cmd_questions(a):
    """After the planner (--wave 1) or the follow-up triage (--wave 2): orders its questions (those
    that could flip the recommendation first, then by priority), keeps at most k, applies the stop rule on wave 1, and
    writes each question's researcher prompt. The stop rule: the planner recommends stopping, at least
    min_iterations iterations ran before this one, and the previous iteration had at least 60% high-confidence answers
    and at most two follow-up questions per answer. On a stop, this iteration runs no questions."""
    P = Path(a.project)
    d = iter_dir(P, a.iter)
    src = d / ("plan.json" if a.wave == 1 else "followups.json")
    plan, errs = load_checked(src, "plan")
    if errs:
        return fail(f"{rel(P, src)}: " + "; ".join(errs[:5]))
    qs = sorted([q for q in plan["questions"] if isinstance(q, dict) and str(q.get("question") or "").strip()],
                key=lambda q: (not truthy(q.get("could_flip_recommendation")), 3 if to_int(q.get("priority")) is None else to_int(q.get("priority"))))[: a.k]
    wf = read_json(wave_file(P, a.iter), {}) or {}
    have = wf.get("questions") or []
    if a.wave == 1:
        have = [q for q in have if q.get("wave") == 1]
    else:
        have = [q for q in have if q.get("wave") == 1] + [q for q in have if q.get("wave") == 2]
    stop = False
    if a.wave == 1:
        if plan.get("blocked_items"):
            bi = P / "knowledge" / "blocked_items.json"
            prev = read_json(bi, {"items": []}) or {"items": []}
            for b in plan["blocked_items"]:
                if b not in prev["items"]:
                    prev["items"].append(b)
            write_json(bi, prev)
        earlier = [n for n, _ in iter_dirs(P) if n < a.iter]
        if plan.get("stop_recommended") and len(earlier) >= a.min_iterations:
            prev = iteration_answers(P, a.iter - 1)
            hi = sum(1 for x in prev if x.get("confidence") == "high")
            fups = sum(len(x.get("follow_up_questions") or []) for x in prev)
            stop = bool(prev) and hi >= 0.6 * len(prev) and fups <= 2 * len(prev)
    new = [] if stop else [{"wave": a.wave, "section": section_id(q.get("section") or q.get("topic")), "question": str(q.get("question")),
                            "why_it_matters": str(q.get("why_it_matters") or ""), "what_would_answer_it": str(q.get("what_would_answer_it") or ""),
                            "could_flip_recommendation": truthy(q.get("could_flip_recommendation"))} for q in qs]
    if a.wave == 2 or not have:
        have = [q for q in have if q.get("wave") == 1] + new if a.wave == 2 else new
        write_json(wave_file(P, a.iter), {"questions": have})
    out = []
    start = len([q for q in have if q.get("wave") == 1]) if a.wave == 2 else 0
    for i, q in enumerate(have):
        if q.get("wave") != a.wave:
            continue
        out.append({"qid": f"q{i + 1:02d}", "section": q["section"], "prompt": researcher_prompt(P, a.iter, i, q),
                    "answer_file": str(answer_file(P, a.iter, i))})
    return done({"ok": True, "stop": stop, "questions": out, "stop_recommended": bool(plan.get("stop_recommended")), "first": start})


def cmd_checks(a):
    """After the researchers: the number checker's prompt for each answer with key numbers that has no check yet."""
    P = Path(a.project)
    out = []
    for i, q in enumerate((read_json(wave_file(P, a.iter), {}) or {}).get("questions") or []):
        if a.wave and q.get("wave") != a.wave:
            continue
        f = answer_file(P, a.iter, i)
        got = read_answer(P, a.iter, i)
        ans = got.get("answer")
        if not isinstance(ans, dict) or not ans.get("key_numbers"):
            continue
        # the prompt lives outside research/, which the checker may not open
        pf = P / "checks" / "research" / f"iter{a.iter:02d}-{f.stem}.md"
        render_prompt(P, "card_checker", pf, key_numbers="\n".join(f"- {k}" for k in ans["key_numbers"][:40]),
                      sources="\n".join(f"- {x}" for x in (ans.get("sources") or [])[:40]))
        out.append({"qid": f.stem, "prompt": str(pf), "check_file": str(f.with_suffix(".check.json"))})
    return done({"ok": True, "checks": out})


def cmd_followups(a):
    """The follow-up triage prompt (prompts/followup_triage.md) from this iteration's answers' follow-up questions."""
    P = Path(a.project)
    rows = []
    for i, q in enumerate((read_json(wave_file(P, a.iter), {}) or {}).get("questions") or []):
        got = read_answer(P, a.iter, i)
        for fq in (got.get("answer") or {}).get("follow_up_questions") or []:
            rows.append(f"- {q['section']} · q{a.iter:02d}_{i + 1:02d} · {fq}")
    if not rows or a.k <= 0:
        return done({"ok": True, "prompt": "", "followups": len(rows)})
    f = iter_dir(P, a.iter) / "followups.prompt.md"
    left = render_prompt(P, "followup_triage", f, iteration=a.iter, iteration_dir=f"iter{a.iter:02d}", k=a.k, followups="\n".join(rows[:80]))
    return done({"ok": True, "prompt": str(f), "followups": len(rows), "unfilled": left})


def answer_md(qid, q, a, check):
    """One research answer as markdown, with its independent number check under it when there is one."""
    txt = (f"# {qid} · {q.get('section', '')}\n\n**Question.** {q.get('question', '')}\n\n**Bottom line ({a.get('confidence')}).** "
           f"{a.get('bottom_line', '')}\n\n{a.get('answer_markdown', '')}\n\n## Key numbers\n" + "\n".join(f"- {k}" for k in a.get("key_numbers") or [])
           + "\n\n## Gaps\n" + "\n".join(f"- {g}" for g in a.get("gaps") or []) + "\n\n## Sources\n" + "\n".join(f"- {s}" for s in a.get("sources") or []) + "\n")
    if isinstance(check, dict):
        txt += "\n".join(["\n## Number check (independent re-derivation from sources)"] + check_lines(check)
                         + [f"- Checker summary: {check.get('summary', '')}"]) + "\n"
    return txt


def cmd_research_apply(a):
    """Answers and their checks -> research/iterNN/qNN.md; the integrator prompt for each section touched."""
    P = Path(a.project)
    it = a.iter
    d = iter_dir(P, it)
    pj = project(P)
    secs = {s["id"]: s for s in sections_of(pj)}
    qs = (read_json(wave_file(P, it), {}) or {}).get("questions") or []
    touched = collections.OrderedDict()
    answers, fups = [], []
    for i, q in enumerate(qs):
        f = answer_file(P, it, i)
        got = read_answer(P, it, i)
        if not isinstance(got.get("answer"), dict):
            continue
        chk = read_json(f.with_suffix(".check.json"))
        if isinstance(chk, dict) and isinstance(chk.get("checks"), list):
            got["check"] = chk
            write_json(f, got)
        md = f.with_suffix(".md")
        write_text(md, answer_md(f"q{it:02d}_{i + 1:02d}", q, got["answer"], got.get("check")))
        answers.append(got["answer"])
        fups += [{"section": q["section"], "question": x, "from": f"q{it:02d}_{i + 1:02d}"} for x in got["answer"].get("follow_up_questions") or []]
        touched.setdefault(q["section"], []).append(rel(P, md))
    out = []
    for sec, files in touched.items():
        s = secs.get(sec) or {"id": sec, "title": sec.replace("_", " "), "question": ""}
        pf = d / f"integrate__{sec}.prompt.md"
        render_prompt(P, "integrator", pf, section=sec, section_title=s.get("title", ""), section_question=s.get("question", ""),
                      iteration=f"iter{it:02d}", answer_files=", ".join(files))
        out.append({"section": sec, "prompt": str(pf), "result_file": str(d / f"integrate__{sec}.json")})
    write_json(d / "followups_raised.json", fups)
    hi = sum(1 for x in answers if x.get("confidence") == "high")
    return done({"ok": True, "iteration": it, "answers": len(answers), "questions": len(qs), "high_confidence": hi,
                 "follow_ups": len(fups), "integrations": out})


def cmd_integrate_apply(a):
    """The integrators' answers -> knowledge/<section>.md (each earlier version kept under knowledge/history/), then the
    open questions and the learning log."""
    P = Path(a.project)
    it = a.iter
    d = iter_dir(P, it)
    K = P / "knowledge"
    done_s, bad, opens = [], [], []
    for f in sorted(x for x in d.glob("integrate__*.json") if not x.name.endswith(".rejected.json")):
        sec = f.stem.split("__", 1)[1]
        o, errs = load_checked(f, "integrate")
        if errs:
            bad.append(f"{sec}: " + "; ".join(errs[:2]))
            continue
        kf = K / f"{sec}.md"
        stamp = f"<!-- integrated iteration {it} -->"
        if kf.is_file() and stamp in kf.read_text(encoding="utf-8"):
            done_s.append(sec)
            opens += [{"section": sec, "question": x, "status": "open", "iteration": it} for x in o["open_questions"]]
            continue
        if kf.is_file():
            write_text(K / "history" / f"{sec}.iter{it - 1:02d}.md", kf.read_text(encoding="utf-8"))
        write_text(kf, o["synthesis_markdown"] + f"\n\n{stamp}\n")
        opens += [{"section": sec, "question": x, "status": "open", "iteration": it} for x in o["open_questions"]]
        done_s.append(sec)
    fups = read_json(d / "followups_raised.json", []) or []
    write_json(K / "open_questions.json", {"iteration": it, "open": opens + [{**x, "status": "open"} for x in fups]})
    learning_log(P)
    return done({"ok": True, "sections": done_s, "rejected": bad})


def learning_log(P):
    """knowledge/learning_log.md: per iteration, the planner's reasoning and what each question learned."""
    L = ["# Learning log: research loop", ""]
    for it, d in iter_dirs(P):
        plan = (read_json(d / "plan.json", {}) or {})
        L += [f"## Iteration {it}", f"**Planner's reasoning.** {plan.get('rationale', '')}", f"*Stop recommended: {plan.get('stop_recommended')}*", ""]
        for i, q in enumerate((read_json(d / "questions.json", {}) or {}).get("questions") or []):
            got = read_answer(P, it, i)
            ans = got.get("answer") or {}
            flip = " (could flip the recommendation)" if q.get("could_flip_recommendation") else ""
            wave = " (follow-up)" if q.get("wave") == 2 else ""
            L.append(f"- **q{it:02d}_{i + 1:02d}**{flip}{wave} [{q.get('section')}]: {str(q.get('question', ''))[:300]}")
            if ans:
                L.append(f"  - Learned ({ans.get('confidence')}): {ans.get('bottom_line', '')}")
                kn = ans.get("key_numbers") or []
                if kn:
                    L.append("  - Key numbers: " + " · ".join(str(k)[:140] for k in kn[:4]))
            else:
                L.append("  - (no answer on disk)")
        L.append("")
    oq = (read_json(P / "knowledge" / "open_questions.json", {}) or {}).get("open", [])
    L += [f"## Still open ({len(oq)} questions)"] + [f"- [{o.get('section')}] {str(o.get('question', ''))[:220]}" for o in oq[:40]]
    write_text(P / "knowledge" / "learning_log.md", "\n".join(L))


