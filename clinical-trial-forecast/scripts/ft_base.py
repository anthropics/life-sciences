#!/usr/bin/env python3
# Copyright 2026 Anthropic PBC
# SPDX-License-Identifier: Apache-2.0
"""Shared parts of the clinical-trial-forecast tools: constants, file helpers, registry records, the project's prompt
variables and the prompt renderer, and a small JSON-schema checker. scripts/ft_tools.py is the command line."""

import argparse
import collections
import csv
import datetime as dt
import hashlib
import io
import json
import math
import re
import shutil
import subprocess
import sys
import urllib.parse
from pathlib import Path

NCT = re.compile(r"NCT[0-9]{8}")  # use with fullmatch: "NCT01234567\n" is not an id
API = "https://clinicaltrials.gov/api/v2/studies"
ACTIVE = ["NOT_YET_RECRUITING", "RECRUITING", "ENROLLING_BY_INVITATION", "ACTIVE_NOT_RECRUITING"]
PATHWAY = "target_pathway"
SKILL = Path(__file__).resolve().parent.parent

# The trial universe: the default tier rules.
UNIVERSE = {
    "study_type": "INTERVENTIONAL",
    "intervention_types": ["DRUG", "BIOLOGICAL", "GENETIC", "COMBINATION_PRODUCT"],
    "A": {"label": "enrollment competitors", "statuses": ACTIVE, "also_completed_since": "2023-01-01"},
    "B": {"label": "endpoint and placebo evidence base", "statuses": ["COMPLETED", "TERMINATED"],
          "start_since": "2015-01-01", "phases": ["PHASE2", "PHASE3"], "require_key_endpoint": True},
    "C": {"label": "strategic analogs"},
}
# Run sizes. "test" is an end-to-end check of a setup (every phase, small amounts). "full" is the most thorough (16
# research iterations of k 8 with no follow-up wave, 3 model rounds, 20,000 simulations, 40 claims, 16 helpers at
# once; the other sizes run 8 at once). Only cs_workflow.py's copy of these sizes reads at_once. 0 means no
# limit.
SIZES = {
    "test": {"max_indications": 1, "per_indication": 15, "per_drug": 5, "max_drugs": 2, "cards": 3,
             "iterations": 1, "questions": 2, "followups": 0, "rounds": 1, "n_sims": 5000, "claims": 5, "batch": 15, "at_once": 8},
    "small": {"max_indications": 2, "per_indication": 40, "per_drug": 10, "max_drugs": 6, "cards": 8,
              "iterations": 2, "questions": 4, "followups": 0, "rounds": 2, "n_sims": 20000, "claims": 12, "batch": 15, "at_once": 8},
    "medium": {"max_indications": 3, "per_indication": 150, "per_drug": 25, "max_drugs": 12, "cards": 15,
               "iterations": 6, "questions": 8, "followups": 0, "rounds": 3, "n_sims": 20000, "claims": 20, "batch": 15, "at_once": 8},
    "full": {"max_indications": 0, "per_indication": 600, "per_drug": 0, "max_drugs": 0, "cards": 25,
             "iterations": 16, "questions": 8, "followups": 0, "rounds": 3, "n_sims": 20000, "claims": 40, "batch": 20, "at_once": 16},
}
TRIAL_KEYS = ["nct_id", "title", "brief_title", "acronym", "status", "phase", "study_type", "conditions", "interventions",
              "sponsor", "collaborators", "enrollment", "start_date", "primary_completion_date", "completion_date",
              "brief_summary", "detailed_description", "eligibility_criteria", "minimum_age", "maximum_age", "sex",
              "healthy_volunteers", "primary_outcomes", "secondary_outcomes", "other_outcomes", "locations", "has_results",
              "url"]
LISTS = {"phase", "conditions", "interventions", "collaborators", "primary_outcomes", "secondary_outcomes",
         "other_outcomes", "locations"}


# ---- small helpers ---------------------------------------------------------------------------------------------------
def now_iso():
    return dt.datetime.now(dt.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def read_json(path, default=None):
    try:
        return json.loads(Path(path).read_text(encoding="utf-8-sig"))  # utf-8-sig: a byte-order mark, if any, is dropped
    except (OSError, ValueError, RecursionError):
        return default


def finite_only(v):
    """v with every NaN and Infinity replaced by None. Python's json reads and writes them, but they are not JSON."""
    if isinstance(v, float) and not math.isfinite(v):
        return None
    if isinstance(v, dict):
        return {k: finite_only(x) for k, x in v.items()}
    if isinstance(v, list):
        return [finite_only(x) for x in v]
    return v


def read_results(folder):
    """folder/results.json, which the model step wrote, as a dict of finite values ({} if it is missing or not an object)."""
    r = read_json(Path(folder) / "results.json")
    return finite_only(r) if isinstance(r, dict) else {}


def copy_pngs(src, dst, replace=False):
    """Copies src's PNG files to dst: regular files only, no links, no other kind of file. Keeps a file dst already has
    unless replace is set."""
    Path(dst).mkdir(parents=True, exist_ok=True)
    for f in sorted(Path(src).glob("*.png")):
        if f.is_file() and not f.is_symlink() and (replace or not (Path(dst) / f.name).exists()):
            shutil.copy2(f, Path(dst) / f.name)


def write_json(path, data):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(data, indent=1, ensure_ascii=False), encoding="utf-8")


def write_text(path, text):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text, encoding="utf-8")


def cell(v):
    if isinstance(v, (list, tuple, set)):
        return "; ".join(str(x) for x in v)
    if isinstance(v, dict):
        return json.dumps(v, ensure_ascii=False)
    return "" if v is None else str(v)


def csv_safe(v):
    """Neutralises spreadsheet formula injection: third-party text (titles, sponsors, sites) goes into the CSVs."""
    s = cell(v)
    return "'" + s if s and s[0] in "=+-@\t\r" else s


def to_int(v):
    """An integer from a model-written value ("180", 180.0, None); None when it is not a number."""
    if isinstance(v, bool):
        return None
    try:
        if isinstance(v, (int, float)):
            return int(v)
        m = re.match(r"^\s*(\d+)", str(v or ""))
        return int(m.group(1)) if m else None
    except (ValueError, OverflowError):  # NaN, Infinity, or more digits than int() takes
        return None


def write_csv(path, rows, cols):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    buf = io.StringIO()
    w = csv.writer(buf)
    w.writerow(cols)
    for r in rows:
        w.writerow([csv_safe(r.get(c)) for c in cols])
    path.write_text(buf.getvalue(), encoding="utf-8")


def md_table(rows, cols, max_len=80):
    if not rows:
        return "(none)"

    def c(v):
        s = cell(v).replace("|", "/").replace("\n", " ")
        return s if len(s) <= max_len else s[: max_len - 1] + "…"

    out = ["| " + " | ".join(cols) + " |", "|" + "---|" * len(cols)]
    out += ["| " + " | ".join(c(r.get(k)) for k in cols) + " |" for r in rows]
    return "\n".join(out)


def slug(s, default="x"):
    s = re.sub(r"[^a-z0-9]+", "_", str(s or "").lower()).strip("_")
    return s[:60] or default


def sha(obj):
    return hashlib.sha256(json.dumps(obj, sort_keys=True).encode()).hexdigest()


def date(s):
    """YYYY-MM-DD, YYYY-MM or YYYY as a date (the first of the month or year when the day or month is missing)."""
    m = re.match(r"^\s*(\d{4})(?:-(\d{2}))?(?:-(\d{2}))?", str(s or ""))
    if not m:
        return None
    try:
        return dt.date(int(m.group(1)), int(m.group(2) or 1), int(m.group(3) or 1))
    except ValueError:
        return None


def truthy(v):
    return v is True or str(v).strip().lower() == "true"


def unique_id(base, seen):
    """base, or base_2, base_3, ... when base is taken; adds the result to seen."""
    out, k = base, 2
    while out in seen:
        out, k = f"{base}_{k}", k + 1
    seen.add(out)
    return out


def done(result):
    result.setdefault("problem", "")
    print(json.dumps(result, ensure_ascii=False))
    return 0


def fail(problem):
    print(json.dumps({"ok": False, "problem": problem}, ensure_ascii=False))
    return 1


def size_of(name):
    if name not in SIZES:
        raise SystemExit(f"unknown size {name!r}; use one of {', '.join(SIZES)}")
    return SIZES[name]


# ---- registry records ------------------------------------------------------------------------------------------------
ROMAN = {"I": "1", "II": "2", "III": "3", "IV": "4"}


def norm_phases(v):
    """Registry phases in the API's form: "Phase 1/Phase 2", "Phase II", "PHASE2" and ["Phase 2", "Phase 3"] all become
    ["PHASE1", "PHASE2"]-style lists; "Early Phase 1" becomes EARLY_PHASE1; "N/A" and unknown text are dropped."""
    items = v if isinstance(v, (list, tuple)) else [v]
    out = []
    for x in items:
        s = str(x or "").upper()
        if re.search(r"EARLY[ _]?PHASE[ _]?(1|I)[AB]?\b", s):
            out.append("EARLY_PHASE1")
            s = re.sub(r"EARLY[ _]?PHASE[ _]?(1|I)[AB]?\b", " ", s)
        for m in re.finditer(r"PHASE[ _]?(IV|III|II|I|[1-4])[AB]?\b|(?<=[/,&])\s*(?:PHASE[ _]?)?(IV|III|II|I|[1-4])[AB]?\b", s):
            num = m.group(1) or m.group(2)
            p = "PHASE" + ROMAN.get(num, num)
            if p not in out:
                out.append(p)
    return out


def normalise_trial(t):
    out = dict(t)
    for k in TRIAL_KEYS:
        if k not in out or out[k] is None:
            out[k] = [] if k in LISTS else None
    out["phase"] = norm_phases(out["phase"])
    out["enrollment"] = to_int(out.get("enrollment"))
    if isinstance(out.get("status"), str):
        out["status"] = out["status"].strip().upper().replace(",", "").replace(" ", "_")
    if not out.get("url") and out.get("nct_id"):
        out["url"] = f"https://clinicaltrials.gov/study/{out['nct_id']}"
    return out


def from_api_v2(study):
    """One ClinicalTrials.gov API v2 study (protocolSection ...) in the record shape."""
    ps = study.get("protocolSection") or {}
    ident, stat, design = ps.get("identificationModule") or {}, ps.get("statusModule") or {}, ps.get("designModule") or {}
    spons, desc, elig = ps.get("sponsorCollaboratorsModule") or {}, ps.get("descriptionModule") or {}, ps.get("eligibilityModule") or {}
    outc, cl = ps.get("outcomesModule") or {}, ps.get("contactsLocationsModule") or {}

    def d(x):
        return (x or {}).get("date")

    def outcomes(xs):
        return [{"measure": o.get("measure"), "time_frame": o.get("timeFrame")} for o in xs or []]

    n = ident.get("nctId")
    return normalise_trial({
        "nct_id": n, "title": ident.get("officialTitle"), "brief_title": ident.get("briefTitle"),
        "acronym": ident.get("acronym"), "status": stat.get("overallStatus"), "phase": design.get("phases") or [],
        "study_type": design.get("studyType"), "conditions": (ps.get("conditionsModule") or {}).get("conditions") or [],
        "interventions": [i.get("name") for i in (ps.get("armsInterventionsModule") or {}).get("interventions") or []],
        "sponsor": (spons.get("leadSponsor") or {}).get("name"),
        "collaborators": [c.get("name") for c in spons.get("collaborators") or []],
        "enrollment": (design.get("enrollmentInfo") or {}).get("count"),
        "start_date": d(stat.get("startDateStruct")), "primary_completion_date": d(stat.get("primaryCompletionDateStruct")),
        "completion_date": d(stat.get("completionDateStruct")), "brief_summary": desc.get("briefSummary"),
        "detailed_description": desc.get("detailedDescription"), "eligibility_criteria": elig.get("eligibilityCriteria"),
        "minimum_age": elig.get("minimumAge"), "maximum_age": elig.get("maximumAge"), "sex": elig.get("sex"),
        "healthy_volunteers": elig.get("healthyVolunteers"), "primary_outcomes": outcomes(outc.get("primaryOutcomes")),
        "secondary_outcomes": outcomes(outc.get("secondaryOutcomes")), "other_outcomes": outcomes(outc.get("otherOutcomes")),
        "locations": [{k: loc.get(k) for k in ("facility", "city", "state", "country", "status")} for loc in cl.get("locations") or []],
        "has_results": bool(study.get("hasResults")), "url": f"https://clinicaltrials.gov/study/{n}" if n else None,
    })


def load_record(P, n):
    """raw/ctgov/<NCT>.json as {nct_id, found, source, fetched_at, trial}. Also accepts a bare API v2 study or a bare
    trial object. None when the file is missing or says found false; the index also treats a record whose nct_id is not
    the file's as missing."""
    rec = read_json(P / "raw" / "ctgov" / f"{n}.json")
    if not isinstance(rec, dict) or rec.get("found") is False:
        return None
    if "protocolSection" in rec:
        res = rec.get("resultsSection")
        rec = {"nct_id": n, "found": True, "source": "api_v2", "trial": from_api_v2(rec)}
        if res:
            rec["trial"]["results"] = res  # the posted results, which the prompts read as part of the record
    elif "trial" not in rec and rec.get("nct_id"):
        rec = {"nct_id": n, "found": True, "source": rec.get("source", ""), "trial": rec}
    if not rec.get("found", True) or not isinstance(rec.get("trial"), dict):
        return None
    t = rec["trial"]
    if "protocolSection" in t:
        res = t.get("resultsSection")
        t = from_api_v2(t)
        if res:
            t["results"] = res
    rec["trial"] = normalise_trial(t)
    rec["trial"]["nct_id"] = rec["trial"].get("nct_id") or n
    return rec



# ---- the project (project.json, written by the scope step) and the prompts' variables ------------------------------------
REPORT_KIND = "probability-of-technical-success forecast"

# What a prompt's {{modeling_goals}} carries: these goals when project.json gives none, and the standing rules always.
DEFAULT_MODELING_GOALS = """MODELING GOALS, in order of importance (when two conflict, the lower-numbered one wins):
1. Accurately reflect confidence and uncertainty. Do not overclaim and do not over-extrapolate. Use all available data to reach the best judgement, and reflect the uncertainty honestly where the data do not help (wider, justified intervals beat sharper numbers that lean on an extrapolation).
2. Get the most accurate forecast.
3. Build a forecast that withstands very heavy scrutiny. The audience is trained in biostatistics and pharma development: use the choices that are idiomatic in pharma modelling (normal approximation / exact tests for responder differences as pre-specified, non-responder imputation, a standard copula or latent-variable construction for correlated co-primaries, mixture or hierarchical priors for translation/attenuation, explicit multiplicity rules), and justify any departure.
4. Build a forecast that is easily communicable. Do not needlessly overcomplicate when a simple choice does the trick; every added mechanism must earn its place by changing the answer or its credibility."""

STANDING_MODEL_RULES = """Two standing rules follow from goals 1-2: price each mechanism once (a discount calibrated on precedent pairs already contains those pairs' winner's curse and design changes; never stack a second term for the same gap), and put caution into the width of a distribution, not into its centre (a judged adjustment without an estimator is centred at no shift: 0 for an additive term, 1 for a multiplicative factor)."""

FETCH_HINT = "web_fetch (if you have it and the person allowed it; otherwise you have no page-fetching tool)"


def project(P):
    pj = read_json(Path(P) / "project.json")
    if not isinstance(pj, dict):
        raise SystemExit("project.json is missing: the scope step must write it first")
    return pj


def sections_of(pj):
    return [s for s in (pj.get("report") or {}).get("sections") or [] if isinstance(s, dict) and s.get("id")]


def sections_text(pj):
    return "\n".join(f"- {s['id']}: {s.get('title', '')} — {s.get('question', '')}" for s in sections_of(pj))


def prompt_vars(P, pj):
    """The placeholders every prompt may use, filled from project.json."""
    ind = pj.get("indication") or {}
    ke = ind.get("key_endpoint") or {}
    a = pj.get("asset") or {}
    R = pj.get("report") or {}
    public = bool(a.get("public"))
    asset = ((f"{a.get('name') or 'the asset'}" + (f" ({a['code']})" if a.get("code") else "")
              + (f", {a['developer']}" if a.get("developer") else "")) if public
             else "the program described in the intake documents (confidential)")
    if public:
        conf = (f"The asset under evaluation is public: {asset}. Name it freely in searches and questions. "
                "Intake documents may still hold non-public notes — never quote anything that is not traceable to a public source.")
    else:
        conf = ("The asset under evaluation is CONFIDENTIAL: the intake documents never leave no-web steps. Write every search, "
                "question and note in PUBLIC terms only — name public drugs, mechanisms, endpoints, geographies, regulators; never the "
                "program, its sponsor, its doses or its target-profile numbers.")
    goals = (pj.get("simulation") or {}).get("goals")
    modeling = (("MODELING GOALS, in order of importance (when two conflict, the lower-numbered one wins):\n"
                 + "\n".join(f"{i}. {g}" for i, g in enumerate(goals, 1))) if goals else DEFAULT_MODELING_GOALS) + "\n" + STANDING_MODEL_RULES
    style = (SKILL / "prompts" / "STYLE_GUIDE.md")
    tgt = pj.get("target") or {}
    reg = pj.get("registry") or {}
    return {
        "indication": ind.get("name", ""),
        "indication_aliases": ", ".join(ind.get("aliases") or []),
        "reference_population": ind.get("reference_population", ""),
        "key_endpoint_label": ke.get("label", ""),
        "key_endpoint_description": ke.get("description", ""),
        "key_endpoint_keywords": ", ".join(ke.get("keywords") or []),
        "key_timepoint_weeks": ke.get("typical_timepoint_weeks") or "not fixed",
        "target": tgt.get("name", ""),
        "target_aliases": ", ".join(tgt.get("aliases") or []),
        "asset": asset,
        "asset_name": (a.get("name") or "the asset") if public else "the asset",
        "asset_public": "yes" if public else "no",
        "asset_modality": a.get("modality", ""),
        "asset_stage": a.get("stage", ""),
        "confidentiality_rule": conf,
        "audience": R.get("audience", ""),
        "scope_note": R.get("scope_note", ""),
        "goal": R.get("goal", ""),
        "decision": R.get("decision", ""),
        "flip_criteria": R.get("flip_criteria") or "a recommendation in the report's options (trial design, endpoint, geography, N/duration, go/no-go reading)",
        "report_kind": REPORT_KIND,
        "tier_legend": "; ".join(f"{k} = {UNIVERSE[k]['label']}" for k in ("A", "B", "C")),
        "fetch_tools": FETCH_HINT,
        "program_horizon": "",
        "local_sources": "",
        "modeling_goals": modeling,
        "style_guide": style.read_text(encoding="utf-8") if style.is_file() else "",
        "sections": sections_text(pj),
        "mechanism_classes": "\n".join(f"- {m}" for m in list(reg.get("mechanism_classes") or []) + ["other", "non-drug", "unknown"]),
        "analog_features": "\n".join(f"- {f}" for f in reg.get("analog_features") or []),
        "target_description": (pj.get("simulation") or {}).get("target_description") or "the trial the report is about",
        "success_definition": (pj.get("simulation") or {}).get("success_definition") or "the trial meets its primary endpoint(s) per the pre-specified analysis",
        "considerations": "\n".join(f"- {c}" for c in (pj.get("simulation") or {}).get("considerations") or []) or "- (none)",
    }


PLACEHOLDER = re.compile(r"\{\{(\w+)\}\}")


def render_prompt(P, name, out, **kw):
    """prompts/<name>.md with every placeholder filled (explicit values win over the project's), written to out.
    Returns the placeholders left unfilled (none, when the call is right)."""
    pj = project(P)
    full = {**prompt_vars(P, pj), **{k: v for k, v in kw.items() if v is not None}}
    txt = source_note(P, name) + (SKILL / "prompts" / f"{name}.md").read_text(encoding="utf-8")
    # one pass, so text pasted in as a value (a critique, a record) is never itself searched for placeholders
    left = sorted({m.group(1) for m in PLACEHOLDER.finditer(txt) if m.group(1) not in full})
    txt = PLACEHOLDER.sub(lambda m: str(full[m.group(1)]) if m.group(1) in full else m.group(0), txt)
    write_text(out, txt)
    return left


def rel(P, path):
    try:
        return str(Path(path).relative_to(P))
    except ValueError:
        return str(path)


# A minimal checker for the JSON schemas in reference/schemas.json (type, properties, required, enum, items).
def schema_errors(obj, schema, where="$", out=None):
    out = [] if out is None else out
    if len(out) >= 12:
        return out
    t = schema.get("type")
    types = t if isinstance(t, list) else ([t] if t else [])
    pyt = {"object": dict, "array": list, "string": str, "boolean": bool, "null": type(None)}

    def is_type(x, ty):
        if ty == "integer":
            return isinstance(x, int) and not isinstance(x, bool)
        if ty == "number":
            return isinstance(x, (int, float)) and not isinstance(x, bool)
        return isinstance(x, pyt.get(ty, object))

    if types and not any(is_type(obj, ty) for ty in types):
        out.append(f"{where}: expected {'/'.join(types)}, got {type(obj).__name__}")
        return out
    if "enum" in schema and obj not in schema["enum"]:
        out.append(f"{where}: {obj!r} is not one of {schema['enum']}")
    if isinstance(obj, dict):
        for k in schema.get("required") or []:
            if k not in obj:
                out.append(f"{where}: missing {k}")
        for k, sub in (schema.get("properties") or {}).items():
            if k in obj:
                schema_errors(obj[k], sub, f"{where}.{k}", out)
    if isinstance(obj, list) and isinstance(schema.get("items"), dict):
        for i, x in enumerate(obj):
            schema_errors(x, schema["items"], f"{where}[{i}]", out)
    return out


def schemas():
    return read_json(SKILL / "reference" / "schemas.json", {}) or {}


def wrong_shape(obj, schema_name):
    """The schema errors that would make a tool fail: a field that should be a list or an object and is not. A missing
    field or an unlisted word is let through."""
    return [e for e in schema_errors(obj, schemas()[schema_name]) if ": expected array" in e or ": expected object" in e]


def load_checked(path, schema_name):
    """(object, errors): the JSON in path checked against reference/schemas.json's schema_name. A file that is not valid
    JSON or does not match is moved aside to <name>.rejected.json, so that on the next launch its step runs again
    instead of being skipped for having a result file."""
    path = Path(path)
    obj = read_json(path)
    errs = [f"{path} is missing or is not valid JSON"] if obj is None else schema_errors(obj, schemas()[schema_name])
    if errs and path.is_file():
        path.replace(path.with_name(path.stem + ".rejected.json"))
        errs.append(f"moved it to {path.stem}.rejected.json, so the step runs again on the next launch")
    return obj, errs


def cmd_setup(a):
    """project.json -> scope.json and scope_public.json (the registry part, in the shape the registry commands read);
    run_settings.json records whether outside sources are allowed."""
    P = Path(a.project)
    write_json(P / "run_settings.json", {"outside": a.outside == "yes"})
    pj = project(P)
    errs = schema_errors(pj, schemas()["project"])
    if errs:
        return fail("project.json does not match the project schema: " + "; ".join(errs[:8]))
    f = (pj.get("focal_trial") or {}).get("nct_id", "")
    if not NCT.fullmatch(str(f)):
        return fail(f"project.json focal_trial.nct_id {f!r} is not an NCT id")
    if a.trial and f.upper() != a.trial.upper():
        return fail(f"project.json is set up for {f}, not for {a.trial}: a project.json packed for another question? Remove it and launch again")
    if not sections_of(pj):
        return fail("project.json has no report sections")
    reg = pj.get("registry") or {}
    scope = {"title": pj.get("title", ""), "scope_label": pj.get("scope_label", ""), "public_thesis": pj.get("public_thesis", ""),
             "horizon": "", "target": pj.get("target") or {"name": "", "aliases": []}, "indications": reg.get("indications") or [],
             "mechanism_classes": reg.get("mechanism_classes") or [], "analog_features": reg.get("analog_features") or [],
             "pathway_interventions": reg.get("pathway_interventions") or [], "focal_nct": f, "notes": pj.get("notes", "")}
    if not scope["indications"]:
        return fail("project.json registry.indications is empty")
    old = read_json(P / "scope.json")
    if isinstance(old, dict) and old.get("indications"):
        # keep the ids plan gave the indications and any indications the research loop added
        scope["indications"] = old["indications"]
    write_json(P / "scope.json", scope)
    write_json(P / "scope_public.json", {k: v for k, v in scope.items() if k != "notes"})
    return done({"ok": True, "focal_nct": f, "sections": [s["id"] for s in sections_of(pj)],
                 "indications": [i.get("name") for i in scope["indications"]]})


def cmd_classify_prompt(a):
    """The classify prompt for one trial: prompts/classify_trial.md with the trial's record, cut by compact_record."""
    P = Path(a.project)
    n = a.nct.strip().upper()
    rec = load_record(P, n)
    if rec is None:
        return fail(f"raw/ctgov/{n}.json is missing or says the trial was not found")
    out = P / "classify" / "prompts" / f"{n}.md"
    left = render_prompt(P, "classify_trial", out, record=json.dumps(compact_record(rec), ensure_ascii=False, indent=1))
    return done({"ok": True, "prompt": str(out), "unfilled": left})


def compact_record(rec, max_chars=3500):
    """The fields of a record the classifier sees, long text cut at max_chars."""
    t = rec.get("trial") or {}

    def c(s):
        s = s or ""
        return s if len(s) <= max_chars else s[:max_chars] + " …[truncated]"

    def outs(key):
        return [{"measure": o.get("measure"), "time_frame": o.get("time_frame")} for o in (t.get(key) or []) if isinstance(o, dict)][:12]

    return {"nct_id": t.get("nct_id"), "title": t.get("title"), "brief_title": t.get("brief_title"), "acronym": t.get("acronym"),
            "status": t.get("status"), "phase": t.get("phase"), "study_type": t.get("study_type"), "conditions": t.get("conditions"),
            "interventions": t.get("interventions"), "sponsor": t.get("sponsor"), "collaborators": t.get("collaborators"),
            "enrollment": t.get("enrollment"), "start_date": t.get("start_date"), "primary_completion_date": t.get("primary_completion_date"),
            "brief_summary": c(t.get("brief_summary")), "detailed_description": c(t.get("detailed_description")),
            "eligibility_criteria": c(t.get("eligibility_criteria")), "minimum_age": t.get("minimum_age"), "maximum_age": t.get("maximum_age"),
            "sex": t.get("sex"), "healthy_volunteers": t.get("healthy_volunteers"), "primary_outcomes": outs("primary_outcomes"),
            "secondary_outcomes": outs("secondary_outcomes"), "n_locations": len(t.get("locations") or [])}


def cut(s, n):
    s = str(s or "")
    return s if len(s) <= n else s[: n - 1] + "…"




# ---- the notes put in front of a step's prompt: the local source library, and what to use without outside sources --
def library_note(pre=""):
    pre_note = f" (put {pre} in front of it from where you work)" if pre else ""
    return f"""LOCAL SOURCE LIBRARY (where present): the documents the person supplied for this run, to search before the web
(or instead of it, in a run without outside sources).
- {pre}sources/INDEX.tsv: one line per document, with columns id, kind, url_or_id, host, first_fetch_utc, chars, path,
  note. Start here. Its path column is relative to the project folder{pre_note}.
- {pre}sources/docs/: the documents themselves (papers, protocols, statistical analysis plans, regulatory documents).
  PDFs are kept as PDF: open them with your file-reading tool, since Grep does not search inside them.
Any other document folder a prompt names does not exist in this run.
Search with Grep (a regular expression over file contents: a drug name, an endpoint, an NCT id, a number) and Glob
(file names), then Read what you need. Cite a document by its path, and by its original URL when the document gives one.

"""


def no_web_note(pre=""):
    return f"""TOOLS IN THIS ENVIRONMENT: you have the project files only — no web search, no connectors (clinical_trials, pubmed,
biorxiv, chembl, fda), no page fetching. Where the task below tells you to search the web, query a connector or fetch a
page, use the project files instead: {pre}trials/<NCT>.json (registry records incl. posted results),
{pre}trials_index.csv, {pre}analyses/, {pre}knowledge/, {pre}research/ and the local source library below. If the answer
needs a source that is not on disk, do not guess: say so in the answer's gaps, and add a paper_requests entry when a
specific paper would settle it.

""" + library_note(pre)


def planner_note(pre=""):
    return f"""SOURCES IN THIS ENVIRONMENT: the researchers who answer your questions have the project files only — no web
search, no connectors, no page fetching. Ask only questions the project files can answer: {pre}trials/<NCT>.json
(registry records incl. posted results), {pre}analyses/, {pre}topics/, {pre}knowledge/, {pre}research/ and the local
source library below. Where a decision needs a document that is not on disk, say so in blocked_items (or as a paper the
project should get), rather than planning a question that needs the web.

""" + library_note(pre)


# Which prompts assume web tools, which plan for researchers, and which get no note (the model builder works only from
# its parameter sheet; the sheet's brief is a file in its folder; the scope step is briefed directly).
WEB_PROMPTS = {"trial_checklist_worker", "researcher", "card_checker", "verify_claim"}
PLANNER_PROMPTS = {"research_planner", "followup_triage"}
NO_NOTE = {"sim_build", "sim_brief", "sim_finalize", "forecast_scope"}


def outside_sources(P):
    """True when the person allowed an outside registry or literature route for this run (setup --outside yes)."""
    return bool((read_json(Path(P) / "run_settings.json", {}) or {}).get("outside"))


def source_note(P, name):
    """The note put in front of a task's prompt: with no outside sources, the no-web note for a task whose
    prompt assumes web tools, the planner note for the planner and the triage, the library alone for any other; with
    outside sources, the library alone."""
    if name in NO_NOTE:
        return ""
    if outside_sources(P):
        return library_note()
    if name in WEB_PROMPTS:
        return no_web_note()
    if name in PLANNER_PROMPTS:
        return planner_note()
    return library_note()


# ---- what earlier launches finished (so a relaunch does not spend a workflow step on each finished helper) ----------
INVENTORY_GLOBS = ("project.json", "raw/ctgov/*.json", "registry/search/*.json", "steps/classify-*.json", "analyses/*/card.md",
                   "analyses/*/check.json", "analyses/*/prompts/check.md", "knowledge/trials_digest.md", "research/iter*/plan.json",
                   "research/iter*/followups.json", "research/iter*/q*.json", "research/iter*/integrate__*.json",
                   "sim/history/round*/params.json", "sim/history/round*/critique.json", "report/out/*.json", "report/claims.json",
                   "checks/claims/c*.json")


def cmd_inventory(a):
    """The result files already on disk, among those the workflow's helpers write, as paths relative to the project
    folder (the script puts the folder in front). The script skips a helper whose result file is listed, as the helper
    itself would, without spending a step. Of the registry records (raw/ctgov/) only the --trial record is listed, since
    it is the only one a helper has as its result file; listing every record would make the list far too long."""
    P = Path(a.project)
    out = set()
    for g in INVENTORY_GLOBS:
        if g.startswith("raw/ctgov/"):
            g = f"raw/ctgov/{a.trial}.json" if NCT.fullmatch(a.trial or "") else ""
        for f in (P.glob(g) if g else []):
            if f.is_file():
                out.add(f.relative_to(P).as_posix())
    return done({"ok": True, "files": sorted(out)})
