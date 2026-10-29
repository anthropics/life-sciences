# Copyright 2026 Anthropic PBC
# SPDX-License-Identifier: Apache-2.0
"""Registry commands of the clinical-trial-forecast tools: the run's files, the query plan, candidates, index, tiers and trial selection."""

import collections
import hashlib
import json
import re
import shutil
import subprocess
import sys
from pathlib import Path

from ft_base import *  # noqa: F401,F403

# ---- plan ------------------------------------------------------------------------------------------------------------
def api_url(params, page_size=1000):
    adv = [params["advanced_query"]] if params.get("advanced_query") else []
    if params.get("study_type"):
        adv.append(f"AREA[StudyType]{params['study_type']}")
    if params.get("phase"):
        adv.append("(" + " OR ".join(f"AREA[Phase]{x}" for x in params["phase"]) + ")")
    q = {"format": "json", "pageSize": str(max(10, min(1000, page_size))), "fields": "NCTId", "countTotal": "true"}
    if params.get("condition"):
        q["query.cond"] = params["condition"]
    if params.get("intervention"):
        q["query.intr"] = params["intervention"]
    if params.get("status"):
        q["filter.overallStatus"] = ",".join(params["status"])
    if adv:
        q["filter.advanced"] = " AND ".join(adv)
    return API + "?" + urllib.parse.urlencode(q)


def cmd_plan(a):
    P = Path(a.project)
    S = size_of(a.size)
    scope = read_json(P / "scope.json")
    if not isinstance(scope, dict) or not scope.get("indications"):
        return fail("scope.json is missing or has no indications: the scope step must write it first")
    U = UNIVERSE
    itypes = " OR ".join(f"AREA[InterventionType]{x}" for x in U["intervention_types"])
    A, B = U["A"], U["B"]
    cap = a.per_indication or S["per_indication"]
    inds = [i for i in scope["indications"] if not i.get("not_pulled")]
    inds = inds[: S["max_indications"]] if S["max_indications"] else inds
    dropped = [i.get("name") for i in scope["indications"] if i not in inds and not i.get("not_pulled")]
    groups, seen_ids = [], {PATHWAY}
    for ind in inds:
        cond = str(ind.get("ctgov_condition_query") or "").strip()
        if not cond:
            return fail(f"indication {ind.get('id')!r} has no ctgov_condition_query; fix scope.json")
        iid = unique_id(slug(ind.get("id") or ind.get("name"), "indication"), seen_ids)
        ind["id"] = iid  # written back below, so the classify step and the index use the same ids
        groups.append({"id": iid, "indication": iid, "cap": cap, "queries": [
            {"family": "A-active", "params": {"condition": cond, "status": A["statuses"], "study_type": U["study_type"], "advanced_query": f"({itypes})"}},
            {"family": "A-completed-recent", "params": {"condition": cond, "status": ["COMPLETED"], "study_type": U["study_type"],
                                                        "advanced_query": f"({itypes}) AND AREA[PrimaryCompletionDate]RANGE[{A['also_completed_since']},MAX]"}},
            {"family": "B-evidence", "params": {"condition": cond, "status": B["statuses"], "phase": B["phases"], "study_type": U["study_type"],
                                                "advanced_query": f"({itypes}) AND AREA[StartDate]RANGE[{B['start_since']},MAX]"}}]})
    drugs = [x.strip() for x in scope.get("pathway_interventions") or [] if isinstance(x, str) and x.strip()]
    drugs = drugs[: S["max_drugs"]] if S["max_drugs"] else drugs
    if drugs:
        per_drug = a.per_drug or S["per_drug"] or max(40, (cap or 600) // len(drugs))
        for drug in drugs:
            iq = f'"{drug}"' if " " in drug else drug
            groups.append({"id": unique_id("pathway_" + slug(drug, "drug"), seen_ids), "indication": PATHWAY, "cap": per_drug, "drug": drug, "queries": [
                {"family": "T-pathway-active", "params": {"intervention": iq, "status": A["statuses"], "study_type": U["study_type"]}},
                {"family": "T-pathway-completed", "params": {"intervention": iq, "status": ["COMPLETED", "TERMINATED"], "study_type": U["study_type"],
                                                             "advanced_query": f"AREA[StartDate]RANGE[{B['start_since']},MAX]"}}]})
    for g in groups:
        for q in g["queries"]:
            q["api_url"] = api_url(q["params"], g["cap"] or 1000)
    plan = {"size": a.size, "rules": "per group: the union of its queries' NCT ids in the order listed, the first cap kept (sizes small "
            "and medium first give each query an equal share of the cap, so completed evidence trials are not crowded out); a trial's "
            "indications are the groups that kept it (pathway groups count as target_pathway). Each api_url asks for at most cap ids, "
            "so one page is enough.", "groups": groups}
    plan["queries_sha"] = sha(groups)
    write_json(P / "scope.json", scope)
    # The helpers that search the web read this copy: the scope prompt asks for public wording in every field but notes.
    write_json(P / "scope_public.json", {k: v for k, v in scope.items() if k != "notes"})
    write_json(P / "registry" / "queries.json", plan)
    return done({"ok": True, "groups": [g["id"] for g in groups], "searches": sum(len(g["queries"]) for g in groups),
                 "indications_left_out": dropped, "pathway_drugs_left_out": [x for x in scope.get("pathway_interventions") or [] if x not in drugs],
                 "summary": f"{len(groups)} search groups, {sum(len(g['queries']) for g in groups)} searches"})


# ---- candidates ------------------------------------------------------------------------------------------------------
# The fields the record shape needs (from_api_v2 reads these), so a fetched record stays short.
RECORD_FIELDS = ("NCTId,BriefTitle,OfficialTitle,Acronym,OverallStatus,Phase,StudyType,Condition,InterventionName,LeadSponsorName,"
                 "CollaboratorName,EnrollmentCount,StartDate,PrimaryCompletionDate,CompletionDate,BriefSummary,EligibilityCriteria,"
                 "MinimumAge,MaximumAge,Sex,HealthyVolunteers,PrimaryOutcome,SecondaryOutcome,LocationFacility,LocationCity,"
                 "LocationState,LocationCountry,LocationStatus,HasResults")


def record_url(n):
    return f"{API}/{n}?" + urllib.parse.urlencode({"format": "json", "fields": RECORD_FIELDS})


def cmd_record_url(a):
    """The record address the classify step fetches for one trial (used by the lead's probe before launch)."""
    n = a.nct.strip().upper()
    if not NCT.fullmatch(n):
        return fail(f"{a.nct!r} is not an NCT id")
    return done({"ok": True, "record_url": record_url(n), "fallback_url": f"{API}/{n}?format=json"})


def cmd_candidates(a):
    P = Path(a.project)
    S = size_of(a.size)
    plan = read_json(P / "registry" / "queries.json")
    if not plan:
        return fail("registry/queries.json is missing: run plan first")
    kept, fams, order, missing, sources, bad, truncated = {}, {}, [], [], collections.Counter(), 0, []
    notes = []
    for g in plan["groups"]:
        got = read_json(P / "registry" / "search" / f"{g['id']}.json")
        local = read_json(P / "registry" / "search" / f"{g['id']}.local.json")
        local_ids = [str(n).upper() for q in ((local or {}).get("queries") or []) if isinstance(q, dict) for n in q.get("nct_ids") or []] if isinstance(local, dict) else []
        if local_ids:
            sources["attached files"] += 1
            for n in local_ids:
                if NCT.fullmatch(n):
                    fams.setdefault(n, [])
                    if "local" not in fams[n]:
                        fams[n].append("local")
        if not isinstance(got, dict):
            if not local_ids:
                missing.append(g["id"])
                continue
            got = {"queries": []}
        by_family = {q.get("family"): q for q in got.get("queries") or [] if isinstance(q, dict)}
        seen, per_family = [], []
        for q in g["queries"]:
            r = by_family.get(q["family"]) or {}
            if got.get("queries"):
                sources[r.get("source") or "missing"] += 1
            ids = r.get("nct_ids") or []
            if isinstance(ids, str):
                ids = re.findall(r"NCT\d{8}", ids.upper())
            total = to_int(r.get("total"))
            if total and len(ids) < min(total, int(g.get("cap") or total)):
                truncated.append(f"{g['id']}|{q['family']}")
            fam_ids = []
            for n in ids:
                n = str(n).strip().upper()
                if not NCT.fullmatch(n):
                    bad += 1
                    continue
                fams.setdefault(n, [])
                if q["family"] not in fams[n]:
                    fams[n].append(q["family"])
                if n not in seen:
                    seen.append(n)
                if n not in fam_ids:
                    fam_ids.append(n)
            per_family.append(fam_ids)
        cap = int(g.get("cap") or 0)
        # An equal share of the cap for each query first, so that a large active-trial search cannot crowd out the
        # completed evidence trials when the cap is small; then the rest of the union in order. With no cap, all of it.
        keep = []
        if cap and plan.get("size") in ("small", "medium"):
            share = -(-cap // max(1, len(per_family)))
            for fam_ids in per_family:
                for n in fam_ids[:share]:
                    if n not in keep and len(keep) < cap:
                        keep.append(n)
            for n in seen:
                if len(keep) >= cap:
                    break
                if n not in keep:
                    keep.append(n)
        else:  # the first cap ids of the union, in query order
            keep = seen[:cap] if cap else seen
        # attached records that match the group come after the searches' ids, within the same cap
        for n in local_ids:
            if NCT.fullmatch(n) and n not in keep and (not cap or len(keep) < cap):
                keep.append(n)
        left_out = sum(1 for n in local_ids if NCT.fullmatch(n) and n not in keep)
        if left_out:
            notes.append(f"{g['id']}: the cap of {cap} left out {left_out} attached records")
        for n in keep:
            if n not in kept:
                kept[n] = []
                order.append(n)
            if g["indication"] not in kept[n]:
                kept[n].append(g["indication"])
    # the trial the forecast is about is always a candidate, first, whatever the searches returned
    scope = read_json(P / "scope.json", {}) or {}
    focal = str(scope.get("focal_nct") or "").upper()
    if NCT.fullmatch(focal):
        first_ind = ((scope.get("indications") or [{}])[0] or {}).get("id") or "focal"
        if focal in kept:
            order.remove(focal)
        else:
            kept[focal] = [first_ind]
            fams[focal] = ["forecast-trial"]
        order.insert(0, focal)
    cands = [{"nct_id": n, "indications": kept[n], "query_families": fams.get(n, []), "record_url": record_url(n)} for n in order]
    write_json(P / "registry" / "candidates.json", {"pulled_at": now_iso(), "queries_sha": plan.get("queries_sha"),
                                                      "sources": dict(sources), "candidates": cands})
    bdir = P / "registry" / "batches"
    bdir.mkdir(parents=True, exist_ok=True)
    for f in bdir.glob("batch-*.json"):
        f.unlink()
    size = a.batch or S["batch"]
    batches = []
    for i in range(0, len(cands), size):
        # the name carries a hash of the batch's trials, so a batch whose contents change gets a new name (and the classify
        # step, whose result file is named after the batch, runs for it again)
        h = hashlib.sha256(" ".join(c["nct_id"] for c in cands[i:i + size]).encode()).hexdigest()[:8]
        f = bdir / f"batch-{i // size + 1:02d}-{h}.json"
        write_json(f, {"batch": i // size + 1, "trials": cands[i:i + size]})
        batches.append(str(f))
    return done({"ok": bool(cands), "candidates": len(cands), "batches": batches, "groups_without_results_file": missing,
                 "searches_with_fewer_ids_than_total": truncated, "ignored_ids": bad, "sources": dict(sources), "step_notes": notes,
                 "problem": "" if cands else "no NCT ids found in registry/search/*.json"})




# ---- files the person attached: registry records and documents -------------------------------------------------------
DOC_KINDS = {".pdf": "pdf", ".txt": "text", ".md": "text", ".html": "html", ".htm": "html", ".csv": "table", ".tsv": "table", ".docx": "docx"}


def records_in(obj):
    """Registry records in one attached JSON file: an API v2 study, an API v2 search result ({"studies": [...]}), a
    record in this skill's shape, a bare trial object, or a list of any of these."""
    if isinstance(obj, list):
        for x in obj:
            try:
                yield from list(records_in(x))
            except Exception:  # one malformed study must not cost the file's other studies
                pass
        return
    if not isinstance(obj, dict):
        return
    if isinstance(obj.get("raw_json"), str):  # a row of a registry dump that keeps each study's API v2 JSON as text
        try:
            yield from records_in(json.loads(obj["raw_json"]))
        except ValueError:
            pass
        return
    if isinstance(obj.get("studies"), list):
        yield from records_in(obj["studies"])
        return
    if "protocolSection" in obj:
        t = from_api_v2(obj)
        t["nct_id"] = str(t.get("nct_id") or "").strip()
        if obj.get("resultsSection"):
            t["results"] = obj["resultsSection"]  # the posted results, which the prompts read as part of the record
        yield {"nct_id": t["nct_id"], "found": True, "source": "attached (ClinicalTrials.gov API v2)", "fetched_at": obj.get("fetched_at"), "trial": t}
        return
    t = obj.get("trial") if isinstance(obj.get("trial"), dict) else (obj if obj.get("nct_id") else None)
    if isinstance(t, dict):
        if "protocolSection" in t:
            res = t.get("resultsSection")
            t = from_api_v2(t)
            if res:
                t["results"] = res
        t = normalise_trial(t)
        t["nct_id"] = str(t.get("nct_id") or obj.get("nct_id") or "").strip()
        yield {"nct_id": t["nct_id"], "found": True, "source": obj.get("source") or "attached",
               "fetched_at": obj.get("fetched_at"), "trial": t}


def cmd_import_local(a):
    """The run's files (folders searched recursively: the files attached to the task, and the data/ folder packed into
    the skill) -> registry records under raw/ctgov/ (never over a record already there) and documents under
    sources/docs/, listed in sources/INDEX.tsv. A missing folder is not an error: it just holds no files."""
    P = Path(a.project)
    recs, docs, skipped, bad = [], [], [], []
    idx_rows = []
    pairs = []
    for s in a.source:
        src = Path(s)
        if src.is_dir():
            found = [x for x in sorted(src.rglob("*")) if x.is_file() and not any(p.startswith(".") for p in x.relative_to(src).parts)]
            pairs += [(src, x) for x in found if not x.is_symlink()]
            bad += [rel(src, x) + " (a symbolic link: not followed)" for x in found if x.is_symlink()]
    packed_project = ""
    for src, f in pairs:
        try:  # one malformed file must not stop the import
            if f.name == "project.json" and f.parent == src:
                # a project setup packed with the data (prepared outside the run): used in place of the scope step's
                if not (P / "project.json").exists():
                    shutil.copy2(f, P / "project.json")
                    packed_project = rel(src, f)
                continue
            if f.suffix.lower() in (".json", ".jsonl", ".ndjson"):
                obj = read_json(f)
                if obj is None and f.suffix.lower() != ".json":  # one study per line
                    obj = []
                    for line in f.read_text(encoding="utf-8", errors="ignore").splitlines():
                        try:
                            obj.append(json.loads(line))
                        except ValueError:
                            pass
                got = list(records_in(obj)) if obj is not None else []
                got = [r for r in got if NCT.fullmatch(str(r.get("nct_id") or ""))]
                if not got:
                    bad.append(rel(src, f))
                    continue
                for r in got:
                    n = r["nct_id"]
                    dst = P / "raw" / "ctgov" / f"{n}.json"
                    old = read_json(dst)
                    if isinstance(old, dict) and old.get("found") is not False:
                        skipped.append(n)
                        continue
                    write_json(dst, r)
                    recs.append(n)
            elif f.suffix.lower() in DOC_KINDS:
                dst = P / "sources" / "docs" / f.name
                k = 2
                while dst.exists() and dst.read_bytes() != f.read_bytes():  # two different files with one name
                    dst = dst.with_name(f"{f.stem}_{k}{f.suffix}")
                    k += 1
                if not dst.exists():
                    dst.parent.mkdir(parents=True, exist_ok=True)
                    shutil.copy2(f, dst)
                chars = len(dst.read_text(encoding="utf-8", errors="ignore")) if DOC_KINDS[f.suffix.lower()] in ("text", "html", "table") else ""
                idx_rows.append([slug(f.stem, "doc"), DOC_KINDS[f.suffix.lower()], "", "attached", now_iso(), str(chars), re.sub(r"[\t\r\n]", " ", rel(P, dst)),
                                 "attached by the person; cite it by this path"])
                docs.append(rel(P, dst))
            else:
                bad.append(rel(src, f))
        except Exception as e:
            bad.append(rel(src, f) + " (" + type(e).__name__ + ")")
    if idx_rows:
        idx = P / "sources" / "INDEX.tsv"
        have = set()
        lines = []
        if idx.is_file():
            lines = idx.read_text(encoding="utf-8").splitlines()
            have = {ln.split("\t")[6] for ln in lines[1:] if len(ln.split("\t")) > 6}
        if not lines:
            lines = ["\t".join(["id", "kind", "url_or_id", "host", "first_fetch_utc", "chars", "path", "note"])]
        lines += ["\t".join(r) for r in idx_rows if r[6] not in have]
        write_text(idx, "\n".join(lines) + "\n")
    return done({"ok": True, "records": len(recs), "records_already_there": len(skipped), "documents": len(docs),
                 "not_used": bad[:40], "record_ids": recs[:60], "project_json": packed_project})


def _terms(ind):
    out = [ind.get("name") or ""] + list(ind.get("aliases") or [])
    q = str(ind.get("ctgov_condition_query") or "")
    out += [x.strip().strip('"()') for x in re.split(r"\bOR\b|\bAND\b", q)]
    return [t.lower() for t in out if len(t.strip()) >= 3]


def local_rank(n, t):
    """Sort key for attached records: active, recently completed, then phase 2/3 evidence, newest first in each."""
    A, B = UNIVERSE["A"], UNIVERSE["B"]
    status, pcd, start = t.get("status"), date(t.get("primary_completion_date")), date(t.get("start_date"))
    if status in A["statuses"]:
        rank = 0
    elif status == "COMPLETED" and pcd and pcd >= date(A["also_completed_since"]):
        rank = 1
    elif status in B["statuses"] and start and start >= date(B["start_since"]) and set(t.get("phase") or []) & set(B["phases"]):
        rank = 2
    else:
        rank = 3
    return rank, -(start.toordinal() if start else 0), n


def cmd_local_search(a):
    """For each group of registry/queries.json, the attached records that match it -> registry/search/<group>.local.json:
    an indication group takes records whose conditions or title name the indication (its name, aliases or the terms of
    its condition query); a pathway group takes records that list the drug among their interventions."""
    P = Path(a.project)
    plan = read_json(P / "registry" / "queries.json")
    if not plan:
        return fail("registry/queries.json is missing: run plan first")
    scope = read_json(P / "scope.json", {}) or {}
    inds = {i.get("id"): i for i in scope.get("indications") or []}
    recs = {}
    for f in sorted((P / "raw" / "ctgov").glob("NCT*.json")):
        r = load_record(P, f.stem)
        if r is not None and str(r.get("source") or "").startswith("attached"):
            recs[f.stem] = r["trial"]
    counts = {}
    for g in plan["groups"]:
        if g["indication"] == PATHWAY:
            drug = str(g.get("drug") or "").lower()
            word = re.compile(r"(?<![a-z0-9])" + re.escape(drug) + r"(?:e?s)?(?![a-z0-9])") if drug else None  # whole word (a plural too): "ate" is not in "acetate"
            ids = [n for n, t in recs.items() if word and any(word.search(str(x).lower()) for x in t.get("interventions") or [])]
        else:
            # whole words or phrases, plurals included: "solid tumor" matches "solid tumors", but "sle" does not match "sleep apnea"
            terms = [re.compile(r"(?<![a-z0-9])" + re.escape(term) + r"(?:e?s)?(?![a-z0-9])") for term in _terms(inds.get(g["indication"]) or {})]
            ids = [n for n, t in recs.items()
                   if any(term.search(" ".join([*(map(str, t.get("conditions") or [])), str(t.get("brief_title") or ""), str(t.get("title") or "")]).lower()) for term in terms)]
        ids.sort(key=lambda n: local_rank(n, recs[n]))
        write_json(P / "registry" / "search" / f"{g['id']}.local.json",
                   {"group": g["id"], "queries": [{"family": "local", "source": "attached files", "total": None, "nct_ids": ids, "note": "records the person attached"}]})
        counts[g["id"]] = len(ids)
    return done({"ok": True, "records_on_disk": len(recs), "matches": counts})
# ---- classification, tiers, index, derived tables -------------------------------------------------------------------
def canon_mech(vocab, label):
    """A free-text mechanism label mapped onto the scope's controlled vocabulary (exact, prefix or token match)."""
    if not label:
        return "unknown"
    lab = str(label).strip().lower()
    for v in vocab:
        if lab == v.lower():
            return v
    for v in vocab:
        vl = v.lower()
        if lab.startswith(vl.split(" (")[0]) or vl.startswith(lab.split(" (")[0]):
            return v
    for v in vocab:
        head = v.lower().split(" (")[0].split(" / ")[0]
        if head and head in lab:
            return v
    return "unknown"


def assign_tiers(t, cls):
    if not cls.get("is_drug_trial"):
        return [], "not a drug trial"
    ke = cls.get("key_endpoint") or {}
    has_ke = bool(ke.get("primary") or ke.get("secondary"))
    if to_int(cls.get("relevance_score")) == 0 and not cls.get("analog_features"):
        return [], "not in any scoped indication"
    status, pcd, start = t.get("status"), date(t.get("primary_completion_date")), date(t.get("start_date"))
    phases = set(t.get("phase") or [])
    A, B = UNIVERSE["A"], UNIVERSE["B"]
    tiers = []
    if status in A["statuses"] or (status == "COMPLETED" and pcd and pcd >= date(A["also_completed_since"])):
        tiers.append("A")
    if (status in B["statuses"] and phases & set(B["phases"]) and start and start >= date(B["start_since"])
            and (has_ke or not B["require_key_endpoint"])):
        tiers.append("B")
    if cls.get("analog_features") and cls.get("novel_investigational"):
        tiers.append("C")
    return tiers, ("" if tiers else "no tier rule matched")


KNOWN_STATUSES = set(ACTIVE) | {"COMPLETED", "TERMINATED", "WITHDRAWN", "SUSPENDED", "UNKNOWN", "NO_LONGER_AVAILABLE",
                                 "AVAILABLE", "TEMPORARILY_NOT_AVAILABLE", "APPROVED_FOR_MARKETING", "WITHHELD"}

INDEX_COLUMNS = ["nct_id", "title", "acronym", "sponsor", "phase", "status", "start_date", "primary_completion_date", "enrollment",
                 "primary_mechanism_class", "interventions", "route", "key_endpoint_primary", "key_endpoint_timepoint_weeks", "indications",
                 "population", "duration_weeks", "tiers", "analog_features", "relevance_score", "n_locations", "n_countries", "top_countries",
                 "has_results", "url", "excluded_reason", "one_line"]


def cmd_index(a):
    P = Path(a.project)
    scope = read_json(P / "scope.json") or {}
    vocab = list(scope.get("mechanism_classes") or []) + ["other", "non-drug", "unknown"]
    cand = read_json(P / "registry" / "candidates.json") or {}
    ind_of = {c["nct_id"]: set(c.get("indications") or []) for c in cand.get("candidates") or []}
    rows, no_record, no_class, recs, odd_status = [], [], [], {}, []
    missing = {"is_drug_trial": None, "primary_mechanism_class": "unknown", "relevance_score": 0, "key_endpoint": {},
               "population": {}, "design": {}, "analog_features": [], "one_line": "", "notes": "classification missing"}

    def add_row(n, rec, cls):
        cls = dict(cls)
        cls["primary_mechanism_class"] = canon_mech(vocab, cls.get("primary_mechanism_class", ""))
        for i in cls.get("interventions") or []:
            if isinstance(i, dict):
                i["mechanism_class"] = canon_mech(vocab, i.get("mechanism_class", ""))
        t = rec["trial"]
        tiers, excl = assign_tiers(t, cls)
        locs = [x for x in t.get("locations") or [] if isinstance(x, dict)]
        countries = collections.Counter(loc.get("country") for loc in locs if loc.get("country"))
        inv = [i for i in cls.get("interventions") or [] if isinstance(i, dict) and i.get("role") == "investigational"]
        pop = cls.get("population") or {}
        ke = cls.get("key_endpoint") or {}
        row = {"nct_id": n, "title": t.get("brief_title") or t.get("title"), "acronym": t.get("acronym"), "sponsor": t.get("sponsor"),
               "phase": "/".join(t.get("phase") or []), "status": t.get("status"), "start_date": t.get("start_date"),
               "primary_completion_date": t.get("primary_completion_date"), "enrollment": t.get("enrollment"),
               "primary_mechanism_class": cls.get("primary_mechanism_class"),
               "interventions": [i.get("canonical_name") for i in inv] or t.get("interventions"),
               "route": "; ".join(sorted({i.get("route") or "" for i in inv} - {""})),
               "key_endpoint_primary": ke.get("primary"), "key_endpoint_timepoint_weeks": ke.get("timepoint_weeks"),
               "indications": sorted(ind_of[n]),
               "population": "; ".join(str(pop[k]) for k in ("severity_criteria", "age_group") if pop.get(k)), "duration_weeks": (cls.get("design") or {}).get("duration_weeks"),
               "tiers": "".join(tiers), "novel_investigational": cls.get("novel_investigational"), "analog_features": cls.get("analog_features"),
               "relevance_score": to_int(cls.get("relevance_score")) or 0, "n_locations": len(locs), "n_countries": len(countries),
               "top_countries": [f"{c} ({k})" for c, k in countries.most_common(5)], "has_results": t.get("has_results"),
               "url": t.get("url"), "excluded_reason": excl, "one_line": cls.get("one_line")}
        write_json(P / "trials" / f"{n}.json", {"nct_id": n, "updated": now_iso(), "source": rec.get("source") or "ClinicalTrials.gov",
                                                "fetched_at": rec.get("fetched_at"), "trial": t, "classification": cls,
                                                "tiers": tiers, "excluded_reason": excl, "index_row": row})
        rows.append(row)
        recs[n] = t

    for n in sorted(ind_of):
        rec = load_record(P, n)
        if rec is None or rec["trial"].get("nct_id") != n:
            no_record.append(n)
            continue
        if str(rec["trial"].get("status")) not in KNOWN_STATUSES:
            odd_status.append(n)
        try:
            add_row(n, rec, read_json(P / "classify" / f"{n}.json"))
        except Exception:  # no classification, or a malformed one: it must not stop the index for every other trial
            no_class.append(n)
            cf = P / "classify" / f"{n}.json"
            if cf.is_file():  # moved aside, as load_checked does, so that the next launch asks for it again
                cf.replace(cf.with_name(f"{n}.rejected.json"))
            try:
                add_row(n, rec, missing)
            except Exception:
                no_record.append(n)
    rows.sort(key=lambda r: (-(r.get("relevance_score") or 0), r["nct_id"]))
    write_json(P / "trials_index.json", rows)
    write_csv(P / "trials_index.csv", rows, INDEX_COLUMNS)
    incl = [r for r in rows if r["tiers"]]
    md = [f"# Trial universe: {scope.get('title', '')}",
          f"{len(rows)} candidates, {len(incl)} in tiers: A={sum('A' in r['tiers'] for r in incl)} "
          f"B={sum('B' in r['tiers'] for r in incl)} C={sum('C' in r['tiers'] for r in incl)}", "",
          "Tiers: A = competing for the same patients and sites now or recently; B = completed or terminated phase 2/3 trials "
          "with the indication's key endpoint (endpoint and placebo evidence); C = novel programs whose design is a precedent. "
          "A trial can be in more than one.", "", "## Included trials (sorted by relevance)", "",
          md_table(incl, ["nct_id", "tiers", "indications", "relevance_score", "primary_mechanism_class", "interventions", "phase",
                          "status", "sponsor", "enrollment", "n_locations", "one_line"], 70),
          "", "## Excluded candidates", "",
          md_table([r for r in rows if not r["tiers"]], ["nct_id", "title", "excluded_reason", "primary_mechanism_class"], 60)]
    write_text(P / "trials_index.md", "\n".join(md))
    tables = derived_tables(P, rows, recs)
    return done({"ok": bool(rows), "candidates": len(rows), "included": len(incl),
                 "tierA": sum("A" in r["tiers"] for r in incl), "tierB": sum("B" in r["tiers"] for r in incl),
                 "tierC": sum("C" in r["tiers"] for r in incl), "n_records_missing": len(no_record), "records_missing": no_record[:40],
                 "n_classifications_missing": len(no_class), "classifications_missing": no_class[:40], "tables": tables,
                 "unknown_status": odd_status[:40],
                 "problem": "" if rows else "no candidate had a record on disk"})


def derived_tables(P, rows, recs):
    D = P / "derived"
    incl = [r for r in rows if r["tiers"]]
    act = [r for r in incl if "A" in r["tiers"] and r["status"] in ACTIVE]

    def locs(n):
        return [x for x in recs[n].get("locations") or [] if isinstance(x, dict)]

    write_csv(D / "competitors_active.csv", sorted(act, key=lambda r: -(r.get("enrollment") or 0)),
              ["nct_id", "indications", "interventions", "primary_mechanism_class", "novel_investigational", "sponsor", "phase", "status",
               "enrollment", "start_date", "primary_completion_date", "n_locations", "n_countries", "top_countries", "relevance_score",
               "one_line", "url"])
    fac = collections.defaultdict(set)
    ctry, ctry_sites = collections.Counter(), collections.Counter()
    for r in act:
        seen_c = set()
        for loc in locs(r["nct_id"]):
            key = (re.sub(r"\s+", " ", (loc.get("facility") or "").strip().lower()), (loc.get("city") or "").strip().lower(),
                   (loc.get("country") or "").strip())
            if key[0]:
                fac[key].add(r["nct_id"])
            if loc.get("country"):
                ctry_sites[loc["country"]] += 1
                if loc["country"] not in seen_c:
                    ctry[loc["country"]] += 1
                    seen_c.add(loc["country"])
    shared = [{"facility": k[0], "city": k[1], "country": k[2], "n_trials": len(v), "trials": sorted(v)} for k, v in fac.items() if len(v) >= 2]
    shared.sort(key=lambda x: (-x["n_trials"], x["facility"]))
    write_csv(D / "sites_shared.csv", shared, ["facility", "city", "country", "n_trials", "trials"])
    write_csv(D / "countries_active.csv", [{"country": c, "n_active_trials": n, "n_sites": ctry_sites[c]} for c, n in ctry.most_common()],
              ["country", "n_active_trials", "n_sites"])
    sponsors = collections.defaultdict(lambda: {"n_trials": 0, "planned_enrollment": 0, "trials": []})
    for r in act:
        s = sponsors[r.get("sponsor") or "?"]
        s["n_trials"] += 1
        s["planned_enrollment"] += r.get("enrollment") or 0
        s["trials"].append(r["nct_id"])
    write_csv(D / "sponsors_active.csv", [{"sponsor": k, **v} for k, v in sorted(sponsors.items(), key=lambda kv: -kv[1]["planned_enrollment"])],
              ["sponsor", "n_trials", "planned_enrollment", "trials"])
    ep = []
    for r in incl:
        for o in recs[r["nct_id"]].get("primary_outcomes") or []:
            if isinstance(o, dict):
                ep.append({"nct_id": r["nct_id"], "indications": r.get("indications"), "tiers": r["tiers"], "phase": r["phase"],
                           "mechanism": r["primary_mechanism_class"], "enrollment": r["enrollment"], "measure": o.get("measure"),
                           "time_frame": o.get("time_frame")})
    write_csv(D / "primary_endpoints.csv", ep, ["nct_id", "indications", "tiers", "phase", "mechanism", "enrollment", "measure", "time_frame"])
    cc, cs = collections.Counter(), collections.Counter()
    for r in incl:
        seen = set()
        for loc in locs(r["nct_id"]):
            c = loc.get("country")
            if c:
                cs[c] += 1
                if c not in seen:
                    cc[c] += 1
                    seen.add(c)
    write_csv(D / "countries_all_tiers.csv", [{"country": c, "n_trials": n, "n_sites": cs[c]} for c, n in cc.most_common()],
              ["country", "n_trials", "n_sites"])
    no_sites = sum(1 for r in act if not locs(r["nct_id"]))
    return (f"{len(act)} active competitors ({no_sites} without site lists), {len(shared)} facilities shared by 2+ of them, "
            f"{len(cc)} countries across tiered trials")


# ---- trial selection -------------------------------------------------------------------------------------------------
def cmd_select(a):
    """Tiered trials at or above the relevance floor, highest priority first (active competitors, then analogs, then
    evidence trials, later phases, relevance, size); at most per_asset trials per asset; round-robin across indications."""
    P = Path(a.project)
    limit = a.limit or size_of(a.size)["cards"]
    idx = [r for r in read_json(P / "trials_index.json", []) or [] if NCT.fullmatch(str(r.get("nct_id")))]
    rows = [r for r in idx if r.get("tiers") and (r.get("relevance_score") or 0) >= a.min_relevance]

    def prio(r):
        t, st, ph = r["tiers"], r.get("status"), r.get("phase") or ""
        s = 0.0
        if "A" in t and st in ACTIVE:
            s += 3
        if "C" in t:
            s += 2
        if "B" in t:
            s += 1
        if "PHASE3" in ph:
            s += 1
        if "PHASE2" in ph:
            s += 1
        s += (r.get("relevance_score") or 0) * 0.5
        s += min((r.get("enrollment") or 0), 5000) / 5000
        return (-s, r["nct_id"])

    rows.sort(key=prio)
    ok, seen = [], {}
    for r in rows:
        key = str((r.get("interventions") or ["?"])[0]).lower()
        if seen.get(key, 0) >= a.per_asset:
            continue
        seen[key] = seen.get(key, 0) + 1
        ok.append(r)
    queues = {}
    for r in ok:
        queues.setdefault((r.get("indications") or ["?"])[0], []).append(r)
    out = []
    while len(out) < limit and any(queues.values()):
        for k in list(queues):
            if queues[k] and len(out) < limit:
                out.append(queues[k].pop(0))
    # the trial the forecast is about always gets a card, first (it needs no tier)
    focal = str((read_json(P / "scope.json", {}) or {}).get("focal_nct") or "").upper()
    byid = {r["nct_id"]: r for r in idx}
    if focal in byid:
        out = [byid[focal]] + [r for r in out if r["nct_id"] != focal][: max(0, limit - 1)]
    sel = [{"nct_id": r["nct_id"], "tiers": r["tiers"], "indications": r.get("indications"), "one_line": r.get("one_line") or ""} for r in out]
    write_json(P / "analyses" / "selected.json", sel)
    return done({"ok": bool(sel), "selected": sel, "eligible": len(rows), "forecast_trial_has_card": focal in byid,
                 "problem": "" if sel else f"no tiered trial with relevance of at least {a.min_relevance}"})


