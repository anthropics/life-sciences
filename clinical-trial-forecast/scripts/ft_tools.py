#!/usr/bin/env python3
# Copyright 2026 Anthropic PBC
# SPDX-License-Identifier: Apache-2.0
"""Deterministic helpers for the clinical-trial-forecast skill. Python 3 standard library only.

Every command takes --project <folder> (the project folder the workflow uses), except record-url. Every command prints
one JSON object as its last line of output, with "ok" and, when ok is false, "problem". Run each as one plain command,
for example:
    python3 <skill folder>/scripts/ft_tools.py plan --project <run folder>/outputs/work/project --size small

Commands that end in a model step write that step's whole prompt to a file (the step's prompt from prompts/ with every
{{placeholder}} filled from project.json and the files on disk), so a helper only reads the file and follows it.

Commands, in the order a run uses them:
  inventory         the result files earlier launches finished (the workflow skips their helpers)
  record-url        the registry record address for one NCT id (the lead's probe before launch)
  import-local      the run's files (attached to the task, or packed into the skill's data/ folder) -> registry
                    records (raw/ctgov/) and the local source library (sources/)
  setup             project.json (from the scope step) -> scope.json, scope_public.json; records whether outside
                    sources are allowed
  plan              scope.json -> registry/queries.json (the registry searches to run, with ready-made API v2 URLs)
  local-search      the attached records that match each search group -> registry/search/<group>.local.json
  classify-prompt   the classify prompt for one trial
  candidates        registry/search/*.json -> registry/candidates.json + classify batches (the forecast trial first)
  index             raw/ctgov/<NCT>.json + classify/<NCT>.json -> tiers, trials/, trials_index.*, derived/*.csv
  select            trials_index.json -> analyses/selected.json (the forecast trial always first)
  card-prompts      one trial -> the checklist workers' prompts (one per worker group)
  card-apply        the workers' answers -> items/, groups/, then the card writer's prompt
  card-finish       the card writer's answer -> analysis.md, card.md, analysis.json, then the number checker's prompt
  card-check-apply  the number check -> appended to card.md and analysis.json
  digest-prompt     the trials digest prompt
  planner-prompt    the research planner's prompt for one iteration
  questions         after the planner (or the follow-up triage): caps and orders its questions, applies the stop
                    rule, writes each question's researcher prompt
  checks            after the researchers: the number checker's prompt for each answer with key numbers
  followups         the follow-up triage prompt for one iteration
  research-apply    research/iterNN/ -> answer files with their checks, integrator prompts, open questions
  integrate-apply   the integrators' answers -> knowledge/<section>.md
  sim-prep          one model round: the parameter-sheet prompt (when this round needs one)
  sim-spec-apply    the parameter sheet -> checks, sim/params.json, the workspace, the builder's prompt
  sim-build-apply   the builder's result -> the round's copy, the critic's prompt
  sim-critic-apply  the critique -> sim/critique_rN.md; says whether another round is needed
  sim-publish       the last built round -> sim/ and sim/status.json, then the finalize prompt
  sim-finalize-apply the finalize step's files -> sim/
  report-prep       the narrative and appendix writers' prompts (with the model's load-bearing ranking)
  report-body-apply their answers -> report/parts/, figures, then the summary editor's prompt
  report-editor-apply the editor's answer -> summary and methods, then the condense prompts for parts over length
  condense-apply    the condensed parts (kept only when they pass the length and content checks), more condense prompts if needed
  copyedit-prep     the copy-edit prompt for each part
  copyedit-apply    the copy-edited parts (kept only when they pass the length and content checks); assembles report/report.md
  freeze            records the sha256 of report/report.md before its claims are extracted
  claims            after extraction: caps and numbers the claims, writes each checker's prompt
  ledger            report/report.md + report/claims.json + checks/claims/*.json -> report/claim_ledger.{json,md}
  ledger-prep       the ledger-editing prompt for each report part
  ledger-apply      the edited parts (kept only when they pass the length and content checks); reassembles the report
  deliver           the deliverables folder: report, claim ledger, the model and its figures, a one-file summary

These commands hold every rule that needs no model: tier rules, query plan, trial selection, checklist applicability,
stop rule, model rounds, load-bearing ranking, word ceilings, the acceptance checks on condensed, copy-edited and
ledger-edited parts, the claim ledger and footnotes."""

import argparse
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))


from ft_base import *  # noqa: F401,F403
from ft_registry import *  # noqa: F401,F403
from ft_research import *  # noqa: F401,F403
from ft_model import *  # noqa: F401,F403
from ft_report import *  # noqa: F401,F403


# ---- main --------------------------------------------------------------------------------------------------------------
class Parser(argparse.ArgumentParser):
    """Usage errors as the same one-line JSON every command prints."""

    def error(self, message):
        print(json.dumps({"ok": False, "problem": f"{self.prog}: {message}"}))
        sys.exit(2)


def main(argv=None):
    ap = Parser(prog="ft_tools.py", description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True, parser_class=Parser)

    def nct_id(typed):
        s = str(typed).strip().upper()
        if not NCT.fullmatch(s):
            raise argparse.ArgumentTypeError(f"{typed!r} is not an NCT id (NCT and eight digits)")
        return s

    def add(name, fn, project=True, size=False, nct=False, it=False, rnd=False):
        p = sub.add_parser(name)
        if project:
            p.add_argument("--project", required=True)
        if size:
            p.add_argument("--size", default="small", choices=sorted(SIZES))
        if nct:
            p.add_argument("--nct", required=True, type=nct_id)
        if it:
            p.add_argument("--iter", type=int, required=True)
        if rnd:
            p.add_argument("--round", type=int, required=True)
        p.set_defaults(fn=fn)
        return p

    add("record-url", cmd_record_url, project=False, nct=True)
    p = add("inventory", cmd_inventory)
    p.add_argument("--trial", default="", help="the NCT id whose registry record a helper writes (the forecast trial), if any")
    p = add("setup", cmd_setup)
    p.add_argument("--outside", required=True, choices=["yes", "no"])
    p.add_argument("--trial", default="", help="the NCT id the run forecasts; setup fails if project.json is for another trial")
    p = add("import-local", cmd_import_local)
    p.add_argument("--from", dest="source", action="append", required=True, help="a folder of files; repeat for several")
    add("local-search", cmd_local_search)
    p = add("plan", cmd_plan, size=True)
    p.add_argument("--per-indication", type=int, default=0)
    p.add_argument("--per-drug", type=int, default=0)
    p = add("candidates", cmd_candidates, size=True)
    p.add_argument("--batch", type=int, default=0)
    add("classify-prompt", cmd_classify_prompt, nct=True)
    add("index", cmd_index)
    p = add("select", cmd_select, size=True)
    p.add_argument("--limit", type=int, default=0)
    p.add_argument("--min-relevance", type=int, default=3)
    p.add_argument("--per-asset", type=int, default=2)
    add("card-prompts", cmd_card_prompts, nct=True)
    add("card-apply", cmd_card_apply, nct=True)
    add("card-finish", cmd_card_finish, nct=True)
    add("card-check-apply", cmd_card_check_apply, nct=True)
    add("digest-prompt", cmd_digest_prompt)
    p = add("planner-prompt", cmd_planner_prompt, it=True)
    p.add_argument("--k", type=int, required=True)
    p = add("questions", cmd_questions, it=True)
    p.add_argument("--k", type=int, required=True)
    p.add_argument("--wave", type=int, default=1, choices=[1, 2])
    p.add_argument("--min-iterations", type=int, default=2)
    p = add("checks", cmd_checks, it=True)
    p.add_argument("--wave", type=int, default=0)
    p = add("followups", cmd_followups, it=True)
    p.add_argument("--k", type=int, required=True)
    add("research-apply", cmd_research_apply, it=True)
    add("integrate-apply", cmd_integrate_apply, it=True)
    p = add("sim-prep", cmd_sim_prep, rnd=True)
    p.add_argument("--rounds", type=int, required=True)
    add("sim-spec-apply", cmd_sim_spec_apply, size=True, rnd=True)
    p = add("sim-build-apply", cmd_sim_build_apply, rnd=True)
    p.add_argument("--rounds", type=int, required=True)
    p = add("sim-critic-apply", cmd_sim_critic_apply, rnd=True)
    p.add_argument("--rounds", type=int, required=True)
    add("sim-publish", cmd_sim_publish, size=True)
    add("sim-finalize-apply", cmd_sim_finalize_apply)
    add("report-prep", cmd_report_prep)
    add("report-body-apply", cmd_report_body_apply)
    add("report-editor-apply", cmd_report_editor_apply)
    add("condense-apply", cmd_condense_apply, rnd=True)
    add("copyedit-prep", cmd_copyedit_prep)
    add("copyedit-apply", cmd_copyedit_apply)
    add("freeze", cmd_freeze)
    p = add("claims", cmd_claims)
    p.add_argument("--max", type=int, required=True)
    add("ledger", cmd_ledger)
    add("ledger-prep", cmd_ledger_prep)
    add("ledger-apply", cmd_ledger_apply)
    p = add("deliver", cmd_deliver)
    p.add_argument("--out", required=True, help="the deliverables folder")
    a = ap.parse_args(argv)
    try:
        return a.fn(a)
    except SystemExit as e:
        if isinstance(e.code, str):
            return fail(e.code)
        raise
    except Exception as e:  # one JSON line, never a bare traceback: the calling step reads only the last line
        return fail(f"{type(e).__name__}: {e}")


if __name__ == "__main__":
    sys.exit(main())
