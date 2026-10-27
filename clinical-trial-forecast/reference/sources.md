# Sources: what each step may use, and how

The person decided at launch what this run may use. There are two kinds of source, and a run may have either or both:
- **The run's files** (attached to the task, or packed into the skill's data/ folder). Registry records are under
  raw/ctgov/ and documents (papers, protocols, regulatory documents) under sources/docs/, listed in sources/INDEX.tsv.
  Every step may use them, and a run on files only uses nothing else. Your prompt starts with a note on the local source library, and, when the run has no outside
  source, with a note saying that you have the project files only: then use the files wherever the prompt says to
  search the web or query a connector, and write what is missing as a gap.
- **Outside tools the person allowed**, quoted in your brief as "Data access the person allowed": a registry or
  literature connector of their organization, or the task's web search and web fetch, only where the person said
  their organization permits them for this work. Use exactly those and no other, even if you have more tools.

A search query is text sent outside, so never put anything from the person's background (intake/) into a query or a
fetched URL, whatever your step is. The trial being forecast and its drug are public, and may be named freely, unless
project.json's asset.public is false; then follow the confidentiality rule your prompt states.

The prompts in prompts/ name particular tools: connectors called clinical_trials and pubmed (and biorxiv, chembl and
fda), and tools called WebSearch and web_fetch. Your task may have none of these. Map each name onto what you are
allowed to use:
- "clinical_trials MCP" or "clinical_trials connector": an allowed registry connector (its tool names may differ);
  otherwise the registry's API v2 (below) through web_fetch, if allowed; otherwise the run's records (from its files,
  or built from web search as described below).
- "pubmed MCP" or "pubmed connector": an allowed literature connector; otherwise web_search / web_fetch of PubMed and
  journal pages, if allowed; otherwise the attached documents.
- "WebSearch": web_search, if allowed. "{{fetch_tools}}" or web_fetch: web_fetch, if allowed.
- biorxiv, chembl, fda connectors: only if allowed.
Read, Grep and Glob in the prompts mean your own file tools (read a file, search file contents, list files by
pattern), or bash with cat, grep and ls if you have no such tools. If a prompt names any other tool you do not have or
may not use, skip it. Never fill a gap from memory as if it came from a source: say what is missing in the gaps field
and, when a specific paper would settle it, in paper_requests.

Never reach websites from code (curl, python, or any other command): outside sources come only through the tools the
person allowed.

## Registry searches (the search step; only in runs with an outside registry route)

registry/queries.json lists each group's queries. Each query has params (condition or intervention, status, phase,
study_type, advanced_query in the registry's AREA[...] syntax) and a ready-made api_url.
1. Connector: call its search tool with the same filters (condition or intervention, the statuses, the phases, the
   date range from advanced_query), asking for as many results as the group keeps (its cap). If it cannot apply a
   filter, say so in the query's note.
2. web_fetch of api_url: it asks only for NCT ids (fields=NCTId), at most as many as the group keeps (pageSize), and
   the total (countTotal=true). The answer is JSON: studies[].protocolSection.identificationModule.nctId and
   totalCount. One page is enough: the run keeps at most that many ids per group, so do not follow nextPageToken.
3. web_search alone cannot list a registry search completely. Use it only as the next section says, when your brief
   tells you to build the list from web search (a run with no registry route).
Record only ids that appeared in a tool result. Put the source you used in "source" and the registry's total in "total".

## Registry records (only in runs with an outside registry route; otherwise the attached records are used)

1. Connector: its get-trial-details tool for the NCT id.
2. web_fetch of the trial's record_url from the batch file: the registry's API v2 record, limited to the fields the
   record shape needs. If that address returns an error, fetch https://clinicaltrials.gov/api/v2/studies/<NCT>?format=json
   (the whole record) instead, and say so in your step's problem field.
Write raw/ctgov/<NCT>.json in the trial_record shape of reference/schemas.json, with fetched_at from the command
date -u +%Y-%m-%dT%H:%M:%SZ run right after the fetch. Copy values exactly. Dates as the
registry gives them (YYYY-MM-DD or YYYY-MM). Status in the registry's upper-case form (RECRUITING, COMPLETED, ...).
Phases as a list (["PHASE2"]). Keep every site's facility, city and country for trials that are recruiting, not yet
recruiting, enrolling by invitation or active, because the competitive-landscape tables come from them. You may save a
fetched API v2 study object unchanged as the file's "trial" value instead; the tools convert it.

## Records and searches built from web search (only in a run with no registry route)

When the run has no registry connector and no web fetch but may use web search, the registry facts come from web
search results instead, as well as they can. Everything built this way is marked, and the report says so.
- A trial's record: search for its NCT id, and its acronym or title once you know them. Fill only the fields a search
  result states (the registry's own page text when a result shows it, then the sponsor's releases and publications),
  copying values as given; leave every other field null or []. Write the record in the trial_record shape with
  "source": "web search (not the registry)", "fetched_at" as usual, and "field_sources": for each field you filled,
  the URL(s) that state it. Write "found": false only if no result shows that the trial exists.
- A group's search: run a few searches for the query's condition or drug with its status and phase (for example
  "<condition> phase 3 completed trial NCT"), and list only NCT ids that appear in the results. Use "source": "web
  search", "total": null, and a note saying the list is incomplete.

## Literature and other documents

Cite every fact inline as the prompts say: [NCT01234567], [PMID 12345678], [DOI 10.xxxx/yyyy] or [URL https://...].
Registry facts cite the NCT id. An attached document with no identifier of its own is cited by its path,
[URL sources/docs/<name>]. A fact with no source you can cite does not go in.
