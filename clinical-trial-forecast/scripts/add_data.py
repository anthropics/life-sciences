#!/usr/bin/env python3
# Copyright 2026 Anthropic PBC
# SPDX-License-Identifier: Apache-2.0
"""Packs a trial's data into this skill's data/ folder, for a run that reads files only. Run it before the run, on
your own machine, then zip the skill folder and upload it as the skill for that one
question (the run reads data/ itself, so no file upload is needed). The packed files go wherever the uploaded skill
is shared, so keep confidential documents out of data/:

    python3 <skill folder>/scripts/add_data.py --out <new skill folder> --name <new skill name> --records <files or folders>
        [--docs <files or folders>] [--only <NCT ids>] [--zip out.zip]

--out copies this skill to a new folder and packs into the copy (without it, this skill's own data/ is used); --name
gives the copy its own skill name (lower case, digits and hyphens), so it does not clash with this skill when uploaded.

--records takes ClinicalTrials.gov API v2 JSON: single studies, search exports ({"studies": [...]}), JSON lists, or
one study per line (.jsonl / .ndjson); rows that keep a study's API v2 JSON as text in a raw_json field; and Parquet
files with a raw_json column (needs pyarrow). Each file is read whole before --only picks trials, so for a
large dump run this next to the data and keep --only short. The records go into data/records/records-NN.json, 500 to
a file. --docs takes papers, protocols and other documents (text, Markdown, HTML, CSV, TSV, PDF, Word .docx), copied to
data/docs/. Symbolic links and hidden files are left out.

The script keeps the packed skill to 200 files and 30 MB uncompressed, and checks both before it writes anything.
Prints what it packed."""

import argparse
import json
import re
import shutil
import sys
import zipfile
from pathlib import Path

sys.dont_write_bytecode = True  # no __pycache__ inside the skill folder, which is uploaded as it is
sys.path.insert(0, str(Path(__file__).resolve().parent))
from ft_registry import DOC_KINDS, NCT, records_in  # noqa: E402

SKILL = Path(__file__).resolve().parent.parent
MAX_FILES = 200               # this script's limits on the packed skill
MAX_BYTES = 30 * 1024 * 1024  # (uncompressed)


def files_of(paths):
    """The files to pack. Symbolic links are not followed (a link could pull in a file from outside the folder you chose),
    and hidden files are left out, as the run itself leaves them out."""
    for p in map(Path, paths):
        found = sorted(x for x in p.rglob("*") if x.is_file() and not any(s.startswith(".") for s in x.relative_to(p).parts)) if p.is_dir() else [p] if p.is_file() else None
        if found is None:
            print(f"Not found: {p}", file=sys.stderr)
            continue
        for x in found:
            if x.is_symlink():
                print(f"Skipping {x}: it is a symbolic link.", file=sys.stderr)
            else:
                yield x


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--records", nargs="*", default=[], help="ClinicalTrials.gov API v2 JSON files or folders")
    ap.add_argument("--docs", nargs="*", default=[], help="documents to copy into data/docs/")
    ap.add_argument("--only", nargs="*", default=[], help="keep only these NCT ids (or files listing them), e.g. to take a few trials from a large dump")
    ap.add_argument("--zip", default="", help="also write the packed skill as this zip")
    ap.add_argument("--replace", action="store_true", help="empty data/ first (it may hold another question's files)")
    ap.add_argument("--out", default="", help="copy the skill to this folder and pack into the copy")
    ap.add_argument("--name", default="", help="the copy's skill name (its folder name and the name in SKILL.md)")
    a = ap.parse_args()
    src_files = [p for p in SKILL.rglob("*") if p.is_file() and "__pycache__" not in p.parts
                 and not any(x.startswith(".") for x in p.relative_to(SKILL).parts) and p.relative_to(SKILL).parts[0] != "data"]
    if a.name and not re.fullmatch(r"[a-z0-9-]{1,64}", a.name):
        print("--name must be lower-case letters, digits and hyphens, at most 64 characters.", file=sys.stderr)
        return 1
    skill = SKILL
    if a.out:
        skill = Path(a.out)
        if a.name and skill.name != a.name:
            skill = skill / a.name
        if skill.exists():
            print(f"{skill} already exists. Choose a new --out folder or --name.", file=sys.stderr)
            return 1
        if SKILL == skill.resolve() or SKILL in skill.resolve().parents:
            print(f"--out {a.out} is inside this skill's folder, so the copy would copy itself. Choose a folder outside {SKILL}.", file=sys.stderr)
            return 1
    out = skill / "data"
    where = f"the copy in {skill}" if a.out else "this skill"
    def overlaps(p):  # p is data/, inside it, or a folder that holds it
        p, d = Path(p).resolve(), out.resolve()
        return p == d or d in p.parents or p in d.parents
    if not a.out and out.exists() and any(out.rglob("*")):
        if not a.replace:
            print(f"{out} already holds files, perhaps for another question, and they would join this run's evidence. "
                  "Pass --replace to empty it first.", file=sys.stderr)
            return 1
        clash = [p for p in a.records + a.docs if overlaps(p)]
        if clash:
            print(f"--replace would empty {out} first, but {clash[0]} is inside it or contains it. Give --records/--docs a "
                  f"folder outside {out} that does not contain it.", file=sys.stderr)
            return 1
    if a.zip and skill.resolve() in Path(a.zip).resolve().parents:
        print(f"--zip {a.zip} is inside the skill folder, so the zip would pack itself. Write it somewhere else.", file=sys.stderr)
        return 1

    # Read everything first; nothing is written until the limits below are checked.
    unusable = []
    packed = {}
    only = set()
    for x in a.only:
        only |= set(re.findall(r"NCT[0-9]{8}", Path(x).read_text(encoding="utf-8-sig").upper())) if Path(x).is_file() else {x.strip().upper()}
    for f in files_of(a.records):
        objs = []
        if f.suffix.lower() == ".parquet":
            try:
                import pyarrow.parquet as pq
            except ImportError:
                print(f"Skipping {f}: reading Parquet needs pyarrow.", file=sys.stderr)
                unusable.append(str(f))
                continue
            try:
                objs = [{"raw_json": x} for x in pq.read_table(f, columns=["raw_json"]).column("raw_json").to_pylist() if x]
            except Exception as e:
                print(f"Skipping {f}: could not read a raw_json column from it ({type(e).__name__}).", file=sys.stderr)
                unusable.append(str(f))
                continue
            text = ""
        else:
            text = f.read_text(encoding="utf-8-sig", errors="ignore")  # utf-8-sig: a byte-order mark (as Windows Notepad saves) is dropped
        try:
            if not text:
                raise ValueError
            objs = [json.loads(text)]
        except (ValueError, RecursionError):
            for line in text.splitlines():  # one study per line
                try:
                    objs.append(json.loads(line))
                except (ValueError, RecursionError):
                    pass
        try:
            got = [r for o in objs for r in records_in(o) if NCT.fullmatch(str(r.get("nct_id") or "")) and (not only or r["nct_id"] in only)]
        except Exception:  # one malformed file must not stop the packing
            got = []
        if not got:
            unusable.append(str(f))
        for r in got:
            packed[r["nct_id"]] = r  # a later file's record of the same trial wins
    if a.records and not packed:
        print("No registry records were found in --records" + (" for the --only ids" if only else "") + ". Nothing was written.",
              file=sys.stderr)
        return 1
    ids = sorted(packed)
    # 500 records to a file, to stay under MAX_FILES
    record_files = {f"records/records-{k // 500 + 1:02d}.json": json.dumps([packed[n] for n in ids[k:k + 500]], ensure_ascii=False, indent=1)
                    for k in range(0, len(ids), 500)}
    docs = {}
    for f in files_of(a.docs):
        if f.suffix.lower() not in DOC_KINDS:
            unusable.append(str(f))
            continue
        name, k = f.name, 2
        while name in docs and docs[name].read_bytes() != f.read_bytes():  # two different documents with one name
            name = f"{f.stem}_{k}{f.suffix}"
            k += 1
        docs.setdefault(name, f)
    n_files = len(src_files) + len(record_files) + len(docs)
    n_bytes = sum(p.stat().st_size for p in src_files) + sum(len(t.encode("utf-8")) for t in record_files.values()) + sum(p.stat().st_size for p in docs.values())
    if n_files > MAX_FILES or n_bytes > MAX_BYTES:
        what = (f"{n_files} files, and this script's limit is {MAX_FILES}" if n_files > MAX_FILES
                else f"{n_bytes / 2**20:.1f} MB, and this script's limit is {MAX_BYTES // 2**20} MB (uncompressed)")
        print(f"If packed as asked, {where} would hold {what}. Nothing was written. Pack fewer or smaller documents (or combine "
              "them), or fewer records with --only, and run this again.", file=sys.stderr)
        return 1

    if a.out:
        shutil.copytree(SKILL, skill, ignore=shutil.ignore_patterns("__pycache__", ".*", "data"))
        if a.name:
            # the copy's own name everywhere it names itself: the frontmatter and the launch text
            f = skill / "SKILL.md"
            f.write_text(f.read_text(encoding="utf-8").replace("clinical-trial-forecast", a.name), encoding="utf-8")
    elif out.exists():
        shutil.rmtree(out)
    for rel_path, text in record_files.items():
        (out / rel_path).parent.mkdir(parents=True, exist_ok=True)
        (out / rel_path).write_text(text, encoding="utf-8")
    for name, f in docs.items():
        (out / "docs").mkdir(parents=True, exist_ok=True)
        shutil.copy2(f, out / "docs" / name)
    plural = lambda n, word: f"{n} {word}" + ("" if n == 1 else "s")
    print(f"Packed {plural(len(packed), 'trial record')} and {plural(len(docs), 'document')} into {out} ({n_files} files in the skill).")
    if unusable:
        print("Not used (no registry record found, or not a document type the run reads): " + "; ".join(unusable[:20]))
    if a.zip:
        with zipfile.ZipFile(a.zip, "w", zipfile.ZIP_DEFLATED) as z:
            for p in sorted(skill.rglob("*")):
                if p.is_file() and "__pycache__" not in p.parts and not any(x.startswith(".") for x in p.relative_to(skill).parts):
                    # the entry's name as Unicode: under a C locale on Linux a non-ASCII file name arrives with escape characters
                    z.write(p, str(Path(skill.name) / p.relative_to(skill)).encode(sys.getfilesystemencoding(), "surrogateescape").decode("utf-8", "replace"))
        print(f"Wrote {a.zip}.")
    return 0

if __name__ == "__main__":
    sys.exit(main())
