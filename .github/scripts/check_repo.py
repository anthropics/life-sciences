#!/usr/bin/env python3
"""Check that the marketplace, skills, Markdown and data files in this repository are valid."""

import json
import re
import sys
from pathlib import Path

import yaml

SKIP_DIRS = {".git", "node_modules", "__pycache__", ".venv"}
FENCE = re.compile(r"^ {0,3}(`{3,}|~{3,})(.*)$")
LINK = re.compile(r"\[[^\]]*\]\(([^)\s]+)(?:\s+\"[^\"]*\")?\)")
INLINE_CODE = re.compile(r"`[^`]*`")
SKILL_NAME = re.compile(r"[a-z0-9]+(-[a-z0-9]+)*")


def repo_files(root, *suffixes):
    for path in sorted(root.rglob("*")):
        if path.is_file() and path.suffix in suffixes and not SKIP_DIRS & set(path.relative_to(root).parts):
            yield path


def split_fences(text):
    """Return the lines outside code blocks, and the line number of a block left open, if any."""
    outside, opened, marker = [], None, None
    for number, line in enumerate(text.splitlines(), 1):
        match = FENCE.match(line)
        if marker is None:
            if match:
                opened, marker = number, match.group(1)
            else:
                outside.append(line)
        elif match and match.group(1)[0] == marker[0] and len(match.group(1)) >= len(marker) and not match.group(2).strip():
            opened, marker = None, None
    return outside, opened


def check_markdown(root):
    problems = []
    for path in repo_files(root, ".md"):
        name = path.relative_to(root)
        outside, opened = split_fences(path.read_text(encoding="utf-8"))
        if opened:
            problems.append(f"{name}:{opened}: code block is never closed")
        for line in outside:
            for target in LINK.findall(INLINE_CODE.sub("", line)):
                target = target.split("#")[0]
                is_path = "/" in target or "." in target  # a bare word such as (url) is a format example
                if is_path and not re.match(r"[a-z][a-z0-9+.-]*:", target) and not (path.parent / target).exists():
                    problems.append(f"{name}: broken link to {target}")
    return problems


def check_data_files(root):
    problems = []
    for path in repo_files(root, ".json"):
        try:
            json.loads(path.read_text(encoding="utf-8"))
        except ValueError as error:
            problems.append(f"{path.relative_to(root)}: invalid JSON: {error}")
    for path in repo_files(root, ".yaml", ".yml"):
        try:
            list(yaml.safe_load_all(path.read_text(encoding="utf-8")))
        except yaml.YAMLError as error:
            problems.append(f"{path.relative_to(root)}: invalid YAML: {error}")
    return problems


def check_skill(root, folder):
    skill = root / folder / "SKILL.md"
    if not skill.is_file():
        return [f"{folder}: listed in the marketplace but has no SKILL.md"]
    match = re.match(r"---\n(.*?)\n---\n", skill.read_text(encoding="utf-8"), re.S)
    if not match:
        return [f"{folder}/SKILL.md: missing frontmatter"]
    try:
        fields = yaml.safe_load(match.group(1))
    except yaml.YAMLError as error:
        return [f"{folder}/SKILL.md: invalid frontmatter: {error}"]
    problems = []
    if not SKILL_NAME.fullmatch(str(fields.get("name", ""))):
        problems.append(f"{folder}/SKILL.md: name must be lower-case letters, digits and hyphens")
    if not fields.get("description"):
        problems.append(f"{folder}/SKILL.md: description is missing")
    if not (root / folder / "LICENSE.txt").is_file():
        problems.append(f"{folder}: LICENSE.txt is missing")
    return problems


def check_marketplace(root):
    marketplace = json.loads((root / ".claude-plugin" / "marketplace.json").read_text(encoding="utf-8"))
    listed = {
        str(Path(skill)) for plugin in marketplace["plugins"] for skill in plugin.get("skills", [])
    }
    problems = []
    for folder in sorted(listed):
        problems += check_skill(root, folder)
    for skill in sorted(root.glob("*/SKILL.md")):
        if skill.parent.name not in listed:
            problems.append(f"{skill.parent.name}: has a SKILL.md but is not in the marketplace")
    return problems


def main(root="."):
    root = Path(root).resolve()
    problems = check_marketplace(root) + check_markdown(root) + check_data_files(root)
    for problem in problems:
        print(problem)
    print(f"{len(problems)} problem(s) found")
    return 1 if problems else 0


if __name__ == "__main__":
    sys.exit(main(*sys.argv[1:]))
