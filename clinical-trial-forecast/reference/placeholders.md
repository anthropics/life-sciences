# How the prompts are filled

Helpers do not fill the placeholders in prompts/ by hand: scripts/ft_tools.py writes each step's whole prompt to a
file, with every {{placeholder}} filled from project.json and the files on disk, and the step's brief names that file.
Read the file and follow it.

Where a prompt's file layout differs from this run's:
- Paths in the prompts are relative to the project folder (your brief gives its full path). A path written data/...
  means the same path in the project folder, without data/.
- The "intake" a prompt mentions is intake/notes/background.md here, when the person supplied background.
- Each checklist worker group's answer goes into analyses/<NCT>/groups.json (see the card step's brief); the tools then
  write analyses/<NCT>/items/*.json and groups/*.md.
- Research answers are research/iterNN/qNN.json (with the number check beside them as qNN.check.json) and qNN.md.
- topics/<section>/SECTION.md is not written in this skill's runs; the integrator prompt says so itself ("may not
  exist").
- The documents the person supplied are in sources/docs/, listed in sources/INDEX.tsv. Any other document folder a
  prompt names does not exist in this run.
- Paper requests go in the paper_requests field as the prompts say; the run lists them in its closing note instead of
  sending them anywhere.
