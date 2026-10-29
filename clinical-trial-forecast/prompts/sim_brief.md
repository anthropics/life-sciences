# Simulation brief

Model: {{target_description}}.
Technical success: {{success_definition}}.
Project: {{indication}} · asset: {{asset}} · report goal: {{goal}} · decision: {{decision}}.

The parameter sheet `params.json` in this directory is authoritative. Implement it literally; where it is ambiguous, choose the simplest faithful interpretation and record the choice in `run_log.txt` and in the report's "Implementation choices" table.

Deliverables (all in this directory):
- `simulate.py` — numpy/scipy/stdlib only; `python3 simulate.py --n-sims {{n_sims}} --seed 1 --out results.json` must run end to end in under 60 minutes, started in the background with bg.py (your bash tool stops any single command after about 5 minutes; your step's rules say how); split a longer run into pieces, each writing its own file; `--scenario <id>` and `--set <param_id>=<value>` overrides; every function short and named for what it does in the trial (sample_epistemic, simulate_trial, analyze_trial, program_success, …).
- `results.json` — see sim_build instructions for the required keys.
- `figures/*.png` — matplotlib, one idea per figure, axis labels with units, titles that state the takeaway.
- `sim_report.md` — the simulation chapter of the report (plain language first, numbers second, code last).
- `run_log.txt` — what you ran, in order, with outputs; implementation choices; the two-seed reproducibility check.
