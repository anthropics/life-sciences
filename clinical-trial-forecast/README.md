# Clinical Trial Forecast

Estimate the probability that a registered clinical trial meets its goal. You give Claude the trial's NCT id, your question and a folder of source files. Claude returns a probability of success, a cited report and the simulation model behind the number. The skill runs in [Claude Science](https://claude.com/docs/claude-science/overview).

## How it works

Separate sub-agents do each step, so no agent checks its own work.

| Phase | What happens |
|---|---|
| 1. Registry | Ranks related trials from your records. |
| 2. Trial cards | Writes a fact card for each key trial from a 24-item checklist. A second agent checks the numbers. |
| 3. Research | Answers questions from your files and re-derives the key numbers. |
| 4. Model | Builds a Monte Carlo simulation, which a critic reviews. |
| 5. Report | Writes the report, then checks and corrects each headline claim. |

## What you need

- **Claude Science on Linux.** macOS is untested. Windows is not supported.
- **Delegation** turned on in the session settings menu.
- **Python 3 with numpy, scipy and pandas.** See [Dependencies](#dependencies).
- **Your own source files:**
  - ClinicalTrials.gov records as API v2 JSON: the trial's own record (required), and records of similar trials and of other drugs in the class.
  - Papers, protocols, analysis plans and regulatory reviews as text, Markdown or HTML. The run accepts PDF and Word files but cannot search them.

## Install

**For yourself.** In Claude Science, go to Settings > Skills > Add skill > Upload a skill. Choose the `clinical-trial-forecast-<version>.zip` file from this repository's releases.

**For your organization.** An admin uploads the same zip in claude.ai under Organization settings > Plugins & skills > Upload a skill. The Skills switch on the Policy tab must be on.

## Run it

**1. Make a run folder and add your files.**

```
mkdir -p ~/trial-forecast/uploads/records ~/trial-forecast/uploads/docs ~/trial-forecast/outputs
```

Copy records into `uploads/records/` and documents into `uploads/docs/`. To download one record:

```
curl -o ~/trial-forecast/uploads/records/NCT01234567.json "https://clinicaltrials.gov/api/v2/studies/NCT01234567?format=json"
```

Copy the files in: the run ignores symbolic links and hidden files. The folder's full path may hold only letters, digits and `. _ - /`. Use a new run folder for each trial, because a run resumes from whatever it finds in `outputs/`.

**2. Make a project** in Claude Science and paste this into its custom instructions:

```
This project runs the clinical-trial-forecast skill.
If your first message begins "HOW THIS STEP IS RUN HERE", you are a helper (a sub-agent): do what that message says and skip the rest of this note.
- Every helper the skill's workflow script hands you is done by a sub-agent, never by you. Treat the loop in SKILL.md part B, step 6, as the task itself.
- Sub-agents reach files by absolute path only.
- Stop only processes you started, by their PID. Never use pkill -f or killall.
- Do not ask me questions. Go straight through the skill's phases.
```

**3. Start a session with Delegation on and send this request**, filled in:

```
Run the clinical-trial-forecast skill: read its SKILL.md, part B, and follow it.

RUN SETTINGS
Run folder: <the full path of your run folder, for example /home/you/trial-forecast>
Size: test
Files: in uploads/: records/ holds <how many> ClinicalTrials.gov records (the trial's own among them); docs/ holds <how many and what kind of> documents
Data access allowed: files only
Concurrency: at most 8 helpers at once, 16 at full size

TRIAL
<the NCT id>

QUESTION (the person's request, verbatim)
Goal: <what you want estimated>
What counts as success: <for example, the trial meets its primary endpoint as pre-specified>
Decision the forecast informs: <what you will do with the number>
Audience: <who reads the report>

BACKGROUND (confidential, supplied by the person; optional)
none
```

When the app asks for the run folder, choose Read & write. Stay nearby for the first few minutes, and refuse any other folder or website you did not expect. A files-only run needs neither.

## What you get back

The run writes these files to `outputs/deliverables/`:

| File | What it holds |
|---|---|
| `forecast_summary.json` | The probability of success, and its Monte Carlo error (how much the number moves from simulation noise alone) |
| `forecast_report.md`, `.html` | The report, with numbered sources |
| `forecast_claim_ledger.md` | Each headline claim and its checker's verdict |
| `model/` | The model's inputs, the code an agent wrote, its results and the critiques |

To skip the report, add `Stop after phase: 4` under RUN SETTINGS. The forecast is then `pts_overall` (probability of technical success) in `outputs/work/project/sim/results.json`.

## What to expect

Start with `test` to check your setup. Larger sizes take longer. A run can take from about an hour to most of a day.

| Size | Trial cards | Research questions | Model rounds |
|---|---|---|---|
| `test` | 3 | 2 | 1 |
| `small` (default) | 8 | up to 8 | up to 2 |
| `medium` | 15 | up to 48 | up to 3 |
| `full` | 25 | up to 128 | up to 3 |

Two runs on the same files give different numbers. For a decision that matters, run the skill more than once and compare. The closing message lists papers the run wanted and could not find.

## Safety and privacy

- **Agents write and run code** on your computer, in the Claude Science sandbox. The sandbox and the folders you grant limit what that code can do. Grant your run folder only.
- **Text in your files can steer the agents.** Use sources you trust, and check claims against their sources. A `project.json` at the top of `uploads/` replaces the run's setup step.
- **Some limits are instructions only.** "Files only" is an instruction to the agents. So is keeping your question and background away from the steps that search the web. Nothing enforces either. To be sure, turn off research connectors in Settings > Connectors, and leave BACKGROUND as none if it is sensitive.
- **Where your data goes.** The skill's scripts write plain, unencrypted files to your run folder. What the agents read is sent to Anthropic ([how Claude Science works with your data](https://claude.com/docs/claude-science/how-claude-science-works-with-your-data)). A web search or connector you allow receives the agents' queries.
- **Background jobs.** When the workflow finishes, it stops its simulation jobs. If a session ends early, a job can keep running. Find it with `ps axww | grep "bg.py _run"` and stop it with `kill -- -PID`, which also stops the simulation it started.

## Dependencies

The skill installs nothing. Its scripts use only the Python 3 standard library. The model step needs numpy, scipy and pandas. Three packages are optional ([requirements.txt](requirements.txt)): matplotlib for figures, markdown for the HTML report, and pyarrow for Parquet input to `scripts/add_data.py`.

## License

Apache 2.0. See [LICENSE.txt](LICENSE.txt).
