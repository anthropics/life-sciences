# Contributing

Thanks for your interest in improving the Life Sciences marketplace. This guide covers how to propose changes and what a pull request needs before it can be reviewed.

## Before you start

- **Fixes and small improvements**: open a pull request directly.
- **New plugins or skills**: open an issue first to describe what you'd like to add and who it's for. We can't accept every proposal, and this saves you building something that won't be merged.
- **Security issues**: don't open a public issue. See [SECURITY.md](SECURITY.md).

## Adding a skill

Each skill lives in its own top-level directory:

```
your-skill/
├── SKILL.md          # Required: frontmatter with name and description, then instructions
├── LICENSE.txt       # Required: Apache-2.0, as used by the existing skills
├── references/       # Optional: detail loaded when needed
└── scripts/          # Optional: helper code
```

Then:

1. Add an entry to `.claude-plugin/marketplace.json`, following the existing skill entries (`"source": "./"`, `"strict": false`, and a `skills` list pointing at your directory).
2. Add the install command to the Quick Start list in `README.md`.

A good skill:

- Has a `description` that says what it does and when to use it.
- Does one job. If it depends on another skill, agent or connector, says so and explains what happens when it's missing.
- Keeps file references relative to the skill directory.
- Cites its sources, and doesn't tell Claude to guess at scientific or clinical facts.

## Adding an MCP server

Please open an issue first, so we can talk through the best way to make it available.

## Before you open a pull request

- Run `claude plugin validate .` from the repository root.
- Test from your local checkout: run `/plugin marketplace add ./path/to/life-sciences`, then `/plugin install your-skill@life-sciences`, and try it on a realistic task.
- Check that you have the right to contribute everything in the change, and that any third-party material is credited and licensed.
- Remove anything that shouldn't be public: credentials, personal or patient data, confidential or unpublished data, internal paths and hostnames, and names of real people or companies that aren't needed.

Everything pushed to this repository is public as soon as it's pushed, including branches and draft pull requests, and removing it later doesn't undo that.

## Review

Pull requests need an approving review before they can be merged. Keep each pull request to one plugin or skill where you can, so it's easy to review and revert.
