# AGENTS.md

This file mirrors [CLAUDE.md](CLAUDE.md) for agents that read `AGENTS.md` by
convention (Codex, Cursor, Aider, etc.). The authoritative guidance is in
`CLAUDE.md` — read that. Key points:

- All physics is computed in nondimensional units; dimensional values are a
  UI-layer convenience.
- The model spec is in `ginzburg_landay_neutrophils.md`; the implementation
  plan is in `docs/PLAN.md`.
- Setup 3 (swarming) cue dynamics are not yet specified — do not invent them.
- No build step; plain ES modules + CDN libraries.
