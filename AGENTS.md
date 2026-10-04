# AGENTS.md

This file mirrors [CLAUDE.md](CLAUDE.md) for agents that read `AGENTS.md` by
convention (Codex, Cursor, Aider, etc.). The authoritative guidance is in
`CLAUDE.md` — read that. The repo map is in [README.md](README.md). Key points:

- All physics is computed in nondimensional units; dimensional values are a
  UI-layer convenience.
- The model spec is in `ginzburg_landay_neutrophils.md`; per-setup physics in
  `docs/physics/`; the implementation plan is in `docs/PLAN.md`.
- Docs lead code: spec → `docs/physics/` → `docs/PLAN.md` → source.
- No build step; plain ES modules + CDN libraries. Serve with `python3 serve.py`.
