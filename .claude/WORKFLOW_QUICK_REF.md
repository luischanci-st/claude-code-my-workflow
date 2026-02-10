# Workflow Quick Reference

**Model:** Contractor (you direct, Claude orchestrates)

---

## The Loop

```
Your instruction
    ↓
[PLAN] (if multi-file or unclear) → Show plan → Your approval
    ↓
[EXECUTE] Implement, verify, done
    ↓
[REPORT] Summary + what's ready
    ↓
Repeat
```

---

## I Ask You When

- **Design forks:** "Option A (fast) vs. Option B (robust). Which?"
- **Code ambiguity:** "Spec unclear on X. Assume Y?"
- **Replication edge case:** "Just missed tolerance. Investigate?"
- **Scope question:** "Also refactor Y while here, or focus on X?"

---

## I Just Execute When

- Code fix is obvious (bug, pattern application)
- Verification (tolerance checks, compilation, Julia/R runs)
- Documentation (logs, commits)
- Plotting (per established standards)

---

## Quality Gates (No Exceptions)

| Score | Action |
|-------|--------|
| >= 80 | Ready to commit |
| < 80  | Fix blocking issues |

---

## Non-Negotiables

- **Path convention:** Relative paths from `banking1/`; `here::here()` for R scripts
- **Seed convention:** Not applicable — GMM is deterministic; wild bootstrap uses Mammen two-point weights
- **Figure standards:** White background, 300 DPI, PDF format for paper, publication-quality
- **Tolerance thresholds:** 1e-6 for point estimates (β, ρ, σ); 1e-4 for SEs (bootstrap variability)

---

## Preferences

**Visual:** PDF figures, white background, journal-standard dimensions (6.5" x 4.5")
**Reporting:** Concise bullets; flag near-misses in replication
**Session logs:** Always (post-plan, incremental, end-of-session)
**Replication:** Strict — flag any deviation; exact match for significance stars

---

## Exploration Mode

For experimental work, use the **Fast-Track** workflow:
- Work in `explorations/` folder
- 60/100 quality threshold (vs. 80/100 for production)
- No plan needed — just a research value check (2 min)
- See `.claude/rules/exploration-fast-track.md`

---

## Next Step

You provide task → I plan (if needed) → Your approval → Execute → Done.
