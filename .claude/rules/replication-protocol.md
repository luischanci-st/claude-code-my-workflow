---
paths:
  - "Chkb/**/*.jl"
  - "Chkb/**/*.R"
  - "scripts/**/*.R"
---

# Estimation Pipeline Verification Protocol

**Core principle:** Verify the estimation pipeline reproduces known results before making changes.

---

## Phase 1: Baseline Verification

Before modifying any estimation code:

- [ ] Record current output values from `Chkb/Outcomes/Includes/`
- [ ] Store baseline in `quality_reports/` as reference:

```markdown
## Baseline Estimates

| Parameter | M1 | M2 | M3 |
|-----------|----|----|-----|
| ρ | [value] | [value] | [value] |
| σ_v | [value] | [value] | [value] |
| σ_u | [value] | [value] | N/A |
| δ_const | N/A | N/A | [value] |
```

---

## Phase 2: Run & Compare

- [ ] Run `julia Chkb/main.jl` (full pipeline)
- [ ] Compare new outputs against baseline
- [ ] Apply tolerance thresholds

### Tolerance Thresholds

| Type | Tolerance | Rationale |
|------|-----------|-----------|
| Point estimates (β, ρ, σ) | < 1e-6 | Numerical precision of Nelder-Mead optimization |
| Standard errors | < 1e-4 | Wild bootstrap MC variability with B=50 |
| Significance stars | Exact match | Discrete classification must agree |
| JLMS inefficiency | < 1e-5 | Computed through spatial filter |
| Bank count (N) | Exact match | 20-bank sample, no reason for difference |

### If Mismatch

**Do NOT proceed with modifications.** Isolate which step introduces the difference:
1. Step 1 OLS β coefficients changed?
2. Step 2 GMM grid search found different ρ?
3. Nelder-Mead converged to different σ_v/σ_u?
4. Bootstrap drew different Mammen weights?

Document the investigation even if unresolved.

---

## Phase 3: Julia-to-R Translation Pitfalls

| Julia | R (via RCall) | Trap |
|-------|---------------|------|
| `Float64` | `numeric` | R may silently convert to lower precision |
| `DataFrame` column types | R factor vs character | Check column types after `@rget` |
| `Matrix{Float64}` | Dense matrix | Sparse matrices lose sparsity through RCall |
| 1-based indexing | 1-based indexing | Same — but watch for 0-based C libraries |
| `missing` | `NA` | RCall handles this, but check edge cases |

---

## Phase 4: Verification Report

Save to `quality_reports/pipeline_verification.md`:

```markdown
# Pipeline Verification Report
**Date:** [YYYY-MM-DD]
**Julia version:** [version]
**R version:** [version]

## Summary
- **Parameters checked / Passed / Failed:** N / M / K
- **Overall:** [VERIFIED / PARTIAL / FAILED]

## Results Comparison

| Parameter | Baseline | Current | Diff | Status |
|-----------|----------|---------|------|--------|

## Environment
- Julia version, key packages
- R version, key packages
- OS, CPU cores used
```
