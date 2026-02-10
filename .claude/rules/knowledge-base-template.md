---
paths:
  - "Edit/**/*.tex"
  - "Chkb/**/*.jl"
  - "Chkb/**/*.R"
  - "scripts/**/*.R"
---

# Research Knowledge Base: Spatial SFA for Banking

## Notation Registry

| Rule | Convention | Example | Anti-Pattern |
|------|-----------|---------|-------------|
| Spatial autocorrelation | ρ (rho) | ρ ∈ (-0.6, 0.03) | Do not use λ for spatial lag |
| Weight matrix | W_t (time-varying) | Row-normalized directed adjacency | Using W without time subscript |
| Noise variance | σ_v | σ_v > 0 | Confusing σ_v with σ_u |
| Inefficiency variance | σ_u (M1/M2) or δ vector (M3) | σ_u > 0 | Negative inefficiency variance |
| Cost function coefficients | β | Within-estimator (Step 1 OLS) | Using pooled OLS |
| Spatial filter | (I - ρW_t)^{-1} | Pre-computed per grid point | Forgetting to invert |

## Symbol Reference

| Symbol | Meaning | Where Used |
|--------|---------|------------|
| W_t | Row-normalized directed adjacency matrix from interbank exposures | Step 2 GMM |
| ρ | Spatial autocorrelation parameter | Grid search in Step 2 |
| σ_v | Noise standard deviation (half-normal) | Step 2 GMM optimization |
| σ_u | Inefficiency standard deviation (M1/M2) | Step 2 GMM optimization |
| δ | Inefficiency determinants vector (M3) | Heterogeneous inefficiency |
| β | Translog cost function coefficients | Step 1 within-estimator |
| JLMS | Jondrow-Lovell-Materov-Schmidt inefficiency estimator | Post-estimation |
| B | Bootstrap replications (default 50) | Wild bootstrap |

## Model Specifications

| Model | Type | Key Parameters | Description |
|-------|------|---------------|-------------|
| M1 | `:basic` | ρ, σ_v, σ_u | Cost frontier with spatial dependence |
| M2 | `:basic` + risk | ρ, σ_v, σ_u | Adds capital adequacy ratio |
| M3 | `:determinants_v1` | ρ, σ_v, δ vector | Ownership and segment dummies replace σ_u |

## Data Pipeline

| Stage | Input | Output | Tool |
|-------|-------|--------|------|
| Raw data | CSV/XLSX in `Data/Raw/` | Panel `df.Rda` + `all_W_df_*.Rda` | R (`dataprocessing.R`) |
| Step 1 OLS | `df.Rda` | β coefficients, residuals | Julia (`Chkb.jl`) |
| Step 2 GMM | Residuals + `all_W_df_Tot.Rda` | ρ, σ_v, σ_u/δ | Julia (`Chkb.jl`) |
| Post-estimation | Step 2 results | JLMS inefficiency | Julia (`Chkb.jl`) |
| Wild bootstrap | Full pipeline | SEs, p-values | Julia (`Chkb.jl`, parallelized) |
| Output | All estimates | `.txt` + `.rds` in `Outcomes/Includes/` | Julia (`Chkb.jl`) |
| Paper | `.txt` files | `performance.tex` via `\input{}` | LaTeX (`pdflatex`) |

## Network Variants

| Variant | Variable | Description |
|---------|----------|-------------|
| Tot | `all_W_df_Tot` | Total interbank obligations (default) |
| FIN | `all_W_df_FIN` | Financial derivatives exposure |
| CAP | `all_W_df_CAP` | Deposits/capital exposure |
| UNSEC | `all_W_df_UNSEC` | Unsecured interbank exposure |

## R/Julia Code Pitfalls

| Bug | Impact | Fix |
|-----|--------|-----|
| RCall data type mismatch | Julia Float64 vs R numeric silently rounds | Explicit type conversion at boundary |
| W matrix not row-normalized | Biased ρ estimates | Verify row sums ≈ 1 after construction |
| igraph directed vs undirected | Wrong adjacency structure | Use `graph_from_data_frame(directed=TRUE)` |
| Nelder-Mead non-convergence | Silent bad optimum | Check convergence flag; try multiple starting values |
| CPI deflation timing | Nominal vs real CLP confusion | Match CPI period to observation date |
| `@everywhere` missing | Worker processes lack definitions | Annotate all shared functions/data |

## Anti-Patterns (Don't Do This)

| Anti-Pattern | What Happened | Correction |
|-------------|---------------|-----------|
| Symmetric W for directed network | Interbank exposures are directional | Use directed adjacency from C18 data |
| Pooled OLS for Step 1 | Ignores bank fixed effects | Use within-group transformation |
| Single starting value for Nelder-Mead | Stuck at local minimum | Grid search over ρ, then optimize σ_v/σ_u |
