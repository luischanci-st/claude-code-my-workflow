# CLAUDE.MD -- Research Paper: Interconnectedness and Banking Performance

**Project:** Interconnectedness and Banking Performance
**Authors:** Chanci, Kumbhakar & Bobadilla
**Institution:** Corporaci&oacute;n Santo Tom&aacute;s / co-authors' affiliations
**Branch:** main

---

## Core Principles

- **Plan first** -- enter plan mode before non-trivial tasks; save plans to `quality_reports/plans/`
- **Verify after** -- compile/run and confirm output at the end of every task
- **Quality gates** -- nothing ships below 80/100
- **[LEARN] tags** -- when corrected, save `[LEARN:category] wrong → right` to MEMORY.md

---

## Folder Structure

```
banking1/
├── CLAUDE.MD                    # This file
├── .claude/                     # Rules, skills, agents, hooks
├── Chkb/                        # Julia/R estimation package
│   ├── main.jl                  # Pipeline orchestrator (entry point)
│   ├── Chkb.jl                  # Core module: all estimation algorithms
│   ├── dataprocessing.R         # Raw CSV → panel DataFrame + adjacency matrices
│   ├── Data/Raw/                # Input: CSV + XLSX files
│   ├── Data/Processed/          # Generated .Rda files (auto-created)
│   └── Outcomes/Includes/       # Output: .txt (LaTeX \input) + .rds DataFrames
├── Edit/                        # Paper LaTeX source
│   ├── performance.tex          # Main manuscript
│   ├── banking.bib              # Bibliography
│   ├── figures/                  # Paper figures (PDF)
│   └── includes/                # LaTeX-ready parameter estimates
├── Figures/                     # Additional figures
├── scripts/                     # Utility scripts (Python, R)
├── quality_reports/             # Plans, session logs, merge reports
├── explorations/                # Research sandbox (see rules)
├── templates/                   # Session log, quality report templates
└── master_supporting_docs/      # Supporting papers and materials
```

---

## Commands

```bash
# Julia estimation pipeline (from repo root)
julia Chkb/main.jl

# R data processing only
Rscript -e "setwd('Chkb'); source('dataprocessing.R')"

# LaTeX paper compilation (3-pass, pdflatex)
cd Edit && pdflatex -interaction=nonstopmode performance.tex
bibtex performance
pdflatex -interaction=nonstopmode performance.tex
pdflatex -interaction=nonstopmode performance.tex

# Quality score
python scripts/quality_score.py Edit/performance.tex
```

---

## Model Specifications

| Model | Type | Estimates | Description |
|-------|------|-----------|-------------|
| M1 | `:basic` | ρ, σ_v, σ_u | Cost frontier with spatial dependence |
| M2 | `:basic` + risk | ρ, σ_v, σ_u | Adds capital adequacy ratio as risk regressor |
| M3 | `:determinants_v1` | ρ, σ_v, δ vector | Heterogeneous inefficiency determinants (ownership, segment) |

---

## Quality Thresholds

| Score | Gate | Meaning |
|-------|------|---------|
| 80 | Commit | Good enough to save |
| 90 | PR | Ready for review |
| 95 | Excellence | Aspirational |

---

## Skills Quick Reference

| Command | What It Does |
|---------|-------------|
| `/compile-latex [file]` | 3-pass pdflatex + bibtex |
| `/proofread [file]` | Grammar/typo/consistency review |
| `/review-r [file]` | R code quality review |
| `/review-paper [file]` | Manuscript review |
| `/validate-bib` | Cross-reference citations |
| `/devils-advocate` | Challenge design decisions |
| `/commit [msg]` | Stage, commit, PR, merge |
| `/lit-review [topic]` | Literature search + synthesis |
| `/research-ideation [topic]` | Research questions + strategies |
| `/interview-me [topic]` | Interactive research interview |
| `/data-analysis [dataset]` | End-to-end R analysis |

---

## Key Output Files

| File Pattern | Location | Content |
|-------------|----------|---------|
| `M{1,2,3}_{param}.txt` | `Chkb/Outcomes/Includes/` | Point estimates for LaTeX `\input{}` |
| `M{1,2,3}_{param}_se.txt` | `Chkb/Outcomes/Includes/` | Bootstrap standard errors |
| `M{1,2,3}_{param}_pval.txt` | `Chkb/Outcomes/Includes/` | Significance stars (`***`/`**`/`*`) |
| `M{1,2,3}_df_out.rds` | `Chkb/Outcomes/Includes/` | Full result DataFrames |
| `performance.tex` | `Edit/` | Main manuscript |
| `banking.bib` | `Edit/` | Bibliography |

---

## Current Project State

| Component | Status | Notes |
|-----------|--------|-------|
| Data processing (R) | Complete | `dataprocessing.R` generates all `.Rda` files |
| M1 estimation | Complete | Basic spatial cost frontier |
| M2 estimation | Complete | With risk regressor |
| M3 estimation | Complete | With inefficiency determinants |
| Paper manuscript | In progress | `Edit/performance.tex` |
