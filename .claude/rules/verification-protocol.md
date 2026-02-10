---
paths:
  - "Edit/**/*.tex"
  - "Chkb/**/*.jl"
  - "Chkb/**/*.R"
  - "scripts/**/*.R"
---

# Task Completion Verification Protocol

**At the end of EVERY task, Claude MUST verify the output works correctly.** This is non-negotiable.

## For LaTeX Paper (Edit/performance.tex):
1. Compile with pdflatex (3-pass + bibtex):
   ```bash
   cd Edit && pdflatex -interaction=nonstopmode performance.tex
   bibtex performance
   pdflatex -interaction=nonstopmode performance.tex
   pdflatex -interaction=nonstopmode performance.tex
   ```
2. Check for overfull hbox warnings
3. Check for undefined citations or missing references
4. Verify PDF was generated with non-zero size

## For Julia Scripts (Chkb/*.jl):
1. Verify module loads without errors:
   ```bash
   julia -e "include(\"Chkb/Chkb.jl\"); using .Chkb"
   ```
2. For full pipeline: run `julia Chkb/main.jl` (caution: long-running)
3. Check that `Chkb/Outcomes/Includes/` contains expected output files
4. Verify `.txt` files contain numeric values (not NaN or errors)
5. Verify `.rds` files have non-zero size

## For R Scripts (Chkb/dataprocessing.R):
1. Run `Rscript Chkb/dataprocessing.R`
2. Verify `.Rda` files created in `Chkb/Data/Processed/`
3. Spot-check: panel dimensions, bank count, time period coverage

## For Bibliography (Edit/banking.bib):
- Check that all `\cite` references in `performance.tex` have entries in `banking.bib`
- Verify no duplicate keys

## Common Pitfalls:
- **RCall interop**: Julia's RCall may fail if R libraries aren't installed system-wide
- **Path separators**: Windows uses `\` but Julia/R often need `/` — use `joinpath()` in Julia
- **Spatial weight matrices**: Must be row-normalized; check `all_W_df_Tot.Rda` dimensions match panel
- **Bootstrap timing**: Wild bootstrap with B=50 takes significant time — don't run unnecessarily

## Verification Checklist:
```
[ ] Output file created successfully
[ ] No compilation/execution errors
[ ] Figures/tables display correctly
[ ] Parameter estimates within expected ranges
[ ] Reported results to user
```
