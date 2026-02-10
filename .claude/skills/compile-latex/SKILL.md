---
name: compile-latex
description: Compile the research paper with pdflatex (3 passes + bibtex). Use when compiling the manuscript.
disable-model-invocation: true
argument-hint: "[filename without .tex extension, default: performance]"
allowed-tools: ["Read", "Bash", "Glob"]
---

# Compile Research Paper LaTeX

Compile the manuscript using pdflatex with full citation resolution.

## Steps

1. **Navigate to Edit/ directory** and compile with 3-pass sequence:

```bash
cd Edit
pdflatex -interaction=nonstopmode $ARGUMENTS.tex
bibtex $ARGUMENTS
pdflatex -interaction=nonstopmode $ARGUMENTS.tex
pdflatex -interaction=nonstopmode $ARGUMENTS.tex
```

**Alternative (latexmk):**
```bash
cd Edit
latexmk -pdf -interaction=nonstopmode $ARGUMENTS.tex
```

2. **Check for warnings:**
   - Grep output for `Overfull \\hbox` warnings
   - Grep for `undefined citations` or `Label(s) may have changed`
   - Report any issues found

3. **Report results:**
   - Compilation success/failure
   - Number of overfull hbox warnings
   - Any undefined citations
   - PDF page count

## Why 3 passes?
1. First pdflatex: Creates `.aux` file with citation keys
2. bibtex: Reads `.aux`, generates `.bbl` with formatted references
3. Second pdflatex: Incorporates bibliography
4. Third pdflatex: Resolves all cross-references with final page numbers

## Important
- **Use pdflatex**, not xelatex (paper uses standard fonts)
- The `.bib` file (`banking.bib`) lives in `Edit/` alongside the `.tex` file
- Parameter estimates are loaded via `\input{}` from `includes/` subdirectory
- Figures are in `figures/` subdirectory (PDF format)
