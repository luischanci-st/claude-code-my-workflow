---
paths:
  - "**/*.R"
  - "Chkb/**/*.R"
  - "scripts/**/*.R"
---

# R Code Standards

**Standard:** Senior Principal Data Engineer + PhD researcher quality

---

## 1. Reproducibility

- All packages loaded at top via `library()` (not `require()`)
- All paths relative to repository root
- `dir.create(..., recursive = TRUE)` for output directories

## 2. Function Design

- `snake_case` naming, verb-noun pattern
- Roxygen-style documentation
- Default parameters, no magic numbers
- Named return values (lists or tibbles)

## 3. Domain Correctness

- Spatial weight matrix construction: directed graph from C18 interbank data, row-normalized
- CPI deflation: match deflator period to observation date; real CLP units
- Bank ID filtering: 20-bank sample (`c("OR","PP","PS",...)`)
- igraph: always use `directed = TRUE` for interbank networks

## 4. Visual Identity

```r
# --- Academic publication palette (neutral) ---
primary_dark  <- "#2c3e50"
accent_blue   <- "#2980b9"
accent_gray   <- "#7f8c8d"
positive_green <- "#27ae60"
negative_red  <- "#c0392b"
```

### Custom Theme
```r
theme_paper <- function(base_size = 11) {
  theme_minimal(base_size = base_size) +
    theme(
      plot.title = element_text(face = "bold"),
      legend.position = "bottom",
      plot.background = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA)
    )
}
```

### Figure Dimensions for Journal
```r
ggsave(filepath, width = 6.5, height = 4.5, bg = "white", dpi = 300)
```

## 5. RDS Data Pattern

**Heavy computations saved as RDS; paper rendering loads pre-computed data.**

```r
saveRDS(result, file.path(out_dir, "descriptive_name.rds"))
```

## 6. Common Pitfalls

| Pitfall | Impact | Prevention |
|---------|--------|------------|
| Missing `bg = "white"` | Transparent bg in paper | Always include in ggsave() |
| Hardcoded paths | Breaks on other machines | Use relative paths or `here::here()` |
| RCall interop issues | Julia-R data type mismatches | Test with small data first |
| igraph undirected default | Wrong network structure | Explicitly set `directed = TRUE` |
| CPI base year mismatch | Incorrect deflation | Verify CPI series matches data period |

## 7. Line Length & Mathematical Exceptions

**Standard:** Keep lines <= 100 characters.

**Exception: Mathematical Formulas** -- lines may exceed 100 chars **if and only if:**

1. Breaking the line would harm readability of the math (matrix ops, formula implementations matching paper equations)
2. An inline comment explains the mathematical operation
3. The line is in a numerically intensive section (estimation routines, data transformations)

**Quality Gate Impact:**
- Long lines in non-mathematical code: minor penalty (-1 to -2 per line)
- Long lines in documented mathematical sections: no penalty

## 8. Code Quality Checklist

```
[ ] Packages at top via library()
[ ] All paths relative
[ ] Functions documented (Roxygen)
[ ] Figures: white bg, 300 DPI, explicit dimensions
[ ] RDS: every computed object saved
[ ] Comments explain WHY not WHAT
```
