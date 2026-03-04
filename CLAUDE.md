# CLAUDE.md — R3D Package

## What This Is

R package implementing regression discontinuity design with distribution-valued outcomes. Two estimators: local polynomial (pointwise quantile-by-quantile) and Frechet regression (global). Multiplier bootstrap for uniform inference. Supports sharp and fuzzy designs.

**Repo:** `Davidvandijcke/r3d` on GitHub
**Package:** `R3D` (R >= 4.0, compiled Fortran via `.Call`)
**Version:** 0.1.0
**Paper:** Van Dijcke (2025), "Regression Discontinuity Design with Distributional Outcomes"

## Project Structure

```
R3D/
├── R/                    # R source (6 files)
│   ├── r3d_main.R        # Main r3d() function — estimation entry point
│   ├── r3d_bootstrap.R   # Multiplier bootstrap for uniform inference
│   ├── r3d_bwselect.R    # Bandwidth selection (MSE-optimal)
│   ├── r3d_methods.R     # S3 methods: plot, print, summary
│   ├── r3d_utils.R       # Internal helpers (kernel weights, data prep)
│   └── R3D_package.R     # Package-level roxygen docs
├── src/                  # Compiled code
│   ├── locweights.f      # Fortran kernel weight computation (performance-critical)
│   ├── init.c            # .Call registration
│   ├── Makevars          # Unix build config (links LAPACK/BLAS)
│   └── Makevars.win      # Windows build config
├── tests/
│   └── testthat/
│       └── test-r3d-features.R   # Integration tests (sharp, fuzzy, bootstrap, bwselect)
├── man/                  # Roxygen-generated .Rd files
├── docs/                 # pkgdown site (auto-deployed to GitHub Pages)
├── stata_r3d/            # Stata port (separate CLAUDE.md inside)
├── DESCRIPTION           # Package metadata, version, dependencies
├── NAMESPACE             # Exports and imports (roxygen-managed)
├── NEWS.md               # Changelog
├── DEVELOPMENT.md        # Dev workflow guide
├── LICENSE / LICENSE.md   # MIT
└── .github/workflows/    # CI: R-CMD-check + pkgdown deploy
```

## Development Environment

```r
# Load in dev mode
devtools::load_all()

# Run tests
devtools::test()

# Full package check (run before every commit)
devtools::check()

# Regenerate docs after roxygen changes
devtools::document()

# Build pkgdown site locally
pkgdown::build_site()
```

**Dependencies:** `devtools`, `roxygen2`, `testthat`, `pkgdown` for development. Runtime: `stats`, `parallel`, `Rearrangement`, `Rdpack`, `Hmisc`.

## Key Architecture

- **Two estimators** in `r3d_main.R`: `method = "simple"` (pointwise local polynomial at each quantile) and `method = "frechet"` (global Frechet regression in Wasserstein space)
- **Fortran backend** (`src/locweights.f`): Computes kernel-weighted local polynomial matrices. Called from R via `.Fortran("locweights", ...)` with registered routines.
- **S3 class** `r3d`: returned by `r3d()`, with `plot.r3d`, `print.r3d`, `summary.r3d` methods
- **Bootstrap** (`r3d_bootstrap.R`): Multiplier bootstrap for uniform confidence bands. Parallelized via `parallel::mclapply`.
- **Bandwidth selection** (`r3d_bwselect.R`): MSE-optimal bandwidth, adapts rdrobust-style approach to distributional setting.

## Known Issues

- **`src/init.c` line 14**: Function named `R_init_YourPackageName` — should be `R_init_R3D` for proper DLL registration
- **`src/R3D.so`** tracked in git despite `src/*.so` in `.gitignore` (committed before rule existed)
- **README.md badge URL** has `yourusername` placeholder in one badge link

## Git Workflow

### Branch Naming

```
feat/<description>    # New features
fix/<description>     # Bug fixes
docs/<description>    # Documentation only
test/<description>    # Test additions/changes
```

### Conventional Commits

```
feat: add weighted regression support
fix: correct bandwidth selection for small samples
docs: update installation instructions
test: add fuzzy RD edge case tests
refactor: simplify kernel weight computation
```

### Before Every Commit

1. `devtools::document()` — regenerate roxygen docs if changed
2. `devtools::test()` — all tests pass
3. `devtools::check()` — no errors, no new warnings/notes

### Pre-PR Checklist

1. **NEWS.md** — Add entry describing the change
2. **DESCRIPTION** — Bump version if releasing
3. **`R CMD check`** — Zero errors on local machine
4. **Roxygen** — `devtools::document()` run, `man/` files current
5. **Tests** — New functionality has tests, existing tests pass
6. **CLAUDE.md** — Update test count, repo structure, Known Issues, architecture decisions if changed.
7. **CODEBASE_MAP.md** — Update if you added/removed/renamed files, classes, or modules.

## Coding Guidelines

**Simplicity first.** Minimum code that solves the problem. No speculative features.

**Surgical changes.** Touch only what you must. Match existing style. Don't refactor adjacent code.

**Goal-driven.** Define what "done" looks like before coding. Run tests to verify.

**R-specific:**
- Use `#' @export` roxygen tags, not manual NAMESPACE edits
- Internal functions: don't export, prefix with `.` if clarity helps
- Prefer base R over tidyverse in package code (fewer dependencies)
- Use `parallel::mclapply` for parallelization (already the pattern here)
- Fortran code changes require `R CMD INSTALL` to recompile

## What NOT to Do

- Don't commit `.Rhistory`, `.RData`, `.Rproj.user/`, compiled `.o`/`.so`/`.dll` files
- Don't manually edit `NAMESPACE` — it's roxygen-managed
- Don't add heavy dependencies without justification (CRAN submission goal)
- Don't edit `docs/` directly — it's pkgdown-generated

---

*Last updated: 2026-03-03*
