## R3D 0.1.1 (development)

### Documentation

* Added `\donttest{}` examples to `r3d_bootstrap()` and `r3d_bwselect()` (fixes CRAN ERROR for missing examples)
* Replaced `\dontrun{}` with `\donttest{}` in `r3d()`, `plot.r3d()`, `print.r3d()`, and `summary.r3d()` examples (fixes CRAN NOTE)
* Added `@noRd` to all internal helpers in `r3d_utils.R` to suppress NOTE about documented non-exported objects; removed orphaned `.Rd` files

### New Features

* Added `xi_mat` parameter to `r3d_bootstrap()` for deterministic bootstrap with pre-generated multiplier draws (enables cross-platform equivalence testing)
* Added Stata port (`stata_r3d/`) with full Mata implementation and 69-test equivalence suite

### Bug Fixes

* Fixed DLL registration: renamed `R_init_YourPackageName` to `R_init_R3D` so `R_registerRoutines` and `R_useDynamicSymbols` are called on load
* Fixed C/Fortran type mismatch: argument 7 of `locweights` changed from `double *KERNELW` to `int *KERNEL_TYPE` to match Fortran `INTEGER` declaration
* Fixed compound quoting bug in `r3d.ado` (lines 39-41) that caused option parsing failures when paths contained spaces

### Monte Carlo

* Added Stata Monte Carlo simulation framework (`02_mc_stata_worker.do`, `02_mc_stata_parallel.sh`, `run_stata_mc_nohup.sh`)
* Added `04_stata_coverage.R` for computing uniform coverage from Stata results
* Switched R MC scripts (`01_mc_R.R`, `03_mc_compare.R`) from pointwise to uniform coverage computation

### Other Changes

* Added lint CI workflow
* Updated `.Rbuildignore` to exclude non-package files (`stata_r3d/`, `.claude/`, etc.)
* Updated `.gitignore` to exclude MC run artifacts and Stata build logs
* Added `.Rbuildignore` patterns to exclude compiled artifacts (`src/*.so`, `src/*.o`, `src/*.dll`) from source tarball (CRAN requirement)

## R3D 0.1.0

### New Features

* Initial release of the R3D package
* Implementation of `r3d()` for distributional regression discontinuity
* Support for both sharp and fuzzy designs
* Two estimation methods: simple (pointwise) and Fréchet (global)
* Automatic bandwidth selection with `r3d_bwselect()`
* Multiplier bootstrap for uniform inference
* Visualization tools via `plot.r3d()`
* Comprehensive summary methods

### Bug Fixes

* None (initial release)

### Other Changes

* None (initial release)