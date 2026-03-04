## R3D 0.1.1 (development)

### New Features

* Added `xi_mat` parameter to `r3d_bootstrap()` for deterministic bootstrap with pre-generated multiplier draws (enables cross-platform equivalence testing)
* Added Stata port (`stata_r3d/`) with full Mata implementation and 69-test equivalence suite

### Other Changes

* Added lint CI workflow
* Updated `.Rbuildignore` to exclude non-package files (`stata_r3d/`, `.claude/`, etc.)

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