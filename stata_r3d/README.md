# R3D for Stata

This is the Stata implementation of the R3D (Regression Discontinuity with Distributional Outcomes) package.

## Installation

### Recommended: Using net install

The easiest way to install R3D is using Stata's `net install` command:

```stata
* Install directly from the GitHub directory
net install r3d, from("https://raw.githubusercontent.com/username/repo/main/stata_r3d/") replace

* OR install from a local directory
net install r3d, from("/path/to/stata_r3d/") replace

* Run post-installation setup
r3d_setup
```

### Alternative: Manual Installation

1. Copy the `.ado` files from `ado/` to your Stata ado-path:
   ```stata
   . adopath
   ```
   
2. Copy the `.mata` file from `mata/` to the same directory.

3. Copy the help files from `doc/` to the same directory.

4. The package will work using pure Mata code.

### Optional: Building the Fortran Plugin for Performance

For better performance with large datasets, you can build the Fortran plugin:

1. Install prerequisites:
   - GNU Fortran compiler (gfortran)
   - LAPACK libraries
   - Stata SDK (for plugin development)

2. Build the plugin:
   ```bash
   cd plugin/
   # Edit Makefile to set correct STATA_INCLUDE path
   make
   ```

3. Copy `r3d_plugin.plugin` to your ado directory

## Basic Usage

```stata
* Sharp RD with distributional outcome
r3d running_var y1 y2 y3 y4 y5, cutoff(0)

* Fuzzy RD with bandwidth selection
r3d score outcome*, cutoff(70) fuzzy(treatment) bwselect

* Frechet method with bootstrap inference
r3d x y*, c(0) method(frechet) bootstrap(999)

* Specific quantiles only
r3d x y1-y10, c(0) quantiles(0.1 0.5 0.9)
```

## Commands

### Main Command
- `r3d` - Main estimation command

### Auxiliary Commands
- `r3d_bwselect` - Standalone bandwidth selection
- `r3d_bootstrap` - Post-estimation bootstrap inference

## Features

- **Sharp and Fuzzy RD**: Supports both designs
- **Two Methods**: Simple (pointwise) and Frechet (with isotonic constraints)
- **Flexible Kernels**: Triangular, Epanechnikov, and Uniform
- **Data-driven Bandwidth Selection**: MSE-optimal and IMSE-optimal
- **Bootstrap Inference**: Uniform confidence bands and hypothesis tests
- **Multiple Outcomes**: Handles repeated measurements naturally

## File Structure

```
stata_r3d/
├── ado/
│   ├── r3d.ado          # Main command
│   ├── r3d_bootstrap.ado # Bootstrap inference
│   └── r3d_bwselect.ado  # Bandwidth selection
├── mata/
│   └── r3d_mata.mata    # Mata functions
├── plugin/
│   ├── r3d_fortran.f90  # Fortran code
│   ├── r3d_plugin.c     # C wrapper
│   ├── stplugin.h       # Plugin header (stub)
│   └── Makefile         # Build configuration
├── doc/
│   └── r3d.sthlp        # Help file
└── test/
    └── test_r3d.do      # Test script
```

## Requirements

- Stata 14.0 or later
- For plugin: GNU Fortran compiler and LAPACK

## Notes

1. The plugin provides significant speed improvements for large datasets but is optional.

2. The `stplugin.h` file is a stub. Replace it with the actual header from the Stata SDK.

3. Multiple outcome variables (y1, y2, ...) are treated as repeated measurements from each unit's distribution.

## Citation

If you use this package, please cite:

Van Dijcke, D. (2025). "Regression Discontinuity Designs with Distributional Outcomes." Working Paper.

## License

This software is provided under the same license as the original R package.