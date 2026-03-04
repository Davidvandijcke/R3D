/*
 * example.do - Simple example for R3D package
 * This can be run after installing the package
 */

version 14.0
clear all
set more off

di as text _n "R3D Package Example"
di as text "{hline 40}"

// Generate example data
set obs 500
set seed 12345

// Running variable centered at 0
gen x = rnormal()

// Treatment (sharp RD)
gen treatment = (x >= 0)

// Generate distributional outcomes
// Each unit has 3 measurements from their outcome distribution
forvalues i = 1/3 {
    gen y`i' = 1 + 0.5*treatment + 0.2*x + rnormal(0, 1 + 0.3*treatment)
}

// Run R3D with default settings
di as text _n "Running R3D with sharp design..."
r3d x y1 y2 y3, cutoff(0) nquantiles(20)

// Show some results
di as text _n "Treatment effects at selected quantiles:"
di as text "  10th percentile: " as result %6.3f _b[q2]
di as text "  50th percentile: " as result %6.3f _b[q10] 
di as text "  90th percentile: " as result %6.3f _b[q18]

// Run with bandwidth selection
di as text _n "Running with automatic bandwidth selection..."
r3d x y*, cutoff(0) bwselect nquantiles(10) nograph

di as text _n "Example completed successfully!"
di as text "For more options, see: " as result "help r3d"