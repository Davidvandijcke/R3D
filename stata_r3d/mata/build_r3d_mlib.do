/*
 * build_r3d_mlib.do - Compile Mata source files into lr3d.mlib
 *
 * Usage:
 *   cd mata/
 *   do build_r3d_mlib.do
 *
 * Output: lr3d.mlib in the current directory
 */

clear all

capture log using build_r3d_mlib, replace

mata:
mata clear
mata set matastrict off
end

// Load all Mata source files
do r3d_mata.mata
do r3d_plugin_interface.mata

// Now create and populate the Mata library
mata:
mata mlib create lr3d, replace

// Functions from r3d_mata.mata
mata mlib add lr3d r3d_kernel()
mata mlib add lr3d r3d_variance()
mata mlib add lr3d r3d_compute_quantiles()
mata mlib add lr3d r3d_estimate_density()
mata mlib add lr3d r3d_kernel_matrices()
mata mlib add lr3d r3d_fit_global_poly()
mata mlib add lr3d r3d_bandwidth_select()
mata mlib add lr3d r3d_prepare_bandwidth_matrix()
mata mlib add lr3d r3d_locpoly()
mata mlib add lr3d r3d_isotonic()
mata mlib add lr3d r3d_estimate_core()
mata mlib add lr3d r3d_simple()
mata mlib add lr3d r3d_frechet()
mata mlib add lr3d r3d_bootstrap()
mata mlib add lr3d r3d_gini_from_quantile()
mata mlib add lr3d r3d_test_gini()
mata mlib add lr3d _r3d_load_csv()
mata mlib add lr3d _r3d_test_locpoly()

// Functions from r3d_plugin_interface.mata
mata mlib add lr3d r3d_plugin_available()
mata mlib add lr3d r3d_plugin_locweights()

end

di as text "lr3d.mlib created successfully."

capture log close
