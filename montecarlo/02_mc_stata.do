/*
 02_mc_stata.do
 Monte Carlo simulations for R3D: Stata side
 Loads R-exported data/xi/bw, runs Stata r3d, exports results.

 Usage: stata-se -b do 02_mc_stata.do
 Expects to be run from montecarlo/ directory.
 Requires: stata_r3d package installed (net install r3d, from("../stata_r3d/"))
*/

clear all
set more off
set maxvar 32767

// Install r3d from local source
capture net install r3d, from("../stata_r3d/") replace force

// Compile Mata if needed
quietly capture mata: mata which r3d_compute_quantiles
if _rc != 0 {
    quietly capture do "../stata_r3d/mata/r3d_mata.mata"
}
quietly capture mata: mata which r3d_plugin_available
if _rc != 0 {
    quietly capture do "../stata_r3d/mata/r3d_plugin_interface.mata"
}

// ============================================================================
// CONFIGURATION
// ============================================================================
local N_STATA    = 100
local BOOT_REPS  = 200
local ALPHA      = 0.05
local P_ORDER    = 2
local KERNEL     = "epanechnikov"
local NQ         = 9

// Build quantile list: (1:9)/10
local q_list ""
forvalues i = 1/`NQ' {
    local q = `i' / 10
    local q_list "`q_list' `q'"
}

local outdir "output"
capture mkdir "`outdir'/results_Stata"

// ============================================================================
// Helper program: run r3d and save results
// ============================================================================
capture program drop mc_run_and_save
program define mc_run_and_save
    syntax, cell(string) sim(integer) nq(integer) outdir(string) ///
        method(string) porder(integer) kernel(string) qlist(string) ///
        bwlist(string) bootreps(integer) xifile(string) ///
        tag(string) tests(string) [fuzzy(string) denbw(string)]

    // Run r3d
    local yvars ""
    forvalues j = 1/`nq' {
        local yvars "`yvars' q`j'"
    }

    if "`fuzzy'" != "" & "`denbw'" != "" {
        capture noisily r3d x `yvars', cutoff(0) method(`method') ///
            fuzzy(`fuzzy') ///
            polynomial(`porder') kernel(`kernel') ///
            quantiles(`qlist') ///
            bandwidth(`bwlist') denbandwidth(`denbw') ///
            bootstrap(`bootreps') ///
            tests(`tests') ///
            ximat("`xifile'") ///
            nograph
    }
    else if "`bwlist'" != "" {
        capture noisily r3d x `yvars', cutoff(0) method(`method') ///
            polynomial(`porder') kernel(`kernel') ///
            quantiles(`qlist') ///
            bandwidth(`bwlist') ///
            bootstrap(`bootreps') ///
            tests(`tests') ///
            ximat("`xifile'") ///
            nograph
    }
    else {
        capture noisily r3d x `yvars', cutoff(0) method(`method') ///
            polynomial(`porder') kernel(`kernel') ///
            quantiles(`qlist') ///
            bwselect ///
            bootstrap(`bootreps') ///
            tests(`tests') ///
            ximat("`xifile'") ///
            nograph
    }

    if _rc != 0 {
        di "  Sim `sim' `tag' FAILED (rc=" _rc ")"
        exit _rc
    }

    // Extract results
    tempname tau_v cb_lo cb_hi bw_v pv
    matrix `tau_v' = e(b)
    matrix `bw_v' = e(bandwidth_num)
    matrix `pv' = e(pvalues)
    capture matrix `cb_lo' = e(cb_lower)
    local has_cb = (_rc == 0)
    if `has_cb' {
        capture matrix `cb_hi' = e(cb_upper)
    }

    // Write tau/cb/bw CSV
    preserve
    quietly {
        clear
        set obs `nq'
        gen q = .
        gen tau = .
        gen cb_lower = .
        gen cb_upper = .
        gen bw = .
        forvalues j = 1/`nq' {
            replace q = `= `j'/10' in `j'
            replace tau = `tau_v'[1,`j'] in `j'
            if `has_cb' {
                replace cb_lower = `cb_lo'[1,`j'] in `j'
                replace cb_upper = `cb_hi'[1,`j'] in `j'
            }
            replace bw = `bw_v'[1,`j'] in `j'
        }
        local outfile "`outdir'/results_Stata/mc_res_`cell'_`tag'_`sim'.csv"
        export delimited using "`outfile'", replace
    }
    restore

    // Write p-value CSV
    preserve
    quietly {
        clear
        set obs 1
        gen pval_nullity = `pv'[1,1]
        capture gen pval_homogeneity = `pv'[1,2]
        local pvfile "`outdir'/results_Stata/mc_pval_`cell'_`tag'_`sim'.csv"
        export delimited using "`pvfile'", replace
    }
    restore
end

// ============================================================================
// MAIN LOOP
// ============================================================================
local dgps "dgp1 dgp2"
local sample_sizes "200 500"
local deltas "0 2"
local methods "simple frechet"

log using "`outdir'/results_Stata/mc_stata.log", replace text

di "=== R3D Monte Carlo Simulation (Stata side) ==="
di "N_STATA: `N_STATA'  BOOT_REPS: `BOOT_REPS'"
di "Sample sizes: `sample_sizes'"
di "Deltas: `deltas'"
di "Methods: `methods'"
di ""

foreach dgp of local dgps {
    foreach n_obs of local sample_sizes {
        foreach delta of local deltas {
            foreach method of local methods {

                local cell "`dgp'_sharp_`n_obs'_d`delta'_`method'"
                di _n "=== Cell: `cell' ==="

                forvalues s = 1/`N_STATA' {

                    local datafile "`outdir'/data/mc_`dgp'_`n_obs'_`delta'_`s'.csv"
                    local xifile   "`outdir'/xi/mc_xi_`n_obs'_`s'.csv"
                    local bwfile   "`outdir'/results_R/mc_bw_`dgp'_`method'_`n_obs'_`delta'_`s'.csv"

                    // Check if data file exists
                    capture confirm file "`datafile'"
                    if _rc != 0 {
                        di "  Sim `s': data file missing, skipping"
                        continue
                    }

                    // Load data
                    quietly import delimited using "`datafile'", clear varnames(1)

                    // --------------------------------------------------------
                    // RUN 1: Shared bandwidths (R bandwidths)
                    // --------------------------------------------------------
                    capture confirm file "`bwfile'"
                    if _rc == 0 {
                        // Read R bandwidths
                        quietly {
                            preserve
                            import delimited using "`bwfile'", clear varnames(1)
                            local bw_list ""
                            local bw_n = _N
                            forvalues j = 1/`bw_n' {
                                local bw_val = bw[`j']
                                local bw_list "`bw_list' `bw_val'"
                            }
                            restore
                        }

                        capture mc_run_and_save, ///
                            cell(`cell') sim(`s') nq(`NQ') outdir(`outdir') ///
                            method(`method') porder(`P_ORDER') kernel(`KERNEL') ///
                            qlist(`q_list') bwlist(`bw_list') ///
                            bootreps(`BOOT_REPS') xifile(`xifile') ///
                            tag(sharedbw) tests(nullity homogeneity)
                    }

                    // --------------------------------------------------------
                    // RUN 2: Own bandwidths (Stata bwselect)
                    // --------------------------------------------------------
                    // Re-load data
                    quietly import delimited using "`datafile'", clear varnames(1)

                    capture mc_run_and_save, ///
                        cell(`cell') sim(`s') nq(`NQ') outdir(`outdir') ///
                        method(`method') porder(`P_ORDER') kernel(`KERNEL') ///
                        qlist(`q_list') bwlist() ///
                        bootreps(`BOOT_REPS') xifile(`xifile') ///
                        tag(ownbw) tests(nullity homogeneity)

                    if mod(`s', 10) == 0 {
                        di "  Completed `s'/`N_STATA' for `cell'"
                    }
                }
            }
        }
    }
}

// ============================================================================
// FUZZY CELL: DGP1, n=200, delta=2, simple
// ============================================================================
di _n "=== Fuzzy cell: dgp1_fuzzy_200_d2_simple ==="

forvalues s = 1/`N_STATA' {

    local datafile "`outdir'/data/mc_dgp1fuzzy_200_2_`s'.csv"
    local xifile   "`outdir'/xi/mc_xi_200_`s'.csv"
    local bwfile   "`outdir'/results_R/mc_bw_dgp1fuzzy_simple_200_2_`s'.csv"
    local bwdenfile "`outdir'/results_R/mc_bwden_dgp1fuzzy_simple_200_2_`s'.csv"

    capture confirm file "`datafile'"
    if _rc != 0 {
        di "  Sim `s': data file missing, skipping"
        continue
    }

    // Load data
    quietly import delimited using "`datafile'", clear varnames(1)

    // Read R bandwidths
    capture confirm file "`bwfile'"
    if _rc == 0 {
        quietly {
            preserve
            import delimited using "`bwfile'", clear varnames(1)
            local bw_list ""
            local bw_n = _N
            forvalues j = 1/`bw_n' {
                local bw_val = bw[`j']
                local bw_list "`bw_list' `bw_val'"
            }
            restore
        }

        // Read denominator bandwidth
        local h_den = 0
        capture confirm file "`bwdenfile'"
        if _rc == 0 {
            quietly {
                preserve
                import delimited using "`bwdenfile'", clear varnames(1)
                local h_den = bw_den[1]
                restore
            }
        }

        // Re-load data
        quietly import delimited using "`datafile'", clear varnames(1)

        capture mc_run_and_save, ///
            cell(dgp1_fuzzy_200_d2_simple) sim(`s') nq(`NQ') outdir(`outdir') ///
            method(simple) porder(`P_ORDER') kernel(`KERNEL') ///
            qlist(`q_list') bwlist(`bw_list') ///
            bootreps(`BOOT_REPS') xifile(`xifile') ///
            tag(sharedbw) tests(nullity) ///
            fuzzy(t_treat) denbw(`h_den')
    }

    if mod(`s', 10) == 0 {
        di "  Completed `s'/`N_STATA' for fuzzy"
    }
}

di _n "=== Stata Monte Carlo complete ==="

log close
