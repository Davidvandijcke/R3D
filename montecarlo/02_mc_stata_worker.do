/*
 02_mc_stata_worker.do
 Parallelizable worker for R3D Monte Carlo simulations (Stata side).
 Generates its own DGP data internally — no dependency on R output.

 Usage: stata-mp -b do 02_mc_stata_worker.do <sim_start> <sim_end>
 Example: stata-mp -b do 02_mc_stata_worker.do 1 7

 Expects to be run from montecarlo/ directory.
 Requires: r3d installed in ~/ado/personal/ with lr3d.mlib
*/

clear all
set more off
set maxvar 32767

// Parse command-line arguments: sim_start and sim_end
local sim_start = `1'
local sim_end   = `2'

if missing(`sim_start') | missing(`sim_end') {
    di as error "Usage: stata-mp -b do 02_mc_stata_worker.do <sim_start> <sim_end>"
    exit 198
}

di "Worker starting: sims `sim_start' to `sim_end'"

// Use absolute path for output directory
local basedir "`c(pwd)'"
local outdir "`basedir'/output"

// ============================================================================
// CONFIGURATION
// ============================================================================
local BOOT_REPS  = 200
local ALPHA      = 0.05
local P_ORDER    = 2
local KERNEL     = "epanechnikov"
local NQ         = 9
local N_I        = 50       // draws per unit (same as R)

// Build quantile list for Stata syntax (space-separated) and Mata (comma-separated)
local q_list ""
local q_mata ""
forvalues i = 1/`NQ' {
    local q = `i' / 10
    local q_list "`q_list' `q'"
    if `i' == 1 {
        local q_mata "`q'"
    }
    else {
        local q_mata "`q_mata', `q'"
    }
}

capture mkdir "`outdir'/results_Stata"

// ============================================================================
// Mata: DGP data generation
// ============================================================================
mata:
mata clear

// Compute empirical quantiles of a vector at given probs
real rowvector emp_quantiles(real colvector y, real rowvector probs) {
    real scalar n, k, j
    real colvector sy
    real rowvector q
    real scalar idx, frac

    n = rows(y)
    k = cols(probs)
    sy = sort(y, 1)
    q = J(1, k, .)

    for (j = 1; j <= k; j++) {
        // Type 7 quantile (R default)
        idx = 1 + (n - 1) * probs[j]
        if (idx <= 1) {
            q[j] = sy[1]
        }
        else if (idx >= n) {
            q[j] = sy[n]
        }
        else {
            frac = idx - floor(idx)
            q[j] = (1 - frac) * sy[floor(idx)] + frac * sy[floor(idx) + 1]
        }
    }
    return(q)
}

// Generate DGP1 data: Normal with constant treatment effect
// Y_ij ~ N(5 + 5*X_i + delta*1(X>=0), |1+X_i|)  (min sd = 0.1)
// Returns (n x (1+NQ)) matrix: [X, Q1, ..., Q_NQ]
real matrix generate_dgp1(real scalar n, real scalar delta, real scalar n_i,
                          real rowvector q_grid) {
    real matrix result
    real colvector X, Y_ij
    real scalar i, mu_i, sigma_i, NQ

    NQ = cols(q_grid)
    result = J(n, 1 + NQ, .)

    X = runiform(n, 1) * 2 :- 1  // U(-1,1)
    result[., 1] = X

    for (i = 1; i <= n; i++) {
        mu_i = 5 + 5 * X[i] + delta * (X[i] >= 0)
        sigma_i = abs(1 + X[i])
        if (sigma_i < 0.1) sigma_i = 0.1
        Y_ij = rnormal(n_i, 1, mu_i, sigma_i)
        result[i, 2..(1+NQ)] = emp_quantiles(Y_ij, q_grid)
    }

    return(result)
}

// Generate DGP2 data: Normal + Exponential mixture
// Y_ij ~ N(2*X_i + delta*1(X>=0), 1) + 2*Exp(lambda_i)
// lambda_i = max(0.5, min(1.5, 1 + 0.5*X_i))
real matrix generate_dgp2(real scalar n, real scalar delta, real scalar n_i,
                          real rowvector q_grid) {
    real matrix result
    real colvector X, Y_ij, exp_draws
    real scalar i, mu_i, lambda_i, NQ

    NQ = cols(q_grid)
    result = J(n, 1 + NQ, .)

    X = runiform(n, 1) * 2 :- 1  // U(-1,1)
    result[., 1] = X

    for (i = 1; i <= n; i++) {
        mu_i = 2 * X[i] + delta * (X[i] >= 0)
        lambda_i = max((0.5, min((1.5, 1 + 0.5 * X[i]))))
        Y_ij = rnormal(n_i, 1, mu_i, 1)
        // Exp(lambda) = -ln(U)/lambda
        exp_draws = -ln(runiform(n_i, 1)) / lambda_i
        Y_ij = Y_ij + 2 * exp_draws
        result[i, 2..(1+NQ)] = emp_quantiles(Y_ij, q_grid)
    }

    return(result)
}

end

// ============================================================================
// Helper program: run r3d and save results
// ============================================================================
capture program drop mc_run_and_save
program define mc_run_and_save
    syntax, cell(string) sim(integer) nq(integer) outdir(string) ///
        method(string) porder(integer) kernel(string) qlist(string) ///
        bootreps(integer) ///
        tag(string) tests(string) [bwlist(string) fuzzy(string) denbw(string)]

    // Run r3d
    local yvars ""
    forvalues j = 1/`nq' {
        local yvars "`yvars' q`j'"
    }

    if "`bwlist'" != "" {
        capture noisily r3d x `yvars', cutoff(0) method(`method') ///
            polynomial(`porder') kernel(`kernel') ///
            quantiles(`qlist') ///
            bandwidth(`bwlist') ///
            bootstrap(`bootreps') ///
            tests(`tests') ///
            nograph
    }
    else {
        capture noisily r3d x `yvars', cutoff(0) method(`method') ///
            polynomial(`porder') kernel(`kernel') ///
            quantiles(`qlist') ///
            bwselect ///
            bootstrap(`bootreps') ///
            tests(`tests') ///
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
local methods "simple frechet"

log using "`outdir'/results_Stata/mc_stata_worker_`sim_start'_`sim_end'.log", replace text

di "=== R3D Monte Carlo Worker (Stata side) ==="
di "Sims: `sim_start' to `sim_end'  BOOT_REPS: `BOOT_REPS'"
di "Sample sizes: `sample_sizes'"
di "Output: `outdir'"
di ""

foreach dgp of local dgps {
    foreach n_obs of local sample_sizes {

        // DGP1: deltas 0, 1, 2; DGP2: delta 0 only
        if "`dgp'" == "dgp1" {
            local deltas "0 1 2"
        }
        else {
            local deltas "0"
        }

        foreach delta of local deltas {
            foreach method of local methods {

                local cell "`dgp'_sharp_`n_obs'_d`delta'_`method'"
                di _n "=== Cell: `cell' ==="

                forvalues s = `sim_start'/`sim_end' {

                    // Set seed to match R: data_seed = 1000 + s
                    local data_seed = 1000 + `s'

                    // Generate data in Mata
                    quietly mata: rseed(`data_seed')
                    quietly mata: _mc_data = generate_`dgp'(`n_obs', `delta', `N_I', (`q_mata'))

                    // Push to Stata dataset
                    quietly {
                        clear
                        set obs `n_obs'
                        gen double x = .
                        forvalues j = 1/`NQ' {
                            gen double q`j' = .
                        }
                        mata: st_store(., "x", _mc_data[., 1])
                        forvalues j = 1/`NQ' {
                            mata: st_store(., "q`j'", _mc_data[., 1 + `j'])
                        }
                        mata: mata drop _mc_data
                    }

                    // Run with Stata's own bandwidth selection
                    capture mc_run_and_save, ///
                        cell(`cell') sim(`s') nq(`NQ') outdir(`outdir') ///
                        method(`method') porder(`P_ORDER') kernel(`KERNEL') ///
                        qlist(`q_list') ///
                        bootreps(`BOOT_REPS') ///
                        tag(ownbw) tests(nullity homogeneity)

                    if _rc != 0 {
                        di "  Sim `s' FAILED (rc=" _rc ")"
                    }
                    else {
                        di "  Sim `s' OK"
                    }
                }
            }
        }
    }
}

di _n "=== Worker `sim_start'-`sim_end' complete ==="

log close
