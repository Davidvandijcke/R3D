*! r3d 1.0.0 - Regression Discontinuity with Distribution-valued Outcomes
*! Author: David Van Dijcke
*! Date: 2025

program define r3d, eclass
    version 18.0
    
    // Load Mata functions: try mlib first, then compile from source
    quietly capture mata: mata which r3d_compute_quantiles
    if _rc != 0 {
        mata: mata mlib index
        quietly capture mata: mata which r3d_compute_quantiles
    }
    if _rc != 0 {
        quietly capture do "`c(sysdir_plus)'r/r3d_mata.mata"
        if _rc != 0 {
            quietly capture do "mata/r3d_mata.mata"
        }
    }
    quietly capture mata: mata which r3d_plugin_available
    if _rc != 0 {
        quietly capture do "`c(sysdir_plus)'r/r3d_plugin_interface.mata"
        if _rc != 0 {
            quietly capture do "mata/r3d_plugin_interface.mata"
        }
    }

    // Load plugin if available (defaulting to the Fortran backend)
    capture program drop r3d_plugin_load
    quietly capture do "`c(sysdir_plus)'r/r3d_plugin_load.ado"
    if _rc == 0 {
        quietly capture r3d_plugin_load
    }
    else {
        global R3D_USE_PLUGIN 0
    }
    
    // Normalize aliases to DiSCo-style option names before parsing
    local cmdline = regexr("`0'", "(?i)nquantiles\(", "nq(")
    local cmdline = regexr("`cmdline'", "(?i)denbandwidth\(", "denband(")
    local 0 = trim("`cmdline'")

    syntax varlist(min=2 numeric) [if] [in], ///
        [ cutoff(real 0) ///
          method(string) ///
          fuzzy(varname) ///
          polynomial(integer 2) ///
          pilot(integer 1) ///
          bandwidth(numlist) ///
          denband(numlist) ///
          bwselect ///
          kernel(string) ///
          quantiles(numlist) ///
          nq(integer 99) ///
          coverage ///
          bootstrap(integer 0) ///
          tests(string) ///
          testranges(string) ///
          level(cilevel) ///
          weights(varname) ///
          ximat(string) ///
          saving(string) ///
          replace ///
          nograph ]

    local denbandwidth `denband'
    local nquantiles = real("`nq'")
    local nograph_flag = cond("`nograph'" == "", 0, 1)

    // Set defaults for options
    local cutoff = real("`cutoff'")
    local method = strlower(strtrim("`method'"))
    if "`method'" == "" local method "simple"
    local polynomial = real("`polynomial'")
    local pilot = real("`pilot'")
    local kernel = strlower(strtrim("`kernel'"))
    if "`kernel'" == "" local kernel "epanechnikov"
    local bootstrap = real("`bootstrap'")
    local level_val = cond("`level'" == "", c(level), real("`level'"))
    local coverage_flag = cond("`coverage'" == "", 0, 1)
    local tests_clean = strtrim(strlower(subinstr("`tests'", ",", " ", .)))
    if "`tests_clean'" == "" local tests_clean "none"
    local tests = "`tests_clean'"
    local testranges = strtrim("`testranges'")
    local nograph = `nograph_flag'

    // Parse variables
    tokenize `varlist'
    local xvar `1'
    local yvars : list varlist - xvar
    
    // Mark sample
    marksample touse
    quietly count if `touse'
    if r(N) == 0 error 2000
    local N = r(N)
    
    if !inlist("`method'", "simple", "frechet") {
        di as error "method() must be either 'simple' or 'frechet'"
        exit 198
    }

    if !inlist("`kernel'", "triangular", "epanechnikov", "uniform") {
        di as error "kernel() must be 'triangular', 'epanechnikov', or 'uniform'"
        exit 198
    }

    local method_code = cond("`method'" == "frechet", 1, 0)

    if `polynomial' < 0 {
        di as error "polynomial() must be non-negative"
        exit 198
    }

    if `pilot' < 0 {
        di as error "pilot() must be non-negative"
        exit 198
    }

    if `bootstrap' < 0 {
        di as error "bootstrap() must be non-negative"
        exit 198
    }

    // Set kernel type for plugin
    if "`kernel'" == "triangular" local kernel_type 1
    else if "`kernel'" == "epanechnikov" local kernel_type 2
    else local kernel_type 3
    
    // Generate quantile grid
    if "`quantiles'" == "" {
        if `nquantiles' <= 0 {
            di as error "nquantiles() must be positive"
            exit 198
        }
        local quantiles
        forvalues i = 1/`nquantiles' {
            local q = `i' / (`nquantiles' + 1)
            local quantiles `quantiles' `q'
        }
    }
    else {
        numlist "`quantiles'", range(>0 <1)
        local quantiles `r(numlist)'
    }
    local nq : word count `quantiles'
    
    // Center running variable at cutoff
    tempvar x_centered
    quietly gen double `x_centered' = `xvar' - `cutoff' if `touse'
    
    // Check if fuzzy RD
    local is_fuzzy 0
    if "`fuzzy'" != "" {
        confirm numeric variable ``fuzzy''
        local is_fuzzy 1
        tempvar t_centered
        quietly gen double `t_centered' = `fuzzy' if `touse'
    }
    
    if "`weights'" != "" {
        confirm numeric variable ``weights''
    }

    // Compute empirical quantiles for each observation
    di as text "Computing empirical quantiles..."
    tempname Qmat
    mata: r3d_compute_quantiles("`yvars'", "`quantiles'", "`touse'", "`Qmat'", "`weights'")
    
    // Determine bandwidths (automatic if requested or none supplied)
    local use_bwselect = ("`bwselect'" != "" | "`bandwidth'" == "")

    tempname HNUM
    tempname HNUM_RC
    scalar `HNUM_RC' = 0

    if `use_bwselect' {
        di as text "Selecting bandwidth..."
        local tvarname = cond(`is_fuzzy', "`t_centered'", "")
        mata: r3d_bandwidth_select("`x_centered'", "`Qmat'", "`quantiles'", "`tvarname'", ///
            `polynomial', `pilot', `kernel_type', "`method'", "`touse'", `is_fuzzy', `coverage_flag', "`weights'")

        matrix `HNUM' = r(h_num)
        mata: st_numscalar("`HNUM_RC'", r3d_prepare_bandwidth_matrix("`HNUM'", `nq', "`method'"))
        if scalar(`HNUM_RC') != 0 {
            di as error "Bandwidth selection failed; inspect data or specify bandwidth() manually"
            exit 498
        }

        if `is_fuzzy' {
            local h_den_scalar = r(h_den)
        }
        else {
            local h_den_scalar = 0
        }

        // Report selection summary
        mata: st_numscalar("__r3d_hmin", min(st_matrix("`HNUM'")))
        mata: st_numscalar("__r3d_hmax", max(st_matrix("`HNUM'")))
        di as text "Selected bandwidth range: " as result %9.4f scalar(__r3d_hmin) ///
            as text " to " as result %9.4f scalar(__r3d_hmax)
        scalar drop __r3d_hmin __r3d_hmax
    }
    else {
        // User-supplied numerator bandwidths
        mata: st_matrix("`HNUM'", (strtoreal(tokens("`bandwidth'")))')
        mata: st_numscalar("`HNUM_RC'", r3d_prepare_bandwidth_matrix("`HNUM'", `nq', "`method'"))
        if scalar(`HNUM_RC') != 0 {
            di as error "bandwidth() must be a scalar or vector of length matching quantiles"
            exit 198
        }

        if `is_fuzzy' {
            if "`denbandwidth'" == "" {
                di as error "For fuzzy designs, supply denbandwidth() when bandwidth() is provided"
                exit 198
            }
            local bwden_ct : word count `denbandwidth'
            if `bwden_ct' != 1 {
                di as error "denbandwidth() must contain a single value"
                exit 198
            }
            local h_den_scalar : word 1 of `denbandwidth'
        }
        else {
            local h_den_scalar = 0
        }
    }
    
    // Main estimation
    di as text "Estimating treatment effects..."

    tempname tau se alpha_plus alpha_minus w_plus w_minus e1 e2 int_plus int_minus alpha_t_plus alpha_t_minus w_t_plus w_t_minus denom cb_lower cb_upper bw_num

    local tvarname = cond(`is_fuzzy', "`t_centered'", "")

    if "`method'" == "simple" {
        mata: r3d_simple("`x_centered'", "`Qmat'", "`tvarname'", "`HNUM'", `h_den_scalar', ///
            `polynomial', `kernel_type', "`tau'", "`se'", "`alpha_plus'", "`alpha_minus'", ///
            "`w_plus'", "`w_minus'", "`e1'", "`e2'", "`int_plus'", "`int_minus'", ///
            "`alpha_t_plus'", "`alpha_t_minus'", "`w_t_plus'", "`w_t_minus'", "`denom'", ///
            "`touse'", `is_fuzzy')
    }
    else {
        mata: r3d_frechet("`x_centered'", "`Qmat'", "`tvarname'", "`HNUM'", `h_den_scalar', ///
            `polynomial', `kernel_type', "`tau'", "`se'", "`alpha_plus'", "`alpha_minus'", ///
            "`w_plus'", "`w_minus'", "`e1'", "`e2'", "`int_plus'", "`int_minus'", ///
            "`alpha_t_plus'", "`alpha_t_minus'", "`w_t_plus'", "`w_t_minus'", "`denom'", ///
            "`touse'", `is_fuzzy')
    }
    
    // Bootstrap inference if requested
    tempname bs_results pvals
    matrix `pvals' = J(1, 3, .)

    if `bootstrap' > 0 {
        di as text "Running bootstrap with " as result `bootstrap' as text " replications..."

        // Load xi_mat from CSV if provided
        local xi_mat_name ""
        if "`ximat'" != "" {
            tempname XI_MAT
            mata: st_matrix("`XI_MAT'", _r3d_load_csv("`ximat'"))
            local xi_mat_name "`XI_MAT'"
        }

        mata: r3d_bootstrap("`x_centered'", "`Qmat'", "`tvarname'", "`HNUM'", `h_den_scalar', ///
            `polynomial', `kernel_type', `bootstrap', `level_val', ///
            "`alpha_plus'", "`alpha_minus'", "`w_plus'", "`w_minus'", ///
            "`e1'", "`e2'", "`int_plus'", "`int_minus'", "`tau'", "`denom'", ///
            `method_code', "`quantiles'", "`tests'", "`testranges'", ///
            "`bs_results'", "`pvals'", "`touse'", `is_fuzzy', "`xi_mat_name'")

        matrix `cb_lower' = r(cb_lower)
        matrix `cb_upper' = r(cb_upper)
        // pvals already populated via st_matrix() in Mata
    }
    
    // Store results
    tempname b V
    matrix `b' = `tau'
    
    if `bootstrap' > 0 {
        mata: st_matrix("`V'", diag((st_matrix("`se'"):^2)'))
    }
    else {
        matrix `V' = J(`nq', `nq', 0)
        forvalues i = 1/`nq' {
            matrix `V'[`i',`i'] = `se'[1,`i']^2
        }
    }
    
    // Column names for matrices
    local colnames
    local i = 1
    foreach q of local quantiles {
        local colnames `colnames' q`i'
        local i = `i' + 1
    }
    matrix colnames `b' = `colnames'
    matrix colnames `V' = `colnames'
    matrix rownames `V' = `colnames'

    matrix `bw_num' = `HNUM'

    // Post results
    ereturn post `b' `V', esample(`touse') obs(`N')

    // Additional stored results
    ereturn local cmd "r3d"
    ereturn local method "`method'"
    ereturn local kernel "`kernel'"
    ereturn scalar cutoff = `cutoff'
    ereturn scalar polynomial = `polynomial'
    ereturn matrix bandwidth_num = `bw_num'
    ereturn scalar coverage = `coverage_flag'
    ereturn scalar pilot = `pilot'
    ereturn scalar method_code = `method_code'
    if `is_fuzzy' ereturn scalar bandwidth_den = `h_den_scalar'
    ereturn scalar nquantiles = `nq'
    tempname qmat
    matrix `qmat' = J(1, `nq', .)
    forvalues i = 1/`nq' {
        matrix `qmat'[1,`i'] = `: word `i' of `quantiles''
    }
    ereturn matrix quantiles = `qmat'
    
    if `bootstrap' > 0 {
        ereturn scalar bootstrap = `bootstrap'
    }

    ereturn matrix pvalues = `pvals'
    if `bootstrap' > 0 {
        ereturn matrix tau_boot = `bs_results'
        ereturn matrix cb_lower = `cb_lower'
        ereturn matrix cb_upper = `cb_upper'
    }

    ereturn matrix alpha_plus = `alpha_plus'
    ereturn matrix alpha_minus = `alpha_minus'
    ereturn matrix w_plus = `w_plus'
    ereturn matrix w_minus = `w_minus'
    ereturn matrix e1 = `e1'
    ereturn matrix int_plus = `int_plus'
    ereturn matrix int_minus = `int_minus'
    if `is_fuzzy' {
        ereturn matrix e2 = `e2'
        ereturn matrix alpha_t_plus = `alpha_t_plus'
        ereturn matrix alpha_t_minus = `alpha_t_minus'
        ereturn matrix w_t_plus = `w_t_plus'
        ereturn matrix w_t_minus = `w_t_minus'
        ereturn scalar denom = scalar(`denom')
    }
    else {
        ereturn scalar denom = 1
    }
    ereturn local xvar "`xvar'"
    if `is_fuzzy' ereturn local tvar "`fuzzy'"

    local h_num_display = `= `HNUM'[1,1]'
    local h_den_display = `h_den_scalar'

    // Display results
    di _n as text "Regression Discontinuity with Distributional Outcomes"
    di as text "{hline 78}"
    di as text "Method: " as result "`method'" _col(40) as text "Kernel: " as result "`kernel'"
    di as text "Polynomial order: " as result `polynomial' _col(40) ///
        as text "Observations: " as result `N'
    di as text "Bandwidth (Y): " as result %9.4f `h_num_display' _col(40) ///
        as text "Cutoff: " as result %9.4f `cutoff'
    if `is_fuzzy' di as text "Bandwidth (T): " as result %9.4f `h_den_display' ///
        _col(40) as text "Design: " as result "Fuzzy"
    else di _col(40) as text "Design: " as result "Sharp"
    
    // Plot results unless suppressed
    if !`nograph' {
        r3d_plot, quantiles(`quantiles') level(`level_val')
    }
    
    // Save results if requested
    if "`saving'" != "" {
        preserve
        clear
        quietly {
            set obs `nq'
            gen quantile = .
            gen tau = .
            gen se = .
            if `bootstrap' > 0 {
                gen ci_lower = .
                gen ci_upper = .
            }
            
            forvalues i = 1/`nq' {
                replace quantile = `: word `i' of `quantiles'' in `i'
                replace tau = `tau'[1,`i'] in `i'
                replace se = `se'[1,`i'] in `i'
                if `bootstrap' > 0 {
                    replace ci_lower = tau - invnormal(1-(1-`level_val'/100)/2)*se in `i'
                    replace ci_upper = tau + invnormal(1-(1-`level_val'/100)/2)*se in `i'
                }
            }
        }
        save "`saving'", `replace'
        restore
    }
end

// Plot results
program define r3d_plot
    syntax, quantiles(numlist) [level(cilevel)]
    
    preserve
    
    tempvar quantile tau ci_lower ci_upper
    
    quietly {
        clear
        local nq : word count `quantiles'
        set obs `nq'
        
        gen `quantile' = .
        gen `tau' = .
        gen `ci_lower' = .
        gen `ci_upper' = .
        
        forvalues i = 1/`nq' {
            replace `quantile' = `: word `i' of `quantiles'' in `i'
            replace `tau' = _b[q`i'] in `i'
            replace `ci_lower' = _b[q`i'] - invnormal(1-(1-`level'/100)/2)*_se[q`i'] in `i'
            replace `ci_upper' = _b[q`i'] + invnormal(1-(1-`level'/100)/2)*_se[q`i'] in `i'
        }
    }
    
    twoway (rarea `ci_lower' `ci_upper' `quantile', color(gs14) lwidth(none)) ///
           (line `tau' `quantile', lcolor(navy) lwidth(medthick)) ///
           (line `quantile' `quantile' if `tau'==0, lcolor(gs8) lpattern(dash)), ///
           ytitle("Treatment Effect") xtitle("Quantile") ///
           legend(order(2 "Point Estimate" 1 "`level'% CI") rows(1)) ///
           title("RD Treatment Effects by Quantile") ///
           scheme(s2color)
    
    restore
end
