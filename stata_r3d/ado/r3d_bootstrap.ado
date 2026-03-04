*! r3d_bootstrap 2.0.0 - Bootstrap inference for R3D
*! Author: David Van Dijcke (Stata port updated)

program define r3d_bootstrap, rclass
    version 18.0

    syntax [anything], ///
        [Reps(integer 999) ///
        Level(cilevel) ///
        Tests(string) ///
        TESTRanges(string) ///
        SEED(integer -1) ///
        XImat(string) ///
        NOIsily]

    // Ensure r3d has been run
    if "`e(cmd)'" != "r3d" {
        di as error "r3d_bootstrap requires prior r3d results"
        exit 301
    }

    // Set RNG seed if requested
    if `seed' != -1 {
        set seed `seed'
    }

    // Defaults and option normalization
    local level_val = cond("`level'" == "", c(level), `level')
    local tests_clean = lower(trim(subinstr("`tests'", ",", " ", .)))
    if "`tests_clean'" == "" local tests_clean "nullity homogeneity"
    local testranges_clean = trim("`testranges'")

    // Extract stored estimation artefacts
    local xvar = e(xvar)
    local method = e(method)
    local kernel = e(kernel)
    local method_code = cond("`method'" == "frechet", 1, 0)
    local kernel_type = cond("`kernel'" == "triangular", 1, ///
        cond("`kernel'" == "epanechnikov", 2, 3))
    local polynomial = e(polynomial)

    tempname HNUM
    matrix `HNUM' = e(bandwidth_num)

    local is_fuzzy = ("`e(tvar)'" != "")
    local h_den_scalar = cond(`is_fuzzy', e(bandwidth_den), 0)

    tempname ALPHA_PLUS ALPHA_MINUS W_PLUS W_MINUS E1 INT_PLUS INT_MINUS TAU_ORIG
    matrix `ALPHA_PLUS' = e(alpha_plus)
    matrix `ALPHA_MINUS' = e(alpha_minus)
    matrix `W_PLUS' = e(w_plus)
    matrix `W_MINUS' = e(w_minus)
    matrix `E1' = e(e1)
    matrix `INT_PLUS' = e(int_plus)
    matrix `INT_MINUS' = e(int_minus)
    matrix `TAU_ORIG' = e(b)

    tempname ALPHA_T_PLUS ALPHA_T_MINUS W_T_PLUS W_T_MINUS E2 DENOM
    if `is_fuzzy' {
        matrix `ALPHA_T_PLUS' = e(alpha_t_plus)
        matrix `ALPHA_T_MINUS' = e(alpha_t_minus)
        matrix `W_T_PLUS' = e(w_t_plus)
        matrix `W_T_MINUS' = e(w_t_minus)
        matrix `E2' = e(e2)
        scalar `DENOM' = e(denom)
    }
    else {
        scalar `DENOM' = 1
    }

    tempvar esample
    quietly gen byte `esample' = e(sample)

    tempvar x_centered
    quietly gen double `x_centered' = `xvar' - e(cutoff) if `esample'

    tempvar t_centered
    local tvarname ""
    if `is_fuzzy' {
        local tvar = e(tvar)
        local tvarname = "`t_centered'"
        quietly gen double `t_centered' = ``tvar'' - 0 if `esample'
    }

    tempname QUANT
    matrix `QUANT' = e(quantiles)
    local quantiles_str ""
    forvalues i = 1/`=colsof(`QUANT')' {
        local quantiles_str `quantiles_str' `= `QUANT'[1,`i']'
    }

    // Prepare containers for Mata call
    tempname bs_results pvals cb_lower cb_upper
    matrix `pvals' = J(1, 3, .)

    // Load xi_mat from CSV if provided
    local xi_mat_name ""
    if "`ximat'" != "" {
        tempname XI_MAT
        mata: st_matrix("`XI_MAT'", _r3d_load_csv("`ximat'"))
        local xi_mat_name "`XI_MAT'"
    }

    local e2_arg ""
    if `is_fuzzy' local e2_arg "`E2'"

    mata: r3d_bootstrap("`x_centered'", "", "`tvarname'", "`HNUM'", `h_den_scalar', ///
        `polynomial', `kernel_type', `reps', `level_val', ///
        "`ALPHA_PLUS'", "`ALPHA_MINUS'", "`W_PLUS'", "`W_MINUS'", ///
        "`E1'", "`e2_arg'", ///
        "`INT_PLUS'", "`INT_MINUS'", "`TAU_ORIG'", "`DENOM'", ///
        `method_code', "`quantiles_str'", "`tests_clean'", "`testranges_clean'", ///
        "`bs_results'", "`pvals'", "`esample'" , `is_fuzzy', "`xi_mat_name'")

    matrix `cb_lower' = r(cb_lower)
    matrix `cb_upper' = r(cb_upper)
    // pvals and bs_results are already populated via st_matrix() in Mata

    // Compute bootstrap standard errors from draws
    tempname boot_se
    mata: st_matrix("`boot_se'", sqrt(diagonal(variance(st_matrix("`bs_results'")))'))

    // Display summary
    di as text _n "Bootstrap inference for R3D"
    di as text "{hline 78}"
    di as text "Bootstrap replications: " as result %9.0f `reps'
    di as text "Confidence level: " as result %9.1f `level_val' "%"
    di as text "Tests: " as result "`tests_clean'"
    
    // Print table of effects
    local nq = e(nquantiles)
    di as text "{hline 78}"
    di as text "Quantile" _col(12) "Estimate" _col(24) "Boot SE" ///
        _col(36) "[`level_val'% Uniform CI]"
    di as text "{hline 78}"

    forvalues j = 1/`nq' {
        local q = `= `QUANT'[1,`j']'
        local tau = `= `TAU_ORIG'[1,`j']'
        local se = `= `boot_se'[1,`j']'
        local lo = `= `cb_lower'[1,`j']'
        local hi = `= `cb_upper'[1,`j']'
        di as result %8.3f `q' _col(12) %9.4f `tau' _col(24) %9.4f `se' ///
            _col(36) %9.4f `lo' _col(48) %9.4f `hi'
    }
    di as text "{hline 78}"

    // Show test p-values
    di as text "Hypothesis Tests"
    di as text "{hline 40}"
    di as text "Nullity (tau = 0):" _col(28) as result %9.4f `= `pvals'[1,1]'
    di as text "Homogeneity:" _col(28) as result %9.4f `= `pvals'[1,2]'
    if strpos("`tests_clean'", "gini") {
        di as text "Gini equality:" _col(28) as result %9.4f `= `pvals'[1,3]'
    }
    di as text "{hline 40}"

    quietly drop `x_centered'
    if `is_fuzzy' quietly drop `t_centered'

    // Store results in r()
    return scalar reps = `reps'
    return scalar level = `level_val'
    return matrix tau_boot = `bs_results'
    return matrix cb_lower = `cb_lower'
    return matrix cb_upper = `cb_upper'
    return matrix pvalues = `pvals'
    return matrix boot_se = `boot_se'
end
