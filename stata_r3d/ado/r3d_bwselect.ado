*! r3d_bwselect 2.0.0 - Bandwidth selection for R3D
*! Author: David Van Dijcke (Stata port updated)

program define r3d_bwselect, rclass
    version 18.0

    syntax varlist(min=2 numeric) [if] [in], ///
        Cutoff(real 0) ///
        [Method(string) ///
        FUZzy(varname numeric) ///
        POLynomial(integer 2) ///
        PILOT(integer 1) ///
        Kernel(string) ///
        Quantiles(numlist >0 <1 sort) ///
        NQuantiles(integer 99) ///
        Coverage ///
        WEights(varname numeric)]

    tokenize `varlist'
    local xvar `1'
    local yvars : list varlist - xvar

    marksample touse

    if "`method'" == "" local method "simple"
    if !inlist("`method'", "simple", "frechet") {
        di as error "method() must be either 'simple' or 'frechet'"
        exit 198
    }

    if "`kernel'" == "" local kernel "epanechnikov"
    if !inlist("`kernel'", "triangular", "epanechnikov", "uniform") {
        di as error "kernel() must be 'triangular', 'epanechnikov', or 'uniform'"
        exit 198
    }

    if `pilot' < 0 {
        di as error "pilot() must be non-negative"
        exit 198
    }

    local coverage_flag = cond("`coverage'" == "", 0, 1)

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
    local nq : word count `quantiles'

    tempvar x_centered
    quietly gen double `x_centered' = `xvar' - `cutoff' if `touse'

    local is_fuzzy = 0
    if "`fuzzy'" != "" {
        local is_fuzzy = 1
        tempvar t_var
        quietly gen double `t_var' = `fuzzy' if `touse'
    }

    tempname Qmat
    di as text "Computing empirical quantiles..."
    mata: r3d_compute_quantiles("`yvars'", "`quantiles'", "`touse'", "`Qmat'", "`weights'")

    local kernel_type = cond("`kernel'" == "triangular", 1, cond("`kernel'" == "epanechnikov", 2, 3))
    local tvarname = cond(`is_fuzzy', "`t_var'", "")

    di as text "Selecting bandwidth..."
    mata: r3d_bandwidth_select("`x_centered'", "`Qmat'", "`quantiles'", "`tvarname'", ///
        `polynomial', `pilot', `kernel_type', "`method'", "`touse'", `is_fuzzy', `coverage_flag', "`weights'")

    // INVARIANT: all r() scalars/matrices from r3d_bandwidth_select must be
    // captured before any subsequent Mata call that could post to r().
    // r3d_prepare_bandwidth_matrix is safe here because it only calls
    // st_matrix() and does not post to r().
    tempname HNUM PILOT_NUM
    matrix `HNUM' = r(h_num)
    matrix `PILOT_NUM' = r(pilot_num)
    local h_den_scalar = cond(`is_fuzzy', r(h_den), .)
    local pilot_den_scalar = cond(`is_fuzzy', r(pilot_den), .)

    tempname bw_rc
    mata: st_numscalar("`bw_rc'", r3d_prepare_bandwidth_matrix("`HNUM'", `nq', "`method'"))
    if scalar(`bw_rc') != 0 {
        di as error "Bandwidth selection failed"
        scalar drop `bw_rc'
        exit 498
    }
    scalar drop `bw_rc'
    tempname bw_min bw_max
    mata: st_numscalar("`bw_min'", min(st_matrix("`HNUM'")))
    mata: st_numscalar("`bw_max'", max(st_matrix("`HNUM'")))

    di as text "Bandwidth range (numerator): " as result %9.4f scalar(`bw_min') ///
        as text " to " as result %9.4f scalar(`bw_max')
    if `is_fuzzy' di as text "Bandwidth (denominator): " as result %9.4f `h_den_scalar'
    scalar drop `bw_min' `bw_max'

    // Return results
    tempname QGRID
    matrix `QGRID' = J(1, `nq', .)
    forvalues i = 1/`nq' {
        matrix `QGRID'[1,`i'] = `: word `i' of `quantiles''
    }

    return matrix h_num = `HNUM'
    return matrix pilot_num = `PILOT_NUM'
    if `is_fuzzy' {
        return scalar h_den = `h_den_scalar'
        return scalar pilot_den = `pilot_den_scalar'
    }
    return scalar method_code = cond("`method'" == "frechet", 1, 0)
    return scalar pilot = `pilot'
    return scalar polynomial = `polynomial'
    return scalar coverage = `coverage_flag'
    return matrix quantiles = `QGRID'
end
