/*
 * r3d_mata.mata - Core Mata functions for R3D package
 * Implements regression discontinuity with distributional outcomes
 *
 * Changes from original:
 *   - matastrict on (all declarations at function top)
 *   - r3d_locpoly: proper intercept weights e_1'(X'WX)^{-1} * basis * K
 *   - r3d_compute_quantiles: R type-7 interpolation
 *   - r3d_estimate_density: h = 1.06 * sd * n^(-1/5) (no IQR)
 */

version 18.0
mata:
mata set matastrict on

// ============================================================================
// UTILITY FUNCTIONS
// ============================================================================

// Kernel function implementations
real scalar r3d_kernel(real scalar u, real scalar kernel_type)
{
    if (kernel_type == 1) {  // Triangular
        return((abs(u) <= 1) ? 1 - abs(u) : 0)
    }
    else if (kernel_type == 2) {  // Epanechnikov
        return((abs(u) <= 1) ? 0.75 * (1 - u^2) : 0)
    }
    else {  // Uniform
        return((abs(u) <= 1) ? 0.5 : 0)
    }
}

// Compute variance of a vector
real scalar r3d_variance(real colvector x)
{
    real scalar n, mean_x
    n = rows(x)
    if (n <= 1) return(0)
    mean_x = mean(x)
    return(sum((x :- mean_x):^2) / (n - 1))
}

// ============================================================================
// QUANTILE COMPUTATION (R type-7 interpolation)
// ============================================================================

// Compute empirical quantiles from distribution data
// Uses R's type-7 quantile algorithm: h = (k-1)*p + 1, linear interpolation
void r3d_compute_quantiles(string scalar yvars, string scalar quantiles_str,
                          string scalar touse, string scalar matname,
                          string scalar weightvar)
{
    real matrix Y, Q, sorted
    real colvector q_grid, w, cumw
    real rowvector yi, wi
    real scalar i, j, n, nq, ny, k, idx_lo, idx_hi, hh, frac
    string rowvector varlist

    // Parse variable list
    varlist = tokens(yvars)
    ny = cols(varlist)

    // Get quantile grid
    q_grid = strtoreal(tokens(quantiles_str))'
    nq = rows(q_grid)

    // Get data
    st_view(Y, ., varlist, touse)
    n = rows(Y)

    // Get weights if specified
    if (weightvar != "") {
        st_view(w, ., weightvar, touse)
    }
    else {
        w = J(n, 1, 1)
    }

    // Initialize quantile matrix
    Q = J(n, nq, .)

    // Compute quantiles for each observation
    for (i = 1; i <= n; i++) {
        // Get non-missing values for this observation
        yi = select(Y[i,], Y[i,] :< .)
        k = cols(yi)

        if (k > 0) {
            if (weightvar != "") {
                // Weighted quantiles: use cumulative weight approach
                wi = J(1, k, w[i]/k)
                sorted = sort((yi' , wi'), 1)
                cumw = runningsum(sorted[,2])
                cumw = cumw / cumw[k]

                for (j = 1; j <= nq; j++) {
                    // Find index where cumulative weight exceeds quantile
                    for (idx_lo = 1; idx_lo <= k; idx_lo++) {
                        if (cumw[idx_lo] >= q_grid[j]) break
                    }
                    if (idx_lo > k) idx_lo = k
                    Q[i,j] = sorted[idx_lo,1]
                }
            }
            else {
                // Unweighted: R type-7 interpolation
                // Sort values
                sorted = sort(yi', 1)
                for (j = 1; j <= nq; j++) {
                    // R type-7: h = (k - 1) * p + 1
                    hh = (k - 1) * q_grid[j] + 1
                    idx_lo = floor(hh)
                    idx_hi = ceil(hh)
                    if (idx_lo < 1) idx_lo = 1
                    if (idx_hi > k) idx_hi = k
                    if (idx_lo == idx_hi) {
                        Q[i,j] = sorted[idx_lo]
                    }
                    else {
                        frac = hh - idx_lo
                        Q[i,j] = sorted[idx_lo] * (1 - frac) + sorted[idx_hi] * frac
                    }
                }
            }
        }
    }

    // Store result
    st_matrix(matname, Q)
}

// ============================================================================
// BANDWIDTH SELECTION
// ============================================================================

// Estimate density at cutoff: h = 1.06 * sd * n^(-1/5) (matches R)
real scalar r3d_estimate_density(real colvector X, real scalar kernel_type)
{
    real scalar n, h, f_hat, sd_x, i
    real colvector K

    n = rows(X)
    sd_x = sqrt(r3d_variance(X))
    if (sd_x <= 0) sd_x = 1
    h = 1.06 * sd_x * n^(-0.2)

    // Estimate density at zero
    K = J(n, 1, 0)
    for (i = 1; i <= n; i++) {
        K[i] = r3d_kernel(X[i] / h, kernel_type)
    }
    f_hat = mean(K) / h

    return(f_hat)
}

// Compute kernel matrices required for bandwidth selection (Appendix A.3)
void r3d_kernel_matrices(real scalar s, real scalar kernel_type,
                         real matrix Gamma_plus, real matrix Gamma_minus,
                         real colvector Lambda_plus, real colvector Lambda_minus,
                         real matrix Psi_plus, real matrix Psi_minus)
{
    real scalar m, step, i, j, dim, u, k_val
    real rowvector grid
    real colvector basis

    dim = s + 1
    m = 1000
    grid = rangen(-1, 1, m)
    step = 2.0 / (m - 1)

    Gamma_plus  = J(dim, dim, 0)
    Gamma_minus = J(dim, dim, 0)
    Psi_plus    = J(dim, dim, 0)
    Psi_minus   = J(dim, dim, 0)
    Lambda_plus  = J(dim, 1, 0)
    Lambda_minus = J(dim, 1, 0)

    for (i = 1; i <= m; i++) {
        u = grid[i]
        k_val = r3d_kernel(u, kernel_type)
        if (k_val <= 0) continue

        basis = J(dim, 1, 1)
        for (j = 2; j <= dim; j++) {
            basis[j] = basis[j-1] * u
        }

        if (u >= 0) {
            Gamma_plus = Gamma_plus + k_val * (basis * basis') * step
            Psi_plus   = Psi_plus   + (k_val^2) * (basis * basis') * step
            Lambda_plus = Lambda_plus + (u^(s+1)) * k_val * basis * step
        }
        else {
            Gamma_minus = Gamma_minus + k_val * (basis * basis') * step
            Psi_minus   = Psi_minus   + (k_val^2) * (basis * basis') * step
            Lambda_minus = Lambda_minus + (u^(s+1)) * k_val * basis * step
        }
    }
}

// Fit global polynomial and extract derivatives
real rowvector r3d_fit_global_poly(real colvector Y, real colvector X, real scalar order,
                                  real scalar resid_var)
{
    real matrix Xpoly, XtX
    real colvector coefs, resid
    real rowvector derivs
    real scalar i, n

    n = rows(X)

    // Build polynomial design matrix
    Xpoly = J(n, order+1, 1)
    for (i = 2; i <= order+1; i++) {
        Xpoly[,i] = X:^(i-1)
    }

    // Fit via OLS
    if (hasmissing(Y) | hasmissing(X)) {
        derivs = J(1, order, 0)
        resid_var = 1
        return(derivs)
    }

    XtX = cross(Xpoly, Xpoly)

    if (det(XtX) != 0) {
        coefs = lusolve(XtX, cross(Xpoly, Y))

        // Extract derivatives adjusted by factorial
        derivs = J(1, order, 0)
        for (i = 1; i <= order; i++) {
            if (i+1 <= rows(coefs)) {
                derivs[i] = coefs[i+1] * factorial(i)
            }
        }

        // Compute residual variance
        resid = Y - Xpoly * coefs
        resid_var = r3d_variance(resid)
    }
    else {
        derivs = J(1, order, 0)
        resid_var = 1
    }

    return(derivs)
}

// Core bandwidth selection routine mirroring the R implementation
void r3d_bandwidth_select(string scalar xvar, string scalar Qmat,
                          string scalar quantiles_str, string scalar tvar,
                          real scalar p, real scalar s, real scalar kernel_type,
                          string scalar method, string scalar touse,
                          real scalar is_fuzzy, real scalar coverage,
                          string scalar weightvar)
{
    real matrix X, Q, T
    real rowvector q_grid, derivs_plus, derivs_minus, deriv_T_plus, deriv_T_minus
    real scalar n, nq, sigma_X, f_X_hat, j, order_pilot, factorial_term
    real matrix Gamma_plus, Gamma_minus, Psi_plus, Psi_minus
    real colvector Lambda_plus, Lambda_minus
    real matrix I_s, Gamma_plus_inv, Gamma_minus_inv
    real colvector idx_plus, idx_minus
    real colvector pilot_derivs_plus, pilot_derivs_minus
    real colvector pilot_vars_plus, pilot_vars_minus
    real colvector Y_plus, X_plus, Y_minus, X_minus
    real scalar var_plus, var_minus
    real scalar pilot_deriv_T_plus, pilot_deriv_T_minus
    real scalar pilot_var_T_plus, pilot_var_T_minus
    real colvector T_plus, T_minus
    real scalar var_T_plus, var_T_minus
    real colvector e0, pilot_h_num
    real scalar bias_plus, bias_minus, C_1_0, C_1_0_prime
    real colvector temp_plus, temp_minus
    real scalar var_plus_j, var_minus_j, ratio
    real scalar pilot_h_den
    real scalar bias_T_plus, bias_T_minus, C_T_0, C_T_0_prime
    real colvector temp_plus_T, temp_minus_T
    real scalar var_T_plus_calc, var_T_minus_calc, ratio_T
    real matrix alpha_plus_pilot, alpha_minus_pilot, w_plus_pilot, w_minus_pilot
    real matrix h_pilot_vec
    real scalar h_scalar
    real scalar rc_plus, rc_minus
    real matrix alphaT_plus_pilot, alphaT_minus_pilot, wT_plus_pilot, wT_minus_pilot
    real matrix h_den_vec
    real scalar rc_tplus, rc_tminus
    real colvector B_plus, B_minus, V_plus, V_minus
    real scalar h_use, weight_sum, k
    real colvector idxp, idxm
    real colvector Xscaled, fitted_vals, resid_vals
    real matrix basis_mat
    real scalar B_plus_den, B_minus_den, V_plus_den, V_minus_den
    real colvector idxp_den, idxm_den
    real scalar weight_sum_den, k_den
    real colvector h_star_num
    real scalar bias_diff, var_sum
    real scalar dq, A_s, B_s, h_val
    real scalar h_star_den, bias_diff_den, var_sum_den, ratio_den
    real scalar shrink

    st_view(X, ., xvar, touse)
    Q = st_matrix(Qmat)

    if (quantiles_str == "") {
        errprintf("Quantile grid not supplied to bandwidth selector\n")
        error(498)
    }

    q_grid = strtoreal(tokens(quantiles_str))
    if (cols(q_grid) != cols(Q)) {
        errprintf("Quantile grid length does not match Q matrix\n")
        error(498)
    }

    n = rows(X)
    nq = cols(Q)

    if (is_fuzzy) {
        st_view(T, ., tvar, touse)
    }

    sigma_X = sqrt(r3d_variance(X))
    if (sigma_X <= 0) sigma_X = 1
    f_X_hat = r3d_estimate_density(X, kernel_type)
    if (f_X_hat <= 0) f_X_hat = 1

    r3d_kernel_matrices(s, kernel_type, Gamma_plus, Gamma_minus,
        Lambda_plus, Lambda_minus, Psi_plus, Psi_minus)

    I_s = I(s+1)
    if (abs(det(Gamma_plus)) < 1e-10)  Gamma_plus  = Gamma_plus  + 1e-8 * I_s
    if (abs(det(Gamma_minus)) < 1e-10) Gamma_minus = Gamma_minus + 1e-8 * I_s

    Gamma_plus_inv = invsym(Gamma_plus)
    Gamma_minus_inv = invsym(Gamma_minus)

    idx_plus = (X :>= 0)
    idx_minus = (X :< 0)

    pilot_derivs_plus = J(nq, 1, 0)
    pilot_derivs_minus = J(nq, 1, 0)
    pilot_vars_plus = J(nq, 1, 1)
    pilot_vars_minus = J(nq, 1, 1)

    order_pilot = s + 1

    for (j = 1; j <= nq; j++) {
        Y_plus = select(Q[,j], idx_plus)
        X_plus = select(X, idx_plus)
        if (rows(Y_plus) > order_pilot) {
            derivs_plus = r3d_fit_global_poly(Y_plus, X_plus, order_pilot, var_plus)
            if (cols(derivs_plus) >= order_pilot) pilot_derivs_plus[j] = derivs_plus[order_pilot]
            pilot_vars_plus[j] = var_plus
        }

        Y_minus = select(Q[,j], idx_minus)
        X_minus = select(X, idx_minus)
        if (rows(Y_minus) > order_pilot) {
            derivs_minus = r3d_fit_global_poly(Y_minus, X_minus, order_pilot, var_minus)
            if (cols(derivs_minus) >= order_pilot) pilot_derivs_minus[j] = derivs_minus[order_pilot]
            pilot_vars_minus[j] = var_minus
        }
    }

    pilot_deriv_T_plus = 0
    pilot_deriv_T_minus = 0
    pilot_var_T_plus = 1
    pilot_var_T_minus = 1

    if (is_fuzzy) {
        T_plus = select(T, idx_plus)
        if (rows(T_plus) > order_pilot) {
            deriv_T_plus = r3d_fit_global_poly(T_plus, X_plus, order_pilot, var_T_plus)
            if (cols(deriv_T_plus) >= order_pilot) pilot_deriv_T_plus = deriv_T_plus[order_pilot]
            pilot_var_T_plus = var_T_plus
        }

        T_minus = select(T, idx_minus)
        if (rows(T_minus) > order_pilot) {
            deriv_T_minus = r3d_fit_global_poly(T_minus, X_minus, order_pilot, var_T_minus)
            if (cols(deriv_T_minus) >= order_pilot) pilot_deriv_T_minus = deriv_T_minus[order_pilot]
            pilot_var_T_minus = var_T_minus
        }
    }

    e0 = J(s+1, 1, 0)
    e0[1] = 1

    pilot_h_num = J(nq, 1, 0)
    factorial_term = factorial(s + 1)

    for (j = 1; j <= nq; j++) {
        bias_plus = (e0' * (Gamma_plus_inv * Lambda_plus)) * (pilot_derivs_plus[j] / factorial_term)
        bias_minus = (e0' * (Gamma_minus_inv * Lambda_minus)) * (pilot_derivs_minus[j] / factorial_term)
        C_1_0 = bias_plus - bias_minus

        temp_plus = Gamma_plus_inv * e0
        temp_minus = Gamma_minus_inv * e0

        var_plus_j = pilot_vars_plus[j] * (temp_plus' * (Psi_plus * temp_plus))
        var_minus_j = pilot_vars_minus[j] * (temp_minus' * (Psi_minus * temp_minus))
        C_1_0_prime = (var_plus_j + var_minus_j) / f_X_hat

        if (abs(C_1_0) > 1e-14 && C_1_0_prime > 0) {
            ratio = C_1_0_prime / (2 * (s + 1) * (C_1_0^2))
            pilot_h_num[j] = ratio^(1 / (2 * s + 3)) * n^(-1 / (2 * s + 3))
        }

        if (pilot_h_num[j] <= 0) {
            pilot_h_num[j] = 1.06 * sigma_X * n^(-1 / (2 * s + 1))
        }
    }

    pilot_h_den = .
    if (is_fuzzy) {
        bias_T_plus = (e0' * (Gamma_plus_inv * Lambda_plus)) * (pilot_deriv_T_plus / factorial_term)
        bias_T_minus = (e0' * (Gamma_minus_inv * Lambda_minus)) * (pilot_deriv_T_minus / factorial_term)
        C_T_0 = bias_T_plus - bias_T_minus

        temp_plus_T = Gamma_plus_inv * e0
        temp_minus_T = Gamma_minus_inv * e0

        var_T_plus_calc = pilot_var_T_plus * (temp_plus_T' * (Psi_plus * temp_plus_T))
        var_T_minus_calc = pilot_var_T_minus * (temp_minus_T' * (Psi_minus * temp_minus_T))
        C_T_0_prime = (var_T_plus_calc + var_T_minus_calc) / f_X_hat

        if (abs(C_T_0) > 1e-14 && C_T_0_prime > 0) {
            ratio_T = C_T_0_prime / (2 * (s + 1) * (C_T_0^2))
            pilot_h_den = ratio_T^(1 / (2 * s + 3)) * n^(-1 / (2 * s + 3))
        }

        if (!(pilot_h_den > 0)) {
            pilot_h_den = 1.06 * sigma_X * n^(-1 / (2 * s + 1))
        }
    }

    h_pilot_vec = pilot_h_num

    if (method == "frechet") {
        h_scalar = mean(h_pilot_vec)
        if (!(h_scalar > 0)) h_scalar = 1.06 * sigma_X * n^(-1 / (2 * s + 1))
        h_pilot_vec = J(nq, 1, h_scalar)
    }

    if (r3d_plugin_available()) {
        rc_plus = r3d_plugin_locweights(X, Q, h_pilot_vec, s, kernel_type, 1, alpha_plus_pilot, w_plus_pilot)
        rc_minus = r3d_plugin_locweights(X, Q, h_pilot_vec, s, kernel_type, 0, alpha_minus_pilot, w_minus_pilot)
        if (rc_plus != 0 | rc_minus != 0) {
            printf("{txt}Note: plugin failed (rc=%g), falling back to Mata\n", max((rc_plus, rc_minus)))
            st_global("R3D_USE_PLUGIN", "0")
            st_global("R3D_PLUGIN_NAME", "")
            r3d_locpoly(X, Q, h_pilot_vec, s, kernel_type, 1, alpha_plus_pilot, w_plus_pilot)
            r3d_locpoly(X, Q, h_pilot_vec, s, kernel_type, 0, alpha_minus_pilot, w_minus_pilot)
        }
    }
    else {
        r3d_locpoly(X, Q, h_pilot_vec, s, kernel_type, 1, alpha_plus_pilot, w_plus_pilot)
        r3d_locpoly(X, Q, h_pilot_vec, s, kernel_type, 0, alpha_minus_pilot, w_minus_pilot)
    }

    if (is_fuzzy) {
        h_den_vec = J(1, 1, pilot_h_den)

        if (r3d_plugin_available()) {
            rc_tplus = r3d_plugin_locweights(X, T, h_den_vec, s, kernel_type, 1, alphaT_plus_pilot, wT_plus_pilot)
            rc_tminus = r3d_plugin_locweights(X, T, h_den_vec, s, kernel_type, 0, alphaT_minus_pilot, wT_minus_pilot)
            if (rc_tplus != 0 | rc_tminus != 0) {
                printf("{txt}Note: plugin failed (rc=%g), falling back to Mata\n", max((rc_tplus, rc_tminus)))
                st_global("R3D_USE_PLUGIN", "0")
                st_global("R3D_PLUGIN_NAME", "")
                r3d_locpoly(X, T, h_den_vec, s, kernel_type, 1, alphaT_plus_pilot, wT_plus_pilot)
                r3d_locpoly(X, T, h_den_vec, s, kernel_type, 0, alphaT_minus_pilot, wT_minus_pilot)
            }
        }
        else {
            r3d_locpoly(X, T, h_den_vec, s, kernel_type, 1, alphaT_plus_pilot, wT_plus_pilot)
            r3d_locpoly(X, T, h_den_vec, s, kernel_type, 0, alphaT_minus_pilot, wT_minus_pilot)
        }
    }

    B_plus = J(nq, 1, 0)
    B_minus = J(nq, 1, 0)
    V_plus = J(nq, 1, 0)
    V_minus = J(nq, 1, 0)

    for (j = 1; j <= nq; j++) {
        B_plus[j] = alpha_plus_pilot[s+1, j] / factorial_term
        B_minus[j] = alpha_minus_pilot[s+1, j] / factorial_term

        h_use = h_pilot_vec[j,1]
        if (!(h_use > 0)) continue

        idxp = (w_plus_pilot[,j] :> 0)
        idxm = (w_minus_pilot[,j] :> 0)

        if (sum(idxp) > 0) {
            Xscaled = select(X, idxp) / h_use
            basis_mat = J(rows(Xscaled), s+1, 1)
            for (k = 2; k <= s+1; k++) basis_mat[,k] = basis_mat[,k-1] :* Xscaled
            fitted_vals = basis_mat * alpha_plus_pilot[,j]
            resid_vals = select(Q[,j], idxp) - fitted_vals
            weight_sum = sum(select(w_plus_pilot[,j], idxp))
            if (weight_sum > 0) V_plus[j] = sum(select(w_plus_pilot[,j], idxp) :* (resid_vals:^2)) / weight_sum
        }

        if (sum(idxm) > 0) {
            Xscaled = select(X, idxm) / h_use
            basis_mat = J(rows(Xscaled), s+1, 1)
            for (k = 2; k <= s+1; k++) basis_mat[,k] = basis_mat[,k-1] :* Xscaled
            fitted_vals = basis_mat * alpha_minus_pilot[,j]
            resid_vals = select(Q[,j], idxm) - fitted_vals
            weight_sum = sum(select(w_minus_pilot[,j], idxm))
            if (weight_sum > 0) V_minus[j] = sum(select(w_minus_pilot[,j], idxm) :* (resid_vals:^2)) / weight_sum
        }
    }

    B_plus_den = 0
    B_minus_den = 0
    V_plus_den = 0
    V_minus_den = 0

    if (is_fuzzy) {
        B_plus_den = alphaT_plus_pilot[s+1,1] / factorial_term
        B_minus_den = alphaT_minus_pilot[s+1,1] / factorial_term

        idxp_den = (wT_plus_pilot :> 0)
        idxm_den = (wT_minus_pilot :> 0)

        if (sum(idxp_den) > 0) {
            Xscaled = select(X, idxp_den) / pilot_h_den
            basis_mat = J(rows(Xscaled), s+1, 1)
            for (k_den = 2; k_den <= s+1; k_den++) basis_mat[,k_den] = basis_mat[,k_den-1] :* Xscaled
            fitted_vals = basis_mat * alphaT_plus_pilot
            resid_vals = select(T, idxp_den) - fitted_vals
            weight_sum_den = sum(select(wT_plus_pilot, idxp_den))
            if (weight_sum_den > 0) V_plus_den = sum(select(wT_plus_pilot, idxp_den) :* (resid_vals:^2)) / weight_sum_den
        }

        if (sum(idxm_den) > 0) {
            Xscaled = select(X, idxm_den) / pilot_h_den
            basis_mat = J(rows(Xscaled), s+1, 1)
            for (k_den = 2; k_den <= s+1; k_den++) basis_mat[,k_den] = basis_mat[,k_den-1] :* Xscaled
            fitted_vals = basis_mat * alphaT_minus_pilot
            resid_vals = select(T, idxm_den) - fitted_vals
            weight_sum_den = sum(select(wT_minus_pilot, idxm_den))
            if (weight_sum_den > 0) V_minus_den = sum(select(wT_minus_pilot, idxm_den) :* (resid_vals:^2)) / weight_sum_den
        }
    }

    h_star_num = J(nq, 1, 0)

    if (method == "simple") {
        for (j = 1; j <= nq; j++) {
            bias_diff = B_plus[j] - B_minus[j]
            var_sum = (V_plus[j] + V_minus[j]) / f_X_hat
            if (abs(bias_diff) > 1e-14 && var_sum > 0) {
                ratio = var_sum / (2 * (s + 1) * (bias_diff^2))
                h_star_num[j] = ratio^(1 / (2 * s + 3)) * n^(-1 / (2 * s + 3))
            }
            if (h_star_num[j] <= 0) h_star_num[j] = 1.06 * sigma_X * n^(-1 / (2 * s + 1))
        }
    }
    else {
        dq = (max(q_grid) - min(q_grid)) / (nq - 1)
        if (dq <= 0) dq = 1.0 / nq
        A_s = sum((B_plus - B_minus):^2) * dq
        B_s = sum((V_plus + V_minus)) * dq / f_X_hat
        if (A_s > 1e-14 && B_s > 0) {
            ratio = B_s / (2 * (s + 1) * A_s)
            h_val = ratio^(1 / (2 * s + 3)) * n^(-1 / (2 * s + 3))
        }
        else {
            h_val = 1.06 * sigma_X * n^(-1 / (2 * s + 1))
        }
        h_star_num = J(nq, 1, h_val)
    }

    h_star_den = .
    if (is_fuzzy) {
        bias_diff_den = B_plus_den - B_minus_den
        var_sum_den = (V_plus_den + V_minus_den) / f_X_hat
        if (abs(bias_diff_den) > 1e-14 && var_sum_den > 0) {
            ratio_den = var_sum_den / (2 * (s + 1) * (bias_diff_den^2))
            h_star_den = ratio_den^(1 / (2 * s + 3)) * n^(-1 / (2 * s + 3))
        }
        if (!(h_star_den > 0)) h_star_den = 1.06 * sigma_X * n^(-1 / (2 * s + 1))
    }

    if (coverage) {
        shrink = n^(-s / ((2 * s + 3) * (s + 3)))
        h_star_num = h_star_num :* shrink
        if (is_fuzzy) h_star_den = h_star_den * shrink
    }

    st_matrix("r(h_num)", h_star_num)
    st_matrix("r(pilot_num)", pilot_h_num)
    st_matrix("r(B_plus)", B_plus)
    st_matrix("r(B_minus)", B_minus)
    st_matrix("r(V_plus)", V_plus)
    st_matrix("r(V_minus)", V_minus)
    st_numscalar("r(f_X_hat)", f_X_hat)

    if (is_fuzzy) {
        st_numscalar("r(h_den)", h_star_den)
        st_numscalar("r(pilot_den)", pilot_h_den)
        st_numscalar("r(B_plus_den)", B_plus_den)
        st_numscalar("r(B_minus_den)", B_minus_den)
        st_numscalar("r(V_plus_den)", V_plus_den)
        st_numscalar("r(V_minus_den)", V_minus_den)
    }
}

// Validate and broadcast user-supplied bandwidth vector
real scalar r3d_prepare_bandwidth_matrix(string scalar matname, real scalar nq,
                                        string scalar method)
{
    real matrix H
    real scalar unique_count

    H = st_matrix(matname)
    if (rows(H) == 1 & cols(H) == 1) {
        H = J(nq, 1, H[1,1])
    }
    else if (rows(H) == 1 & cols(H) == nq) {
        H = H'
    }
    else if (!(rows(H) == nq & cols(H) == 1)) {
        return(198)
    }

    if (method == "frechet") {
        unique_count = rows( uniqrows(H) )
        if (unique_count > 1) {
            return(198)
        }
        H = J(nq, 1, H[1,1])
    }

    st_matrix(matname, H)
    return(0)
}

// ============================================================================
// LOCAL POLYNOMIAL REGRESSION
// ============================================================================

// Core local polynomial regression with proper intercept weights
// Computes: alpha = (X'WX)^{-1} X'WY (coefficients)
//           weights = e_1' (X'WX)^{-1} * basis_i * K(u_i)  (intercept weights)
// This matches the Fortran plugin (locweights.f90 step 5)
void r3d_locpoly(real colvector X, real matrix Y, real matrix h_mat,
                real scalar p, real scalar kernel_type, real scalar side,
                real matrix alpha, real matrix weights)
{
    real scalar n, nq, i, j, k
    real matrix XWX, XWY
    real colvector basis, e1, first_row_inv
    real scalar u, w_i, w_int
    real scalar h_j, rcond_val
    real colvector IPIV

    n = rows(X)
    nq = cols(Y)

    // Normalise bandwidth input to column vector of length nq
    if (rows(h_mat) == 1 & cols(h_mat) == 1) {
        h_mat = J(nq, 1, h_mat[1,1])
    }
    else if (rows(h_mat) == 1 & cols(h_mat) == nq) {
        h_mat = h_mat'
    }
    else if (!(rows(h_mat) == nq & cols(h_mat) == 1)) {
        errprintf("r3d_locpoly(): bandwidth vector has incorrect dimensions\n")
        return
    }

    // Initialize outputs
    alpha = J(p+1, nq, 0)
    weights = J(n, nq, 0)

    // e_1 vector for extracting first row of inverse
    e1 = J(p+1, 1, 0)
    e1[1] = 1

    for (j = 1; j <= nq; j++) {
        h_j = h_mat[j,1]
        if (h_j <= 0) continue

        // Initialize matrices
        XWX = J(p+1, p+1, 0)
        XWY = J(p+1, 1, 0)

        // Compute weighted matrices
        for (i = 1; i <= n; i++) {
            // Skip observations on wrong side
            if (side == 1 & X[i] < 0) continue
            if (side == 0 & X[i] >= 0) continue

            // Skip if Y is missing
            if (Y[i,j] >= .) continue

            // Compute kernel weight
            u = X[i] / h_j
            w_i = r3d_kernel(u, kernel_type)

            if (w_i > 0) {
                // Compute polynomial basis
                basis = J(p+1, 1, 1)
                for (k = 2; k <= p+1; k++) {
                    basis[k] = basis[k-1] * (X[i]/h_j)
                }

                // Accumulate
                XWX = XWX + w_i * basis * basis'
                XWY = XWY + w_i * basis * Y[i,j]
            }
        }

        // Solve for coefficients (use reciprocal condition number for robust singularity check)
        rcond_val = 1/cond(XWX)
        if (rcond_val > epsilon(1)) {
            alpha[,j] = lusolve(XWX, XWY)

            // Compute first row of (X'WX)^{-1} = solve(XWX, e1)
            first_row_inv = lusolve(XWX, e1)

            // Compute intercept weights for each observation
            // w_int_i = (first_row_inv' * basis_i) * K(u_i)
            // This matches Fortran locweights.f90 lines 131-192
            for (i = 1; i <= n; i++) {
                if (side == 1 & X[i] < 0) continue
                if (side == 0 & X[i] >= 0) continue
                if (Y[i,j] >= .) continue

                u = X[i] / h_j
                w_i = r3d_kernel(u, kernel_type)

                if (w_i > 0) {
                    basis = J(p+1, 1, 1)
                    for (k = 2; k <= p+1; k++) {
                        basis[k] = basis[k-1] * (X[i]/h_j)
                    }

                    // Dot product of first row of inverse with basis, times kernel weight
                    w_int = 0
                    for (k = 1; k <= p+1; k++) {
                        w_int = w_int + first_row_inv[k] * basis[k]
                    }
                    weights[i,j] = w_int * w_i
                }
            }
        }
        else {
            alpha[,j] = J(p+1, 1, 0)
        }
    }
}

// ============================================================================
// ISOTONIC REGRESSION
// ============================================================================

// Pool adjacent violators algorithm for isotonic regression
real colvector r3d_isotonic(real colvector y)
{
    real scalar n, i, j, sum_y, len, avg, k
    real colvector y_iso

    n = rows(y)
    if (n <= 1) return(y)

    y_iso = y

    // Pool adjacent violators algorithm
    i = 1
    while (i < n) {
        if (y_iso[i] > y_iso[i+1]) {
            // Find violating block
            j = i + 1
            sum_y = y_iso[i] + y_iso[i+1]
            len = 2

            // Extend block while violations continue
            while (j < n & sum_y/len > y_iso[j+1]) {
                j++
                sum_y = sum_y + y_iso[j]
                len++
            }

            // Pool the block
            avg = sum_y / len
            for (k = i; k <= j; k++) {
                y_iso[k] = avg
            }

            // Backtrack to check previous values
            if (i > 1) {
                i--
            }
        }
        else {
            i++
        }
    }

    return(y_iso)
}

// ============================================================================
// MAIN ESTIMATION FUNCTIONS
// ============================================================================

// Core estimator covering both simple and Frechet methods
void r3d_estimate_core(string scalar method,
                      string scalar xvar, string scalar Qmat, string scalar tvar,
                      string scalar hnum_matname, real scalar h_den,
                      real scalar p, real scalar kernel_type,
                      string scalar tau_name, string scalar se_name,
                      string scalar alpha_plus_name, string scalar alpha_minus_name,
                      string scalar w_plus_name, string scalar w_minus_name,
                      string scalar e1_name, string scalar e2_name,
                      string scalar int_plus_name, string scalar int_minus_name,
                      string scalar alpha_t_plus_name, string scalar alpha_t_minus_name,
                      string scalar w_t_plus_name, string scalar w_t_minus_name,
                      string scalar denom_name,
                      string scalar touse, real scalar is_fuzzy)
{
    real matrix X, Q, T
    real scalar n, nq, j, k, i
    real matrix h_num_vec
    real matrix alpha_plus, alpha_minus, w_plus, w_minus
    real scalar rc_plus, rc_minus
    real matrix alpha_plus_t, alpha_minus_t, w_plus_t, w_minus_t
    real matrix h_den_vec
    real scalar rc_tplus, rc_tminus
    real matrix e2
    real scalar denominator
    real matrix Eplus, Eminus, Eplus_final, Eminus_final
    real scalar h_j
    real colvector idxp, idxm
    real colvector Xscaled, fitted_vals
    real matrix basis_mat
    real colvector rowfit
    real matrix e1
    real rowvector int_plus, int_minus, tau, se
    real scalar n_plus, n_minus, var_plus, var_minus

    st_view(X, ., xvar, touse)
    Q = st_matrix(Qmat)

    n = rows(X)
    nq = cols(Q)

    h_num_vec = st_matrix(hnum_matname)

    if (rows(h_num_vec) == 1 & cols(h_num_vec) == 1) {
        h_num_vec = J(nq, 1, h_num_vec[1,1])
    }
    else if (rows(h_num_vec) == 1 & cols(h_num_vec) == nq) {
        h_num_vec = h_num_vec'
    }

    if (rows(h_num_vec) != nq | cols(h_num_vec) != 1) {
        errprintf("Bandwidth vector dimension mismatch in r3d_estimate_core()\n")
        return
    }

    if (r3d_plugin_available()) {
        rc_plus = r3d_plugin_locweights(X, Q, h_num_vec, p, kernel_type, 1, alpha_plus, w_plus)
        rc_minus = r3d_plugin_locweights(X, Q, h_num_vec, p, kernel_type, 0, alpha_minus, w_minus)
        if (rc_plus != 0 | rc_minus != 0) {
            printf("{txt}Note: plugin failed (rc=%g), falling back to Mata\n", max((rc_plus, rc_minus)))
            st_global("R3D_USE_PLUGIN", "0")
            st_global("R3D_PLUGIN_NAME", "")
            r3d_locpoly(X, Q, h_num_vec, p, kernel_type, 1, alpha_plus, w_plus)
            r3d_locpoly(X, Q, h_num_vec, p, kernel_type, 0, alpha_minus, w_minus)
        }
    }
    else {
        r3d_locpoly(X, Q, h_num_vec, p, kernel_type, 1, alpha_plus, w_plus)
        r3d_locpoly(X, Q, h_num_vec, p, kernel_type, 0, alpha_minus, w_minus)
    }

    denominator = 1
    e2 = J(n, nq, 0)

    if (is_fuzzy) {
        st_view(T, ., tvar, touse)

        h_den_vec = J(1, 1, h_den)

        if (r3d_plugin_available()) {
            rc_tplus = r3d_plugin_locweights(X, T, h_den_vec, p, kernel_type, 1, alpha_plus_t, w_plus_t)
            rc_tminus = r3d_plugin_locweights(X, T, h_den_vec, p, kernel_type, 0, alpha_minus_t, w_minus_t)
            if (rc_tplus != 0 | rc_tminus != 0) {
                printf("{txt}Note: plugin failed (rc=%g), falling back to Mata\n", max((rc_tplus, rc_tminus)))
                st_global("R3D_USE_PLUGIN", "0")
                st_global("R3D_PLUGIN_NAME", "")
                r3d_locpoly(X, T, h_den_vec, p, kernel_type, 1, alpha_plus_t, w_plus_t)
                r3d_locpoly(X, T, h_den_vec, p, kernel_type, 0, alpha_minus_t, w_minus_t)
            }
        }
        else {
            r3d_locpoly(X, T, h_den_vec, p, kernel_type, 1, alpha_plus_t, w_plus_t)
            r3d_locpoly(X, T, h_den_vec, p, kernel_type, 0, alpha_minus_t, w_minus_t)
        }

        denominator = alpha_plus_t[1,1] - alpha_minus_t[1,1]
    }

    // Predicted values on each side
    Eplus = J(n, nq, .)
    Eminus = J(n, nq, .)

    for (j = 1; j <= nq; j++) {
        h_j = h_num_vec[j,1]
        if (!(h_j > 0)) continue

        idxp = (w_plus[,j] :> 0)
        idxm = (w_minus[,j] :> 0)

        if (sum(idxp) > 0) {
            Xscaled = select(X, idxp) / h_j
            basis_mat = J(rows(Xscaled), p+1, 1)
            for (k = 2; k <= p+1; k++) {
                basis_mat[,k] = basis_mat[,k-1] :* Xscaled
            }
            fitted_vals = basis_mat * alpha_plus[,j]
            Eplus[selectindex(idxp), j] = fitted_vals
        }

        if (sum(idxm) > 0) {
            Xscaled = select(X, idxm) / h_j
            basis_mat = J(rows(Xscaled), p+1, 1)
            for (k = 2; k <= p+1; k++) {
                basis_mat[,k] = basis_mat[,k-1] :* Xscaled
            }
            fitted_vals = basis_mat * alpha_minus[,j]
            Eminus[selectindex(idxm), j] = fitted_vals
        }
    }

    Eplus_final = Eplus
    Eminus_final = Eminus

    if (method == "frechet") {
        for (i = 1; i <= n; i++) {
            if (sum(w_plus[i,] :> 0)) {
                rowfit = Eplus[i,]'
                rowfit = r3d_isotonic(rowfit)
                Eplus_final[i,] = rowfit'
            }
            if (sum(w_minus[i,] :> 0)) {
                rowfit = Eminus[i,]'
                rowfit = r3d_isotonic(rowfit)
                Eminus_final[i,] = rowfit'
            }
        }
    }

    // Residuals for outcome
    e1 = J(n, nq, 0)
    for (j = 1; j <= nq; j++) {
        idxp = (w_plus[,j] :> 0)
        idxm = (w_minus[,j] :> 0)
        if (sum(idxp) > 0) {
            e1[selectindex(idxp), j] = select(Q[,j], idxp) - Eplus_final[selectindex(idxp), j]
        }
        if (sum(idxm) > 0) {
            e1[selectindex(idxm), j] = select(Q[,j], idxm) - Eminus_final[selectindex(idxm), j]
        }
    }

    // Residuals for treatment (replicated across quantiles)
    if (is_fuzzy) {
        idxp = (w_plus_t :> 0)
        idxm = (w_minus_t :> 0)
        if (sum(idxp) > 0) {
            Xscaled = select(X, idxp) / h_den
            basis_mat = J(rows(Xscaled), p+1, 1)
            for (k = 2; k <= p+1; k++) {
                basis_mat[,k] = basis_mat[,k-1] :* Xscaled
            }
            fitted_vals = basis_mat * alpha_plus_t
            e2[selectindex(idxp),] = (select(T, idxp) - fitted_vals) * J(1, nq, 1)
        }
        if (sum(idxm) > 0) {
            Xscaled = select(X, idxm) / h_den
            basis_mat = J(rows(Xscaled), p+1, 1)
            for (k = 2; k <= p+1; k++) {
                basis_mat[,k] = basis_mat[,k-1] :* Xscaled
            }
            fitted_vals = basis_mat * alpha_minus_t
            e2[selectindex(idxm),] = (select(T, idxm) - fitted_vals) * J(1, nq, 1)
        }
    }

    // Intercepts and rearrangement/isotonic adjustments
    int_plus = alpha_plus[1,]
    int_minus = alpha_minus[1,]

    if (method == "simple") {
        int_plus = sort(int_plus', 1)'
        int_minus = sort(int_minus', 1)'
    }
    else {
        int_plus = r3d_isotonic(int_plus')'
        int_minus = r3d_isotonic(int_minus')'
    }
    alpha_plus[1,] = int_plus
    alpha_minus[1,] = int_minus

    // Quantile treatment effect
    tau = int_plus - int_minus

    if (is_fuzzy) {
        if (abs(denominator) > 1e-8) {
            tau = tau / denominator
        }
        else {
            tau = J(1, nq, .)
        }
    }

    // Standard errors (approximate)
    se = J(1, nq, .)
    for (j = 1; j <= nq; j++) {
        idxp = (w_plus[,j] :> 0)
        idxm = (w_minus[,j] :> 0)
        n_plus = sum(idxp)
        n_minus = sum(idxm)
        if (n_plus > p+1 & n_minus > p+1) {
            var_plus = r3d_variance(select(e1[,j], idxp))
            var_minus = r3d_variance(select(e1[,j], idxm))
            if (var_plus >= 0 & var_minus >= 0) {
                se[j] = sqrt(var_plus/n_plus + var_minus/n_minus)
                if (is_fuzzy & abs(denominator) > 1e-8) {
                    se[j] = se[j] / abs(denominator)
                }
            }
        }
    }

    // Store results
    st_matrix(tau_name, tau)
    st_matrix(se_name, se)
    st_matrix(alpha_plus_name, alpha_plus)
    st_matrix(alpha_minus_name, alpha_minus)
    st_matrix(w_plus_name, w_plus)
    st_matrix(w_minus_name, w_minus)
    st_matrix(e1_name, e1)
    st_matrix(int_plus_name, int_plus)
    st_matrix(int_minus_name, int_minus)
    st_numscalar(denom_name, denominator)

    if (is_fuzzy) {
        st_matrix(e2_name, e2)
        if (alpha_t_plus_name != "") st_matrix(alpha_t_plus_name, alpha_plus_t)
        if (alpha_t_minus_name != "") st_matrix(alpha_t_minus_name, alpha_minus_t)
        if (w_t_plus_name != "") st_matrix(w_t_plus_name, w_plus_t)
        if (w_t_minus_name != "") st_matrix(w_t_minus_name, w_minus_t)
    }
    else {
        st_matrix(e2_name, J(n, nq, 0))
        if (alpha_t_plus_name != "") st_matrix(alpha_t_plus_name, J(p+1, 1, .))
        if (alpha_t_minus_name != "") st_matrix(alpha_t_minus_name, J(p+1, 1, .))
        if (w_t_plus_name != "") st_matrix(w_t_plus_name, J(n, 1, 0))
        if (w_t_minus_name != "") st_matrix(w_t_minus_name, J(n, 1, 0))
    }
}

// Simple method wrapper
void r3d_simple(string scalar xvar, string scalar Qmat, string scalar tvar,
               string scalar hnum_matname, real scalar h_den, real scalar p,
               real scalar kernel_type, string scalar tau_name, string scalar se_name,
               string scalar alpha_plus_name, string scalar alpha_minus_name,
               string scalar w_plus_name, string scalar w_minus_name,
               string scalar e1_name, string scalar e2_name,
               string scalar int_plus_name, string scalar int_minus_name,
               string scalar alpha_t_plus_name, string scalar alpha_t_minus_name,
               string scalar w_t_plus_name, string scalar w_t_minus_name,
               string scalar denom_name,
               string scalar touse, real scalar is_fuzzy)
{
    r3d_estimate_core("simple", xvar, Qmat, tvar, hnum_matname, h_den, p, kernel_type,
        tau_name, se_name, alpha_plus_name, alpha_minus_name,
        w_plus_name, w_minus_name, e1_name, e2_name,
        int_plus_name, int_minus_name,
        alpha_t_plus_name, alpha_t_minus_name,
        w_t_plus_name, w_t_minus_name,
        denom_name,
        touse, is_fuzzy);
}

// Frechet method wrapper
void r3d_frechet(string scalar xvar, string scalar Qmat, string scalar tvar,
                string scalar hnum_matname, real scalar h_den, real scalar p,
                real scalar kernel_type, string scalar tau_name, string scalar se_name,
                string scalar alpha_plus_name, string scalar alpha_minus_name,
                string scalar w_plus_name, string scalar w_minus_name,
                string scalar e1_name, string scalar e2_name,
                string scalar int_plus_name, string scalar int_minus_name,
                string scalar alpha_t_plus_name, string scalar alpha_t_minus_name,
                string scalar w_t_plus_name, string scalar w_t_minus_name,
                string scalar denom_name,
                string scalar touse, real scalar is_fuzzy)
{
    r3d_estimate_core("frechet", xvar, Qmat, tvar, hnum_matname, h_den, p, kernel_type,
        tau_name, se_name, alpha_plus_name, alpha_minus_name,
        w_plus_name, w_minus_name, e1_name, e2_name,
        int_plus_name, int_minus_name,
        alpha_t_plus_name, alpha_t_minus_name,
        w_t_plus_name, w_t_minus_name,
        denom_name,
        touse, is_fuzzy);
}

// ============================================================================
// BOOTSTRAP INFERENCE
// ============================================================================

// Multiplier bootstrap for uniform confidence bands and hypothesis tests
void r3d_bootstrap(string scalar xvar, string scalar Qmat, string scalar tvar,
                  string scalar hnum_matname, real scalar h_den, real scalar p,
                  real scalar kernel_type, real scalar B, real scalar level,
                  string scalar alpha_plus_name, string scalar alpha_minus_name,
                  string scalar w_plus_name, string scalar w_minus_name,
                  string scalar e1_name, string scalar e2_name,
                  string scalar int_plus_name, string scalar int_minus_name,
                  string scalar tau_name, string scalar denom_name,
                  real scalar method_code, string scalar quantiles_str,
                  string scalar tests_str, string scalar test_ranges_str,
                  string scalar bs_results_name, string scalar pvals_name,
                  string scalar touse, real scalar is_fuzzy,
                  | string scalar xi_mat_path)
{
    real matrix X
    real scalar n, nq, b, j, r
    real matrix h_num_vec
    real matrix alpha_plus, alpha_minus, w_plus, w_minus
    real matrix e1_mat, e2_mat
    real rowvector int_plus, int_minus, tau_orig, num_diff
    real scalar denominator
    real matrix e1_w_plus, e1_w_minus, e2_w_plus, e2_w_minus
    real scalar is_vector_h_num
    string rowvector tests_tokens, filtered
    real scalar do_nullity, do_homogeneity, do_gini
    real matrix tau_boot, nu_plus_store, nu_minus_store
    real colvector xi
    real rowvector plus_sums, minus_sums, out_sharp
    real scalar denom_term
    real rowvector top
    real colvector supvals, sup_sorted
    real matrix xi_mat_data
    real scalar use_xi_mat
    real scalar prob, alpha_level, index_c, cval
    real rowvector cb_lower, cb_upper, q_grid
    real matrix range_pairs
    real rowvector range_tokens
    real scalar n_pairs, idx
    real scalar range_lo, range_hi
    real colvector range_idx
    real matrix pvals
    real scalar test_stat_null
    real colvector supvals_null, sup_null_sorted
    real scalar idx_null, crit_null, p_null
    real scalar test_stat_homo
    real rowvector range_tau, draw
    real scalar mbar, draw_mean
    real colvector supvals_homo, sup_homo_sorted
    real scalar idx_homo, crit_homo, p_homo
    real scalar gini_above, gini_below, gini_diff, test_stat_gini
    real colvector gini_boot, gini_sorted
    real rowvector boot_plus, boot_minus
    real scalar idx_gini, crit_gini, p_gini

    st_view(X, ., xvar, touse)
    n = rows(X)

    // Load pre-generated xi matrix if path provided
    use_xi_mat = 0
    xi_mat_data = J(0, 0, .)
    if (args() >= 28 & xi_mat_path != "") {
        xi_mat_data = st_matrix(xi_mat_path)
        if (rows(xi_mat_data) == n) {
            use_xi_mat = 1
            B = cols(xi_mat_data)
        }
    }

    h_num_vec = st_matrix(hnum_matname)
    if (rows(h_num_vec) == 1 & cols(h_num_vec) == 1) {
        h_num_vec = J(cols(st_matrix(tau_name)), 1, h_num_vec[1,1])
    }
    else if (rows(h_num_vec) == 1 & cols(h_num_vec) > 1) {
        h_num_vec = h_num_vec'
    }

    alpha_plus = st_matrix(alpha_plus_name)
    alpha_minus = st_matrix(alpha_minus_name)
    w_plus = st_matrix(w_plus_name)
    w_minus = st_matrix(w_minus_name)

    e1_mat = st_matrix(e1_name)

    int_plus = st_matrix(int_plus_name)
    int_minus = st_matrix(int_minus_name)
    tau_orig = st_matrix(tau_name)

    denominator = st_numscalar(denom_name)

    nq = cols(alpha_plus)

    e1_w_plus = e1_mat :* w_plus
    e1_w_minus = e1_mat :* w_minus

    if (is_fuzzy) {
        e2_mat = st_matrix(e2_name)
        e2_w_plus = e2_mat :* w_plus
        e2_w_minus = e2_mat :* w_minus
    }

    num_diff = int_plus - int_minus

    is_vector_h_num = (method_code == 0)

    // Parse tests
    tests_str = strlower(strtrim(subinstr(tests_str, ",", " ", .)))
    tests_tokens = tokens(tests_str)

    filtered = J(1, 0, "")
    for (j = 1; j <= cols(tests_tokens); j++) {
        if (tests_tokens[j] != "" & tests_tokens[j] != "none") {
            filtered = filtered , tests_tokens[j]
        }
    }
    if (cols(filtered) > 0) {
        tests_tokens = filtered
    }
    else {
        tests_tokens = J(1, 0, "")
    }

    do_nullity = 0
    do_homogeneity = 0
    do_gini = 0
    if (cols(tests_tokens) > 0) {
        do_nullity = any(tests_tokens :== "nullity")
        do_homogeneity = any(tests_tokens :== "homogeneity")
        do_gini = any(tests_tokens :== "gini")
    }

    tau_boot = J(B, nq, 0)

    if (do_gini) {
        nu_plus_store = J(B, nq, 0)
        nu_minus_store = J(B, nq, 0)
    }

    for (b = 1; b <= B; b++) {
        xi = (use_xi_mat ? xi_mat_data[,b] : rnormal(n, 1, 0, 1))
        plus_sums = xi' * e1_w_plus
        minus_sums = xi' * e1_w_minus
        out_sharp = plus_sums - minus_sums

        if (do_gini) {
            nu_plus_store[b,] = plus_sums
            nu_minus_store[b,] = minus_sums
        }

        if (!is_fuzzy) {
            tau_boot[b,] = out_sharp
        }
        else {
            denom_term = (xi' * e2_w_plus[,1]) - (xi' * e2_w_minus[,1])
            if (abs(denominator) > 1e-8) {
                top = denominator * out_sharp - num_diff :* denom_term
                tau_boot[b,] = top / (denominator^2)
            }
            else {
                tau_boot[b,] = J(1, nq, 0)
            }
        }
    }

    // Uniform confidence bands
    supvals = J(B, 1, 0)
    for (b = 1; b <= B; b++) {
        supvals[b] = max(abs(tau_boot[b,]))
    }

    sup_sorted = sort(supvals, 1)

    prob = level / 100
    alpha_level = 1 - prob
    index_c = ceil(prob * (B + 1))
    if (index_c < 1) index_c = 1
    if (index_c > B) index_c = B

    cval = sup_sorted[index_c]

    cb_lower = tau_orig :- cval
    cb_upper = tau_orig :+ cval

    // Prepare quantile grid and test ranges
    q_grid = strtoreal(tokens(quantiles_str))

    range_tokens = strtoreal(tokens(test_ranges_str))
    if (cols(range_tokens) < 2) {
        range_pairs = (min(q_grid), max(q_grid))
    }
    else {
        n_pairs = floor(cols(range_tokens) / 2)
        range_pairs = J(n_pairs, 2, .)
        idx = 1
        for (r = 1; r <= n_pairs; r++) {
            range_pairs[r,1] = range_tokens[idx]
            range_pairs[r,2] = range_tokens[idx+1]
            idx = idx + 2
        }
    }

    // For simplicity, use the first range if multiple provided
    range_lo = range_pairs[1,1]
    range_hi = range_pairs[1,2]

    range_idx = J(0, 1, .)
    for (j = 1; j <= nq; j++) {
        if (q_grid[j] >= range_lo & q_grid[j] <= range_hi) {
            range_idx = range_idx \ j
        }
    }

    if (rows(range_idx) == 0) {
        range_idx = (1..nq)'
    }

    pvals = J(1, 3, .)

    if (do_nullity) {
        test_stat_null = max(abs(tau_orig[1, range_idx']))

        supvals_null = J(B, 1, 0)
        for (b = 1; b <= B; b++) {
            supvals_null[b] = max(abs(tau_boot[b, range_idx']))
        }

        sup_null_sorted = sort(supvals_null, 1)
        idx_null = ceil(prob * (B + 1))
        if (idx_null < 1) idx_null = 1
        if (idx_null > B) idx_null = B
        crit_null = sup_null_sorted[idx_null]

        p_null = mean(supvals_null :>= test_stat_null)
        pvals[1,1] = p_null
    }

    if (do_homogeneity) {
        range_tau = tau_orig[1, range_idx']
        mbar = mean(range_tau')
        test_stat_homo = max(abs(range_tau :- mbar))

        supvals_homo = J(B, 1, 0)
        for (b = 1; b <= B; b++) {
            draw = tau_boot[b, range_idx']
            draw_mean = mean(draw')
            supvals_homo[b] = max(abs(draw :- draw_mean))
        }

        sup_homo_sorted = sort(supvals_homo, 1)
        idx_homo = ceil(prob * (B + 1))
        if (idx_homo < 1) idx_homo = 1
        if (idx_homo > B) idx_homo = B

        crit_homo = sup_homo_sorted[idx_homo]

        p_homo = mean(supvals_homo :>= test_stat_homo)
        pvals[1,2] = p_homo
    }

    if (do_gini) {
        gini_above = r3d_gini_from_quantile(q_grid, int_plus)
        gini_below = r3d_gini_from_quantile(q_grid, int_minus)
        gini_diff = gini_above - gini_below

        gini_boot = J(B, 1, 0)
        for (b = 1; b <= B; b++) {
            boot_plus = int_plus + nu_plus_store[b,]
            boot_minus = int_minus + nu_minus_store[b,]
            gini_boot[b] = r3d_gini_from_quantile(q_grid, boot_plus) - ///
                           r3d_gini_from_quantile(q_grid, boot_minus)
        }

        test_stat_gini = abs(gini_diff)

        gini_sorted = sort(abs(gini_boot), 1)
        idx_gini = ceil(prob * (B + 1))
        if (idx_gini < 1) idx_gini = 1
        if (idx_gini > B) idx_gini = B

        crit_gini = gini_sorted[idx_gini]

        p_gini = mean(abs(gini_boot) :>= test_stat_gini)
        pvals[1,3] = p_gini
    }

    st_matrix(bs_results_name, tau_boot)
    st_matrix(pvals_name, pvals)
    st_matrix("r(cb_lower)", cb_lower)
    st_matrix("r(cb_upper)", cb_upper)
}

// ============================================================================
// GINI COEFFICIENT
// ============================================================================

// Calculate Gini coefficient from quantile function
real scalar r3d_gini_from_quantile(real rowvector quantiles, real rowvector qfunc)
{
    real scalar nq, i, mean_val, area, gini, cumulative
    real rowvector lorenz, q, qf

    nq = cols(quantiles)
    if (nq < 2) return(.)

    q = quantiles
    qf = qfunc

    // Add endpoint at 0 if missing (Lorenz curve starts at origin)
    if (q[1] > 0) {
        q  = (0, q)
        qf = (0, qf)
    }
    // Add endpoint at 1 if missing
    if (q[cols(q)] < 1) {
        q  = (q, 1)
        qf = (qf, qf[cols(qf)])
    }
    nq = cols(q)

    // Approximate mean from quantile function
    mean_val = 0
    for (i = 1; i < nq; i++) {
        mean_val = mean_val + 0.5 * (qf[i+1] + qf[i]) * (q[i+1] - q[i])
    }

    if (mean_val <= 0) return(0)

    // Calculate Lorenz curve
    lorenz = J(1, nq, 0)
    cumulative = 0

    for (i = 1; i < nq; i++) {
        cumulative = cumulative + 0.5 * (qf[i+1] + qf[i]) * (q[i+1] - q[i])
        lorenz[i+1] = cumulative / mean_val
    }

    // Area under Lorenz curve
    area = 0
    for (i = 1; i < nq; i++) {
        area = area + 0.5 * (lorenz[i] + lorenz[i+1]) * (q[i+1] - q[i])
    }

    // Gini coefficient
    gini = 1 - 2 * area

    return(max((0, min((1, gini)))))
}

// Test for Gini coefficient difference
void r3d_test_gini(string scalar xvar, string scalar Qmat,
                   real scalar cutoff, string scalar touse,
                   string scalar pval_name)
{
    real matrix X, Q
    real scalar n, nq, gini_plus, gini_minus, diff, i, j
    real rowvector quantiles, qfunc_plus, qfunc_minus
    real colvector idx_plus, idx_minus

    // Get data
    st_view(X, ., xvar, touse)
    Q = st_matrix(Qmat)
    n = rows(X)
    nq = cols(Q)

    // Create quantile grid
    quantiles = J(1, nq, 0)
    for (i = 1; i <= nq; i++) {
        quantiles[i] = i / (nq + 1)
    }

    // Compute average quantile functions above and below cutoff
    idx_plus = (X :>= 0)
    idx_minus = (X :< 0)

    qfunc_plus = J(1, nq, 0)
    qfunc_minus = J(1, nq, 0)

    for (j = 1; j <= nq; j++) {
        qfunc_plus[j] = mean(select(Q[,j], idx_plus))
        qfunc_minus[j] = mean(select(Q[,j], idx_minus))
    }

    // Calculate Gini coefficients
    gini_plus = r3d_gini_from_quantile(quantiles, qfunc_plus)
    gini_minus = r3d_gini_from_quantile(quantiles, qfunc_minus)

    // Difference
    diff = gini_plus - gini_minus

    // Store results (p-value would need bootstrap for proper inference)
    st_numscalar(pval_name, 0.5)  // Placeholder
    st_numscalar("r(gini_plus)", gini_plus)
    st_numscalar("r(gini_minus)", gini_minus)
    st_numscalar("r(gini_diff)", diff)
}

// ============================================================================
// CSV LOADER (for xi_mat)
// ============================================================================

// Load a CSV file (with header row) into a numeric matrix
real matrix _r3d_load_csv(string scalar filepath)
{
    real scalar fh, nc
    string scalar line
    string rowvector tokens_str
    real matrix data
    real rowvector row_data

    fh = fopen(filepath, "r")

    // Skip header
    line = fget(fh)

    // Read all lines
    nc = 0
    data = J(0, 0, .)

    while ((line = fget(fh)) != J(0, 0, "")) {
        tokens_str = tokens(subinstr(line, ",", " ", .))
        if (nc == 0) {
            nc = cols(tokens_str)
            data = J(0, nc, .)
        }
        row_data = strtoreal(tokens_str)
        data = data \ row_data
    }

    fclose(fh)
    return(data)
}

// ============================================================================
// TEST HELPER: locpoly wrapper for test_equivalence.do
// ============================================================================

// Runs r3d_locpoly for both sides and stores results in Stata matrices
// __AP (alpha_plus), __AM (alpha_minus), __WP (w_plus), __WM (w_minus)
void _r3d_test_locpoly(string scalar xvar, string scalar Qmatname,
                       string scalar h_matname, real scalar p, real scalar ktype)
{
    real colvector X
    real matrix Y, h_mat
    real matrix alpha_plus, alpha_minus, w_plus, w_minus

    X = st_data(., xvar)
    Y = st_matrix(Qmatname)
    h_mat = st_matrix(h_matname)

    // Plus side (side=1)
    r3d_locpoly(X, Y, h_mat, p, ktype, 1, alpha_plus, w_plus)
    // Minus side (side=0)
    r3d_locpoly(X, Y, h_mat, p, ktype, 0, alpha_minus, w_minus)

    st_matrix("__AP", alpha_plus)
    st_matrix("__AM", alpha_minus)
    st_matrix("__WP", w_plus)
    st_matrix("__WM", w_minus)
}

end
