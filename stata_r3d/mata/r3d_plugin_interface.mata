/*
 * r3d_plugin_interface.mata - Interface to Fortran plugin for fast computation
 */

version 18.0
mata:

// Check if plugin is available
real scalar r3d_plugin_available()
{
    if (st_global("R3D_FORCE_PLUGIN") == "1") {
        return(1)
    }
    if (st_global("R3D_USE_PLUGIN") == "1") {
        return(1)
    }
    return(0)
}

// Call Fortran plugin for local polynomial regression
real scalar r3d_plugin_locweights(real colvector X, real matrix Y, real matrix h_mat,
                                 real scalar p, real scalar kernel_type, real scalar side,
                                 real matrix alpha, real matrix weights)
{
    // F8: NOT thread-safe: writes and drops __r3d_pi_* variables from the live Stata dataset. Do not call concurrently.

    real scalar n, nq, i, j
    string scalar cmd

    // F3: Validate that X is a column vector and neither X nor Y contains missing values
    assert(cols(X)==1)
    if (hasmissing(X) | hasmissing(Y)) {
        errprintf("r3d_plugin_call: X and Y must not contain missing values\n")
        return(198)
    }

    n = rows(X)
    nq = cols(Y)

    // Normalise bandwidth input to a column vector of length nq
    if (rows(h_mat) == 1 & cols(h_mat) == 1) {
        h_mat = J(nq, 1, h_mat[1,1])
    }
    else if (rows(h_mat) == 1 & cols(h_mat) == nq) {
        h_mat = h_mat'
    }
    else if (!(rows(h_mat) == nq & cols(h_mat) == 1)) {
        errprintf("r3d_plugin_locweights(): bandwidth vector has incorrect dimensions\n")
        return(198)
    }

    // F4: Validate bandwidth positivity
    if (min(h_mat) <= 0) {
        errprintf("r3d_plugin_call: bandwidth must be positive\n")
        return(198)
    }

    // F6: Check available variable slots before creating temporaries
    real scalar _need_vars
    _need_vars = 1 + nq + nq*(p+1) + nq
    if (c("k") + _need_vars > c("maxvar")) {
        errprintf("r3d_plugin_call: insufficient variable slots (need %g, have %g free)\n", _need_vars, c("maxvar")-c("k"))
        return(198)
    }

    // F5: Use unique prefix __r3d_pi_* to avoid colliding with user variables;
    //     use capture so the drop does not error when variables are absent
    stata("quietly capture drop __r3d_pi_*", 1)

    // Create X variable
    stata("generate double __r3d_pi_x = .")
    st_store(., "__r3d_pi_x", X)

    // Create Y matrix variables
    for (j = 1; j <= nq; j++) {
        stata("generate double __r3d_pi_y" + strofreal(j) + " = .")
        st_store(., "__r3d_pi_y" + strofreal(j), Y[,j])
    }

    // Create output variables for alpha coefficients
    for (i = 1; i <= p+1; i++) {
        for (j = 1; j <= nq; j++) {
            stata("generate double __r3d_pi_a" + strofreal(i) + "_" + strofreal(j) + " = 0")
        }
    }

    // Create output variables for weights
    for (j = 1; j <= nq; j++) {
        stata("generate double __r3d_pi_w" + strofreal(j) + " = 0")
    }

    // Build variable list for plugin call
    string scalar varlist
    varlist = "__r3d_pi_x"
    for (j = 1; j <= nq; j++) {
        varlist = varlist + " __r3d_pi_y" + strofreal(j)
    }
    for (i = 1; i <= p+1; i++) {
        for (j = 1; j <= nq; j++) {
            varlist = varlist + " __r3d_pi_a" + strofreal(i) + "_" + strofreal(j)
        }
    }
    for (j = 1; j <= nq; j++) {
        varlist = varlist + " __r3d_pi_w" + strofreal(j)
    }

    string scalar plugin_name
    plugin_name = st_global("R3D_PLUGIN_NAME")
    if (plugin_name == "") plugin_name = "r3d_plugin"

    // F1: Use integer-format serialization to prevent scientific notation (e.g. '1e+02')
    //     which the C plugin arg parser would reject or misparse
    assert(p < 100)
    cmd = "plugin call " + plugin_name + " " + varlist + ", " +
          strofreal(p, "%12.0f") + " " + strofreal(side, "%12.0f") + " " +
          strofreal(kernel_type, "%12.0f") + " " + strofreal(nq)

    // Add bandwidth values
    for (j = 1; j <= nq; j++) {
        cmd = cmd + " " + strofreal(h_mat[j,1])
    }

    real scalar rc
    // F7: Do not wrap plugin call in quietly — plugin errors must be visible
    rc = _stata(cmd)

    if (rc != 0) {
        errprintf("r3d_plugin_call: plugin returned rc=%g; cmd was: %s\n", rc, cmd)
        // Drop temporary variables before exiting
        stata("quietly capture drop __r3d_pi_*")
        return(rc)
    }

    // Read results back into matrices
    alpha = J(p+1, nq, 0)
    weights = J(n, nq, 0)

    // F2: The plugin writes all coefficient values into row 1 only (broadcast convention).
    //     If the plugin instead writes per-row, change st_data(1,...) to st_data(j,...).
    //     Sanity check: confirm the first alpha variable has the expected number of rows.
    assert(rows(st_data(., "__r3d_pi_a1_1"))==n)
    for (i = 1; i <= p+1; i++) {
        for (j = 1; j <= nq; j++) {
            // Get first value (coefficients are constant)
            real scalar val
            val = st_data(1, "__r3d_pi_a" + strofreal(i) + "_" + strofreal(j))
            alpha[i,j] = val
        }
    }

    for (j = 1; j <= nq; j++) {
        weights[,j] = st_data(., "__r3d_pi_w" + strofreal(j))
    }

    // Clean up temporary variables
    stata("quietly capture drop __r3d_pi_*")

    return(0)
}

end
