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
    real scalar n, nq, i, j
    string scalar cmd

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

    // Create temporary variables in Stata for plugin interface
    stata("quietly drop _plugin_*", 1)  // Allow failure if they don't exist

    // Create X variable
    stata("quietly generate double _plugin_x = .")
    st_store(., "_plugin_x", X)

    // Create Y matrix variables
    for (j = 1; j <= nq; j++) {
        stata("quietly generate double _plugin_y" + strofreal(j) + " = .")
        st_store(., "_plugin_y" + strofreal(j), Y[,j])
    }

    // Create output variables for alpha coefficients
    for (i = 1; i <= p+1; i++) {
        for (j = 1; j <= nq; j++) {
            stata("quietly generate double _plugin_a" + strofreal(i) + "_" + strofreal(j) + " = 0")
        }
    }

    // Create output variables for weights
    for (j = 1; j <= nq; j++) {
        stata("quietly generate double _plugin_w" + strofreal(j) + " = 0")
    }

    // Build variable list for plugin call
    string scalar varlist
    varlist = "_plugin_x"
    for (j = 1; j <= nq; j++) {
        varlist = varlist + " _plugin_y" + strofreal(j)
    }
    for (i = 1; i <= p+1; i++) {
        for (j = 1; j <= nq; j++) {
            varlist = varlist + " _plugin_a" + strofreal(i) + "_" + strofreal(j)
        }
    }
    for (j = 1; j <= nq; j++) {
        varlist = varlist + " _plugin_w" + strofreal(j)
    }

    string scalar plugin_name
    plugin_name = st_global("R3D_PLUGIN_NAME")
    if (plugin_name == "") plugin_name = "r3d_plugin"

    // Call plugin directly
    cmd = "plugin call " + plugin_name + " " + varlist + ", " +
          strofreal(p) + " " + strofreal(side) + " " +
          strofreal(kernel_type) + " " + strofreal(nq)

    // Add bandwidth values
    for (j = 1; j <= nq; j++) {
        cmd = cmd + " " + strofreal(h_mat[j,1])
    }

    real scalar rc
    rc = _stata("quietly " + cmd)

    if (rc != 0) {
        // Drop temporary variables before exiting
        stata("quietly drop _plugin_*")
        return(rc)
    }

    // Read results back into matrices
    alpha = J(p+1, nq, 0)
    weights = J(n, nq, 0)

    for (i = 1; i <= p+1; i++) {
        for (j = 1; j <= nq; j++) {
            // Get first value (coefficients are constant)
            real scalar val
            val = st_data(1, "_plugin_a" + strofreal(i) + "_" + strofreal(j))
            alpha[i,j] = val
        }
    }

    for (j = 1; j <= nq; j++) {
        weights[,j] = st_data(., "_plugin_w" + strofreal(j))
    }

    // Clean up temporary variables
    stata("quietly drop _plugin_*")

    return(0)
}

end
