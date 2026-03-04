/*
 * r3d_plugin.c - STATA plugin interface for R3D Fortran routines
 * This provides the C interface that STATA expects for plugins
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "stplugin.h"

// Global Stata plugin pointer (required for plugin interface)
ST_plugin *_stata_;

// Fortran routine declaration (with underscore for Fortran name mangling)
extern void stata_locweights_(double *x, double *ymat, int *n, int *p,
                             double *h, int *side, int *kernel_type,
                             double *alpha, double *wint, int *info, int *nq);

// Plugin initialization (required)
STDLL pginit(ST_plugin *p)
{
    _stata_ = p;
    return(SD_PLUGINVER);
}

// Main plugin entry point
STDLL stata_call(int argc, char *argv[])
{
    ST_retcode rc;
    int n, p, side, kernel_type, nq, info;
    double *x, *ymat, *h, *alpha, *wint;
    int i, j;
    
    // Check minimum arguments
    if (argc < 9) {
        SF_error("r3d_plugin: insufficient arguments\n");
        SF_error("Usage: plugin call r3d_plugin xvar p side kernel nq h_1...h_nq\n");
        return 198;
    }
    
    // Parse arguments
    // argv[0]: plugin name
    // argv[1]: X variable name
    // argv[2]: polynomial order
    // argv[3]: side (0=minus, 1=plus)
    // argv[4]: kernel type (1=triangular, 2=epanechnikov, 3=uniform)
    // argv[5]: number of quantiles
    // argv[6..6+nq-1]: bandwidth values
    // argv[6+nq]: Y variable prefix
    
    // Get sample size
    if ((rc = SF_in1()) != 0) return rc;
    n = SF_nobs();
    if ((rc = SF_in2()) != 0) return rc;
    
    // Get parameters
    p = atoi(argv[2]);
    if (p < 0 || p > 9) {
        SF_error("r3d_plugin: polynomial order must be between 0 and 9\n");
        return 198;
    }
    
    side = atoi(argv[3]);
    kernel_type = atoi(argv[4]);
    nq = atoi(argv[5]);
    
    if (nq <= 0 || nq > 100) {
        SF_error("r3d_plugin: number of quantiles must be between 1 and 100\n");
        return 198;
    }
    
    // Check we have enough arguments for bandwidths
    if (argc < 6 + nq + 1) {
        SF_error("r3d_plugin: insufficient bandwidth arguments\n");
        return 198;
    }
    
    // Allocate memory
    x = (double *) malloc(n * sizeof(double));
    ymat = (double *) malloc(n * nq * sizeof(double));
    h = (double *) malloc(nq * sizeof(double));
    alpha = (double *) malloc((p+1) * nq * sizeof(double));
    wint = (double *) malloc(n * nq * sizeof(double));
    
    if (!x || !ymat || !h || !alpha || !wint) {
        SF_error("r3d_plugin: memory allocation failed\n");
        if (x) free(x);
        if (ymat) free(ymat);
        if (h) free(h);
        if (alpha) free(alpha);
        if (wint) free(wint);
        return 198;
    }
    
    // Read X variable (first variable in varlist)
    for (i = SF_in1(); i <= SF_in2(); i++) {
        if (SF_ifobs(i)) {
            if ((rc = SF_vdata(1, i, &x[i-SF_in1()])) != 0) {
                SF_error("r3d_plugin: cannot read X variable\n");
                goto cleanup;
            }
        }
    }
    
    // Read bandwidth values
    for (j = 0; j < nq; j++) {
        h[j] = atof(argv[6 + j]);
        if (h[j] <= 0) {
            SF_error("r3d_plugin: bandwidth must be positive\n");
            rc = 198;
            goto cleanup;
        }
    }
    
    // Read Y matrix (variables 2 through nq+1)
    for (j = 0; j < nq; j++) {
        for (i = SF_in1(); i <= SF_in2(); i++) {
            if (SF_ifobs(i)) {
                if ((rc = SF_vdata(j+2, i, &ymat[(i-SF_in1()) + j*n])) != 0) {
                    SF_error("r3d_plugin: cannot read Y variable\n");
                    goto cleanup;
                }
            }
        }
    }
    
    // Call Fortran routine
    stata_locweights_(x, ymat, &n, &p, h, &side, &kernel_type, 
                     alpha, wint, &info, &nq);
    
    if (info != 0) {
        SF_error("r3d_plugin: locweights computation failed\n");
        rc = 198;
        goto cleanup;
    }
    
    // Store results back to Stata variables
    // Store alpha coefficients in variables nq+2 through nq+2+(p+1)*nq-1
    int var_idx = nq + 2;
    for (j = 0; j < nq; j++) {
        for (i = 0; i < p+1; i++) {
            for (int k = SF_in1(); k <= SF_in2(); k++) {
                if (SF_ifobs(k)) {
                    double val = (k == SF_in1()) ? alpha[i + j*(p+1)] : 0.0;
                    if ((rc = SF_vstore(var_idx, k, val)) != 0) {
                        SF_error("r3d_plugin: cannot store alpha coefficient\n");
                        goto cleanup;
                    }
                }
            }
            var_idx++;
        }
    }
    
    // Store weight matrix in remaining variables
    for (j = 0; j < nq; j++) {
        for (i = SF_in1(); i <= SF_in2(); i++) {
            if (SF_ifobs(i)) {
                if ((rc = SF_vstore(var_idx + j, i, wint[(i-SF_in1()) + j*n])) != 0) {
                    SF_error("r3d_plugin: cannot store weight\n");
                    goto cleanup;
                }
            }
        }
    }
    
    rc = 0;
    
cleanup:
    free(x);
    free(ymat);
    free(h);
    free(alpha);
    free(wint);
    
    return rc;
}