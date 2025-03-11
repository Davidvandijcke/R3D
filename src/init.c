#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>

/* FORTRAN routines */
void F77_NAME(locweights)(double *X, double *YMAT, int *N, int *P, 
              double *H, int *SIDE, double *KERNELW,
              double *ALPHA, double *WINT, int *INFO, int *NQ);

static const R_FortranMethodDef FortranEntries[] = {
  {"locweights", (DL_FUNC) &F77_NAME(locweights), 11},
  {NULL, NULL, 0}
};

void R_init_YourPackageName(DllInfo *dll) {
  R_registerRoutines(dll, NULL, NULL, FortranEntries, NULL);
  R_useDynamicSymbols(dll, FALSE);
}
