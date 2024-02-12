#include <R.h>
#include <R_ext/RS.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* .Fortran calls */
extern void qfc (double *, double *, int *, int *, double *, double *, int *, double *, double * , int *, double *);
extern void F77_NAME(gausq2)(void *, void *, void *, void *, void *);
extern void F77_NAME(mybnorm)(void *, void *, void *, void *, void *);
extern void F77_NAME(mytnorm)(void *, void *, void *, void *, void *, void *);
extern void F77_NAME(myunorm)(void *, void *, void *);


static const R_FortranMethodDef FortranEntries[] = {
  {"gausq2", (DL_FUNC) &F77_NAME(gausq2), 5},
  {"mybnorm", (DL_FUNC) &F77_NAME(mybnorm), 5},
  {"mytnorm", (DL_FUNC) &F77_NAME(mytnorm), 6},
  {"myunorm", (DL_FUNC) &F77_NAME(myunorm), 3},
  {NULL, NULL, 0}
};

static const R_CMethodDef CEntries[] = {
  {"qfc", (DL_FUNC) &qfc, 11},
  {NULL, NULL, 0}
};

void R_init_gc(DllInfo *dll)
{
    R_registerRoutines(dll, CEntries, NULL, FortranEntries, NULL);
    R_useDynamicSymbols(dll, FALSE);
}


