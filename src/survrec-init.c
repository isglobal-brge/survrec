#include <R_ext/RS.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
Check these declarations against the C/Fortran source code.
*/

/* .Fortran calls */
extern void F77_NAME(bootmedian)(int *, int *, double *, double *, int *,
                     double *, double *, int *, 
                     int *, double *, double *);
extern void F77_NAME(distinctfailed)(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(emalgo)(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(mlevalue)(void *, void *, void *, void *);
extern void F77_NAME(searchforseed)(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(wc2)(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(wcPLE)(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(estalpha)(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);

static const R_FortranMethodDef FortranEntries[] = {
  {"bootmedian",     (DL_FUNC) &F77_NAME(bootmedian),     11},
  {"distinctfailed", (DL_FUNC) &F77_NAME(distinctfailed), 10},
  {"emalgo",         (DL_FUNC) &F77_NAME(emalgo),         11},
  {"mlevalue",       (DL_FUNC) &F77_NAME(mlevalue),        4},
  {"searchforseed",  (DL_FUNC) &F77_NAME(searchforseed),  12},
  {"wc2",            (DL_FUNC) &F77_NAME(wc2),            17},
  {"wcPLE",          (DL_FUNC) &F77_NAME(wc2),            11},
  {"estalpha",       (DL_FUNC) &F77_NAME(wc2),            11},
  {NULL, NULL, 0}
};

void R_init_survrec(DllInfo *dll)
{
  R_registerRoutines(dll, NULL, NULL, FortranEntries, NULL);
  R_useDynamicSymbols(dll, FALSE);
}