#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
Check these declarations against the C/Fortran source code.
*/

/* .C calls */
extern void pstable_c(void *, void *, void *, void *, void *, void *, void *);
extern void stable(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);

static const R_CMethodDef CEntries[] = {
  {"pstable_c", (DL_FUNC) &pstable_c,  7},
  {"stable",    (DL_FUNC) &stable,    10},
  {NULL, NULL, 0}
};

void R_init_stable(DllInfo *dll)
{
  R_registerRoutines(dll, CEntries, NULL, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
}