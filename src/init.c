#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
 Check these declarations against the C/Fortran source code.
 */

/* .C calls */
extern void cov2_R(void *, void *, void *);
extern void entropy_R(void *, void *, void *);
extern void ess_R(void *, void *, void *);
extern void inverse_cdf_weights_R(void *, void *, void *, void *, void *);
extern void multinomial_resample_R(void *, void *, void *, void *);
extern void renormalize_R(void *, void *, void *);
extern void rep2id_R(void *, void *, void *);
extern void residual_resample_R(void *, void *, void *, void *, void *);
extern void stratified_resample_R(void *, void *, void *, void *);
extern void systematic_resample_R(void *, void *, void *, void *);

static const R_CMethodDef CEntries[] = {
  {"cov2_R",                 (DL_FUNC) &cov2_R,                 3},
  {"entropy_R",              (DL_FUNC) &entropy_R,              3},
  {"ess_R",                  (DL_FUNC) &ess_R,                  3},
  {"inverse_cdf_weights_R",  (DL_FUNC) &inverse_cdf_weights_R,  5},
  {"multinomial_resample_R", (DL_FUNC) &multinomial_resample_R, 4},
  {"renormalize_R",          (DL_FUNC) &renormalize_R,          3},
  {"rep2id_R",               (DL_FUNC) &rep2id_R,               3},
  {"residual_resample_R",    (DL_FUNC) &residual_resample_R,    5},
  {"stratified_resample_R",  (DL_FUNC) &stratified_resample_R,  4},
  {"systematic_resample_R",  (DL_FUNC) &systematic_resample_R,  4},
  {NULL, NULL, 0}
};

void R_init_smcUtils(DllInfo *dll)
{
  R_registerRoutines(dll, CEntries, NULL, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
}