#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

extern SEXP CircularDDM_rvm(SEXP, SEXP, SEXP);
extern SEXP CircularDDM_besselzero(SEXP, SEXP, SEXP);
extern SEXP CircularDDM_logLik_resp(SEXP, SEXP);
extern SEXP CircularDDM_logLik_dt(SEXP, SEXP, SEXP);
extern SEXP CircularDDM_dcddm(SEXP, SEXP, SEXP);
extern SEXP CircularDDM_rcddm1(SEXP, SEXP, SEXP);
extern SEXP CircularDDM_rcddm2_internal(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
  {"CircularDDM_rvm",         (DL_FUNC) &CircularDDM_rvm, 3},
  {"CircularDDM_besselzero",  (DL_FUNC) &CircularDDM_besselzero, 3},
  {"CircularDDM_logLik_resp", (DL_FUNC) &CircularDDM_logLik_resp, 2},
  {"CircularDDM_logLik_dt",   (DL_FUNC) &CircularDDM_logLik_dt, 3},
  {"CircularDDM_dcddm",   (DL_FUNC) &CircularDDM_dcddm, 3},
  {"CircularDDM_rcddm1",  (DL_FUNC) &CircularDDM_rcddm1, 3},
  {"CircularDDM_rcddm2_internal",  (DL_FUNC) &CircularDDM_rcddm2_internal, 6},
  {NULL, NULL, 0}
};

void R_init_CircularDDM(DllInfo *dll) {
  R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
}

