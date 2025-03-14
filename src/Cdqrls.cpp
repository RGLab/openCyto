/*
 * from stats package
 * src/library/stats/src/lm.c
 */
#include <cpp11.hpp>
#include <R.h>
#include <Rinternals.h>
#include <R_ext/Applic.h>



#ifdef ENABLE_NLS
#include <libintl.h>
#define _(String) dgettext ("stats", String)
#else
#define _(String) (String)
#endif

[[cpp11::register]]
SEXP Cdqrls(SEXP x, SEXP y, SEXP tol)
{
  SEXP ans;
  SEXP qr, coefficients, residuals, effects, pivot, qraux;
  int n, ny = 0, p, rank, nprotect = 4, pivoted = 0;
  double rtol = Rf_asReal(tol), *work;

  ans = Rf_getAttrib(x, R_DimSymbol);
  int *dims = INTEGER(ans);
  n = dims[0]; p = dims[1];
  if(n) ny = (int)(XLENGTH(y)/n); /* y :  n x ny, or an n - vector */

  /* These lose attributes, so do after we have extracted dims */
  if (TYPEOF(x) != REALSXP) {
    PROTECT(x = Rf_coerceVector(x, REALSXP));
    nprotect++;
  }
  if (TYPEOF(y) != REALSXP) {
    PROTECT(y = Rf_coerceVector(y, REALSXP));
    nprotect++;
  }
  
  double *rptr = REAL(x);
  for (R_xlen_t i = 0 ; i < XLENGTH(x) ; i++)
    if(!R_FINITE(rptr[i])) Rf_error(_("NA/NaN/Inf in '%s'"), "x");
    
    rptr = REAL(y);
    for (R_xlen_t i = 0 ; i < XLENGTH(y) ; i++)
      if(!R_FINITE(rptr[i])) Rf_error(_("NA/NaN/Inf in '%s'"), "y");
      
      const char *ansNms[] = {"qr", "coefficients", "residuals", "effects",
                              "rank", "pivot", "qraux", "tol", "pivoted", ""};
      PROTECT(ans = Rf_mkNamed(VECSXP, ansNms));
      SET_VECTOR_ELT(ans, 0, qr = Rf_shallow_duplicate(x));
      coefficients = (ny > 1) ? Rf_allocMatrix(REALSXP, p, ny) : Rf_allocVector(REALSXP, p);
      PROTECT(coefficients);
      SET_VECTOR_ELT(ans, 1, coefficients);
      SET_VECTOR_ELT(ans, 2, residuals = Rf_shallow_duplicate(y));
      SET_VECTOR_ELT(ans, 3, effects = Rf_shallow_duplicate(y));
      PROTECT(pivot = Rf_allocVector(INTSXP, p));
      int *ip = INTEGER(pivot);
      for(int i = 0; i < p; i++) ip[i] = i+1;
      SET_VECTOR_ELT(ans, 5, pivot);
      PROTECT(qraux = Rf_allocVector(REALSXP, p));
      SET_VECTOR_ELT(ans, 6, qraux);
      SET_VECTOR_ELT(ans, 7, tol);
      
      work = (double *) R_alloc(2 * p, sizeof(double));
      F77_CALL(dqrls)(REAL(qr), &n, &p, REAL(y), &ny, &rtol,
               REAL(coefficients), REAL(residuals), REAL(effects),
               &rank, INTEGER(pivot), REAL(qraux), work);
      SET_VECTOR_ELT(ans, 4, Rf_ScalarInteger(rank));
      for(int i = 0; i < p; i++)
        if(ip[i] != i+1) { pivoted = 1; break; }
        SET_VECTOR_ELT(ans, 8, Rf_ScalarLogical(pivoted));
        UNPROTECT(nprotect);
        
        return ans;
}