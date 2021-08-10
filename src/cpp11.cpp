// Generated by cpp11: do not edit by hand
// clang-format off


#include "cpp11/declarations.hpp"

// misc.cpp
cpp11::doubles_matrix collapseData(cpp11::list mat_list, cpp11::strings colnames);
extern "C" SEXP _openCyto_collapseData(SEXP mat_list, SEXP colnames) {
  BEGIN_CPP11
    return cpp11::as_sexp(collapseData(cpp11::as_cpp<cpp11::decay_t<cpp11::list>>(mat_list), cpp11::as_cpp<cpp11::decay_t<cpp11::strings>>(colnames)));
  END_CPP11
}
// solve_LSAP.cpp
std::vector<int> solve_LSAP_cpp(cpp11::doubles_matrix mat);
extern "C" SEXP _openCyto_solve_LSAP_cpp(SEXP mat) {
  BEGIN_CPP11
    return cpp11::as_sexp(solve_LSAP_cpp(cpp11::as_cpp<cpp11::decay_t<cpp11::doubles_matrix>>(mat)));
  END_CPP11
}
// unlockNamespace.cpp
cpp11::logicals unlockNamespace(cpp11::sexp env);
extern "C" SEXP _openCyto_unlockNamespace(SEXP env) {
  BEGIN_CPP11
    return cpp11::as_sexp(unlockNamespace(cpp11::as_cpp<cpp11::decay_t<cpp11::sexp>>(env)));
  END_CPP11
}

extern "C" {
/* .Call calls */
extern SEXP _openCyto_collapseData(SEXP, SEXP);
extern SEXP _openCyto_solve_LSAP_cpp(SEXP);
extern SEXP _openCyto_unlockNamespace(SEXP);

static const R_CallMethodDef CallEntries[] = {
    {"_openCyto_collapseData",    (DL_FUNC) &_openCyto_collapseData,    2},
    {"_openCyto_solve_LSAP_cpp",  (DL_FUNC) &_openCyto_solve_LSAP_cpp,  1},
    {"_openCyto_unlockNamespace", (DL_FUNC) &_openCyto_unlockNamespace, 1},
    {NULL, NULL, 0}
};
}

extern "C" void R_init_openCyto(DllInfo* dll){
  R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
  R_forceSymbols(dll, TRUE);
}
