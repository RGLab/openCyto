#include <cpp11.hpp>
#include <iostream>
#include <vector>
using namespace std;
SEXP Cdqrls(SEXP x, SEXP y, SEXP tol);

cpp11::list lm_wfit_cpp(cpp11::doubles_matrix<> x, cpp11::doubles y,
                        cpp11::doubles w) {
  using namespace cpp11::literals;  // so we can use ""_nm syntax

  cpp11::writable::list fit_res;
  // weights with zero-values are not possible based on current singletGate use
  // case due to its all-1 default values of x[,1] (which is used as the initial
  // weight input) for fast_rlm
  int n = w.size();
  for (int i = 0; i < w.size(); i++) {
    if (w[i] == 0)
      cpp11::stop("weights with zero-values are currently not supported");
  }
  int nx = x.ncol();

  //
  cpp11::writable::doubles wts(n), y1(n);
  cpp11::writable::doubles_matrix<> x1(n, nx);
  for (int i = 0; i < n; i++) {
    wts[i] = sqrt(w[i]);
    y1[i] = y[i] * wts[i];
    for (int j = 0; j < nx; j++) x1(i, j) = x(i, j) * wts[i];
  }
  cpp11::writable::doubles tol(1);
  tol[0] = 1.0e-7;

  auto z = cpp11::as_cpp<cpp11::writable::list>(Cdqrls(x1, y1, tol));

  auto coef = cpp11::as_cpp<cpp11::writable::doubles>(z["coefficients"]);
  auto pivot = cpp11::as_cpp<cpp11::integers>(z["pivot"]);
  int rank = cpp11::as_cpp<cpp11::integers>(z["rank"])[0];

  if (rank < nx) {
    for (int i = rank; i < nx; i++) coef[i] = NA_REAL;
  }
  cpp11::writable::doubles coef1 = coef;
  for (int i = 0; i < pivot.size(); i++) {
    int idx = pivot[i] - 1;
    double val = coef1[idx];
    coef[i] = val;
  }

  auto residuals = cpp11::as_cpp<cpp11::writable::doubles>(z["residuals"]);
  cpp11::writable::doubles fitted_values(n);
  for (int i = 0; i < n; i++) {
    residuals[i] /= wts[i];
    fitted_values[i] = y[i] - residuals[i];
  }

  fit_res.push_back("fitted.values"_nm = fitted_values);
  fit_res.push_back("weights"_nm = w);
  fit_res.push_back("residuals"_nm = residuals);
  fit_res.push_back("coefficients"_nm = coef);
  fit_res.push_back("rank"_nm = rank);

  return fit_res;
}
