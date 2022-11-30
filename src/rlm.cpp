#if defined(WIN32) || defined(_WIN32) || defined(__WIN32__) || defined(__NT__)
#define BOOST_DISABLE_THREADS
#endif
#include <cpp11.hpp>
#include <vector>
using namespace std;
#include "boost/math/statistics/univariate_statistics.hpp"
#include "boost/accumulators/statistics/median.hpp"
namespace bstacc = boost::accumulators;


cpp11::doubles psi_huber(cpp11::doubles u)
{
    auto k = 1.345;
    int n = u.size();
    cpp11::writable::doubles res(n);
    for(int i = 0; i < n; i++)
    {
      double x = k / abs(u[i]);
      res[i] = min(1.0, x);
      }
  return res;
}
double leastSquareDiff(cpp11::doubles residual_old, cpp11::doubles residual_new){
  double diff_square = 0, old_square = 0;
  
  for(int i = 0; i < residual_new.size(); i++)
  {
    auto x = residual_new[i] - residual_old[i];
    diff_square += x * x;
    old_square += residual_old[i] * residual_old[i];
    
    }
  return sqrt(diff_square/max(1e-20, old_square));
}

[[cpp11::register]]
cpp11::list rlm_cpp(cpp11::doubles_matrix<> x, cpp11::doubles y, int maxit){
  using namespace cpp11::literals; // so we can use ""_nm syntax
  auto lm_wfit = cpp11::package("stats")["lm.wfit"];
  
  cpp11::writable::list fit_res;
  
  //initial model fit
  int n = y.size();
  cpp11::writable::doubles w(n);//weight
  for(int i = 0; i < n; i++)
    w[i] = 1;
  fit_res = cpp11::as_cpp<cpp11::writable::list>(lm_wfit(x, y, w, "method"_nm = "qr"));
  cpp11::writable::doubles residual_old = cpp11::as_cpp<cpp11::doubles>(fit_res["residuals"]);
    
  //update fitted model iteratively with re-weighted least squares 
  vector<double> residual_abs(n);
  bool done = false;  
  double scale;
  while(maxit-- > 0)
  {
    //compute scale
    for(int i = 0; i < n; i++)
      residual_abs[i] = abs(residual_old[i]);
    scale = boost::math::statistics::median(residual_abs)/0.6745;
     
    if(scale == 0) {
      done = true;
      break;
    }
    
    //compute huber weight
    for(int i = 0; i < n; i++)
      w[i] = residual_old[i]/scale;
    w = psi_huber(w);
    
    //fit lm module
    fit_res = cpp11::as_cpp<cpp11::writable::list>(lm_wfit(x, y, w, "method"_nm = "qr"));
    auto residual_new = cpp11::as_cpp<cpp11::doubles>(fit_res["residuals"]);
    
    //check if converge
    auto conv = leastSquareDiff(residual_old, residual_new);
    done = conv <= 1e-4;
    
    if(done)
      break;
    else
      residual_old = residual_new;
  }
  
  if(!done)
    cpp11::warning("'rlm' failed to converge in" + to_string(maxit) + "  steps");
    
  fit_res.push_back("done"_nm = done);
  fit_res.push_back("w"_nm = w);
  fit_res.push_back("scale"_nm = scale);
  
  // cpp11::writable::doubles fitted_y(n);
  //TODO: to handle cases where x has extra columns of ssc-a and ssc-a/fsc-a
  // basically need to do x %*% coef in C
  // auto coefficients = cpp11::as_cpp<cpp11::doubles>(fit_res["coefficients"]);
  // double slope = coefficients[1];
  // double intercept = coefficients[0];
  // for(int i = 0; i < n; i++)
  //   fitted_y[i] = x(i, 1) * slope + intercept;
  // fit_res.push_back("fitted"_nm = fitted_y);
  //   
  return fit_res;
}

