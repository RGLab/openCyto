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
//   sqrt(sum((old - new)^2)/max(1e-20, sum(old^2)))
double irls_delta(cpp11::doubles residual_old, cpp11::doubles residual_new){
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
cpp11::list rlm_cpp(cpp11::doubles_matrix x, cpp11::doubles y, cpp11::writable::doubles residual_old, int maxit){
  using namespace cpp11::literals; // so we can use ""_nm syntax

  cpp11::writable::list res;
  
  int n = residual_old.size();
  vector<double> residual_abs(n);
  cpp11::writable::doubles w(n);
  bool done = false;  
  while(maxit-- > 0) {
    //compute scale
    for(int i = 0; i < n; i++)
      residual_abs[i] = abs(residual_old[i]);
    double scale = boost::math::statistics::median(residual_abs)/0.6745;
    
     
    if(scale == 0) {
      done = true;
      break;
    }
    
    //compte huber weight
    for(int i = 0; i < n; i++)
      w[i] = residual_old[i]/scale;
    w = psi_huber(w);
    
    //fit lm module
    auto lm_wfit = cpp11::package("stats")["lm.wfit"];
    res = cpp11::as_cpp<cpp11::writable::list>(lm_wfit(x, y, w, "method"_nm = "qr"));
    auto residual_new = cpp11::as_cpp<cpp11::doubles>(res["residuals"]);
    
    //check if converge
    auto conv = irls_delta(residual_old, residual_new);
    double acc = 1e-4;
    done = conv <= acc;
    res.push_back("done"_nm = done);
    res.push_back("w"_nm = w);

    if(done)
      break;
    else
      residual_old = residual_new;
  }
  return res;
}

