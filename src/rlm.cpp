#include <cpp11.hpp>
#include <vector>
using namespace std;

[[cpp11::register]]
cpp11::list rlm_cpp(cpp11::doubles_matrix x, cpp11::doubles y, cpp11::doubles w){
  using namespace cpp11::literals; // so we can use ""_nm syntax
  
  auto lm_wfit = cpp11::package("stats")["lm.wfit"];
  auto res = lm_wfit(x, y, w, "method"_nm = "qr");
  return cpp11::as_cpp<cpp11::list>(res);
}