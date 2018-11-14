//This file is part of faust, faster annotation using shape-constrained trees.     

//faust is free software: you can redistribute it and/or modify
//it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.

//faust is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.

//You should have received a copy of the GNU General Public License
//along with faust.  If not, see <http://www.gnu.org/licenses/>.

#include <Rcpp.h>
#include <vector>
#include "faust.h"

/*
 *
 *
 *cppApprox is a restricted port of the base R fucntion approx.
 *provides lienar interpolation.
 *
 */
std::vector<double> cppApprox(std::vector<double>& x, std::vector<double>& y, std::vector<double>& v) {
  long nxy = x.size();
  long nv = v.size();
  long i, j, ij;
  double cv, tmp;
  std::vector<double> rvec(v.size(),0.0);
  for (long k=0; k < nv; k++) {
    i = 0;
    j = nxy - 1;
    cv = v[k];
    if (cv < x[i]) {
      rvec[k] = y[0];
      continue;
    }
    if (cv > x[j]) {
      rvec[k]= y[nxy-1];
      continue;
    }
    /* find the correct interval by bisection */
    while(i < j - 1) { /* x[i] <= v <= x[j] */
      ij = (i + j)/2; /* i+1 <= ij <= j-1 */
      if(cv < x[ij]) j = ij; else i = ij;
      /* still i < j */
    }
    /* provably have i == j-1 */
    /* interpolation */
    if(cv == x[j]) {
      rvec[k] = y[j];
    }
    else if (cv == x[i]) {
      rvec[k] = y[i];
    }
    else {
  /* impossible: if(x[j] == x[i]) return y[i]; */
      tmp= y[i] + (y[j] - y[i]) * ((cv - x[i])/(x[j] - x[i]));
      rvec[k] = tmp;
    }
  }
  return(rvec);
}





