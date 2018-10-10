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

/*
 *This file is a modified copy of the file dip.c from diptest R library. 
 *I belive that makes it a derived work according to GPL-3.
 *The orginal work is found in version 0.75-7 of the diptest library on CRAN.
 *The maintainer of that library, Martin Maechler <maechler at stat.math.ethz.ch>, is unaffiliated with the faust project.
 *Minor changes have been made for debugging within faust. 
 *--Evan Greene <egreene@fredhutch.org>, 3/15/2018
 */
#include <Rcpp.h>
#include <vector>
#include <iostream>
#include "faust.h"


void cppdip(const double x[],
	    const int *n_,
	    double *dip,
	    std::vector<int>& lo_hi,
	    int *ifault,
	    int *gcm,
	    int *lcm,
	    int *mn,
	    int *mj,
	    const int *min_is_0,
	    const int *debug) {
  int low = lo_hi[0];
  int high = lo_hi[1];
  int l_gcm = lo_hi[2];
  int l_lcm = lo_hi[3];

  const int n = *n_;
  int mnj, mnmnj, mjk, mjmjk, ig, ih, iv, ix,  i, j, k;
  double dip_l, dip_u, dipnew;
  /* Parameter adjustments, so I can do "as with index 1" : x[1]..x[n] */
  --mj;
  --mn;
  --lcm;
  --gcm;
  --x;

  /*-------- Function Body ------------------------------ */
  *ifault = 1;
  if (n <= 0) {
    return;
  }
  *ifault = 0;

  /* Check that X is sorted --- if not, return with  ifault = 2*/
  *ifault = 2;
  for (k = 2; k <= n; ++k) {
    if (x[k] < x[k - 1]) {
      return;
    }
  }
  *ifault = 0;

  /* Check for all values of X identical, */
  /*     and for 1 <= n < 4. */

  /* LOW contains the index of the current estimate  of the lower end
     of the modal interval, HIGH contains the index for the upper end.
  */
  low = 1;
  high = n; /*-- IDEA:  *xl = x[low];    *xu = x[high]; --*/

  /* M.Maechler -- speedup: it saves many divisions by n when we just work with
   * (2n * dip) everywhere but the very end! */
  *dip = (*min_is_0) ? 0. : 1.;
  if (n < 2 || x[n] == x[1]) {
    *dip /= (2*n);
    return;
  }
  /* Establish the indices   mn[1..n]  over which combination is necessary
     for the convex MINORANT (GCM) fit.
  */
  mn[1] = 1;
  for (j = 2; j <= n; ++j) {
    mn[j] = j - 1;
    while(1) {
      mnj = mn[j];
      mnmnj = mn[mnj];
      if (mnj == 1 ||
	  ( x[j]  - x[mnj]) * (mnj - mnmnj) <
	  (x[mnj] - x[mnmnj]) * (j - mnj)) break;
      mn[j] = mnmnj;
    }
  }

  /* Establish the indices   mj[1..n]  over which combination is necessary
     for the concave MAJORANT (LCM) fit.
  */
  mj[n] = n;
  for (k = n - 1; k >= 1; k--) {
    mj[k] = k + 1;
    while(1) {
      mjk = mj[k];
      mjmjk = mj[mjk];
      if (mjk == n ||
	  ( x[k]  - x[mjk]) * (mjk - mjmjk) <
	  (x[mjk] - x[mjmjk]) * (k - mjk)) break;
      mj[k] = mjmjk;
    }
  }

  bool stillDipping = true;
  /* ----------------------- Start the cycling. ------------------------------- */
  while (stillDipping) {
    /* Collect the change points for the GCM from HIGH to LOW. */
    gcm[1] = high;
    for(i = 1; gcm[i] > low; i++)
	gcm[i+1] = mn[gcm[i]];
    ig = l_gcm = i; // l_gcm == relevant_length(GCM)
    ix = ig - 1; //  ix, ig  are counters for the convex minorant.

    /* Collect the change points for the LCM from LOW to HIGH. */
    lcm[1] = low;
    for(i = 1; lcm[i] < high; i++)
	lcm[i+1] = mj[lcm[i]];
    ih = l_lcm = i; // l_lcm == relevant_length(LCM)
    iv = 2; //  iv, ih  are counters for the concave majorant.
    /*	Find the largest distance greater than 'DIP' between the GCM and
     *	the LCM from LOW to HIGH. */

    // FIXME: <Rconfig.h>  should provide LDOUBLE or something like it
    long double d = 0.;// <<-- see if this makes 32-bit/64-bit difference go..
    if (l_gcm != 2 || l_lcm != 2) {
      do { /* gcm[ix] != lcm[iv]  (after first loop) */
	long double dx;
	int gcmix = gcm[ix],
	  lcmiv = lcm[iv];
	if (gcmix > lcmiv) {
	  /* If the next point of either the GCM or LCM is from the LCM,
	   * calculate the distance here. */
	  int gcmi1 = gcm[ix + 1];
	  dx = (lcmiv - gcmi1 + 1) -
	    ((long double) x[lcmiv] - x[gcmi1]) * (gcmix - gcmi1)/(x[gcmix] - x[gcmi1]);
	  ++iv;
	  if (dx >= d) {
	    d = dx;
	    ig = ix + 1;
	    ih = iv - 1;
	    //if(*debug >= 2) Rprintf(" L(%d,%d)", ig,ih);
	  }
	}
	else {
	  /* If the next point of either the GCM or LCM is from the GCM,
	   * calculate the distance here. */
	  int lcmiv1 = lcm[iv - 1];
	  /* Fix by Yong Lu {symmetric to above!}; original Fortran: only ")" misplaced! :*/
	  dx = ((long double)x[gcmix] - x[lcmiv1]) * (lcmiv - lcmiv1) /
	    (x[lcmiv] - x[lcmiv1])- (gcmix - lcmiv1 - 1);
	  --ix;
	  if (dx >= d) {
	    d = dx;
	    ig = ix + 1;
	    ih = iv;
	  }
	}
	if (ix < 1)	ix = 1;
	if (iv > l_lcm)	iv = l_lcm;
      } while (gcm[ix] != lcm[iv]);
    }
    else { /* l_gcm or l_lcm == 2 */
      d = (*min_is_0) ? 0. : 1.;
    }
    if (d < *dip) {
      stillDipping = false;
    }
    else {
      /*     Calculate the DIPs for the current LOW and HIGH. */
      int j_best, j_l = -1, j_u = -1;
      /* The DIP for the convex minorant. */
      dip_l = 0.;
      for (j = ig; j < l_gcm; ++j) {
	double max_t = 1.;
	int j_ = -1, jb = gcm[j + 1], je = gcm[j];
	if (je - jb > 1 && x[je] != x[jb]) {
	  double C = (je - jb) / (x[je] - x[jb]);
	  for (int jj = jb; jj <= je; ++jj) {
	    double t = (jj - jb + 1) - (x[jj] - x[jb]) * C;
	    if (max_t < t) {
	      max_t = t; j_ = jj;
	    }
	  }
	}
	if (dip_l < max_t) {
	  dip_l = max_t; j_l = j_;
	}
      }

      /* The DIP for the concave majorant. */
      dip_u = 0.;
      for (j = ih; j < l_lcm; ++j) {
	double max_t = 1.;
	int j_ = -1, jb = lcm[j], je = lcm[j + 1];
	if (je - jb > 1 && x[je] != x[jb]) {
	  double C = (je - jb) / (x[je] - x[jb]);
	  for (int jj = jb; jj <= je; ++jj) {
	    double t = (x[jj] - x[jb]) * C - (jj - jb - 1);
	    if (max_t < t) {
	      max_t = t; j_ = jj;
	    }
	  }
	}
	if (dip_u < max_t) {
	  dip_u = max_t; j_u = j_;
	}
      }
      /* Determine the current maximum. */
      if(dip_u > dip_l) {
	dipnew = dip_u; j_best = j_u;
      } else {
	dipnew = dip_l; j_best = j_l;
      }
      if (*dip < dipnew) {
	*dip = dipnew;
	if(*debug) {
	  Rcpp::Rcout << "-> new larger dip " << dipnew << "(j_best = " << j_best << ")" << std::endl;
	}
      }
      else if(*debug) {
	Rcpp::Rcout << std::endl;
      }
      /*--- The following if-clause is NECESSARY  (may loop infinitely otherwise)!
	--- Martin Maechler, Statistics, ETH Zurich, July 30 1994 ---------- */
      if (low == gcm[ig] && high == lcm[ih]) {
	stillDipping = false;
	if(*debug) {
	  Rcpp::Rcout << "No improvement in  low = " << low << "high = " << high << " --> END" << std::endl;
	}
      } else {
	low  = gcm[ig];
	high = lcm[ih];	   
      }
    }
  }
  /* do this in the caller :
   *   *xl = x[low];  *xu = x[high];
   * rather return the (low, high) indices -- automagically via lo_hi[]  */
  *dip /= (2*n);
  return;
}
