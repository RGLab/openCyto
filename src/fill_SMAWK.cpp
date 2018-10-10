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

/* fill_SMAWK.cpp --- a divide-and-conquer algorithm to compute a
 *   row in the dynamic programming matrix in O(n) time.
 *
 * Copyright (C) 2016 Mingzhou Song
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as
 * published by the Free Software Foundation, either version 3 of the
 * License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program. If not, see <http://www.gnu.org/licenses/>.
 */
//
// Joe Song
// Created: July 15, 2016
// Modified:
//   July 28, 2016. Used stack to implement REDUCE algorithm.
//     Now faster than the log-linear implementation.
//   August 9, 2016. Removed the stacks and perform REDUCE
//     in place without additional data structures
//   August 10, 2016. Removed the relative index which is not necessary
//     to save during the computation
//   August 13, 2016. File name changed from fill_AKMSW.cpp to
//     fill_SMAWK.cpp

/*
 * This program was modified by Evan Greene <egreene@fredhutch.org> on August 1, 2017/March 15, 2018.
 * The modification significantly reduced the functionality of this program.
 * The reduction in functionality was to integrate this routine into the faust package.
 * I belive that makes it a derived work according to GPL-3.
 * For the full functionality, download the R package "Ckmeans.1d.dp", by Joe Song and Haizhou Wang from CRAN.
 * Their comments are left in full at the head of the file.
 */


#include <vector>
#include <algorithm>
#include <iostream>
#include <cassert>
#include <stack>

#include "kMedDP.h"

void reduce_in_place(int imin, int imax, int istep, int q,
                     const std::vector<size_t> & js,
                     std::vector<size_t> & js_red,
                     const std::vector< std::vector<ldouble> > & S,
                     const std::vector< std::vector<size_t> > & J,
                     const std::vector<ldouble> & sum_x,
                     const std::vector<ldouble> & sum_x_sq,
                     const std::vector<ldouble> & sum_w,
                     const std::vector<ldouble> & sum_w_sq,
                     const enum DISSIMILARITY criterion)
{
  int N = (imax - imin) / istep + 1;

  js_red = js;

  if(N >= js.size()) {
    return;
  }

  // Two positions to move candidate j's back and forth
  int left = -1; // points to last favorable position / column
  int right = 0; // points to current position / column

  size_t m = js_red.size();

  while(m > N) { // js_reduced has more than N positions / columns

    int p = (left + 1);

    int i = (imin + p * istep);
    size_t j = (js_red[right]);
    ldouble Sl = (S[q-1][j-1] +
      dissimilarity(criterion, j, i, sum_x, sum_x_sq, sum_w, sum_w_sq));

    size_t jplus1 = (js_red[right+1]);
    ldouble Slplus1 = (S[q-1][jplus1-1] +
      dissimilarity(criterion, jplus1, i, sum_x, sum_x_sq, sum_w, sum_w_sq));

    if(Sl < Slplus1 && p < N-1) {
      js_red[ ++ left ] = j; // i += istep;
      right ++; // move on to next position / column p+1
    } else if(Sl < Slplus1 && p == N-1) {
      js_red[ ++ right ] = j; // delete position / column p+1
      m --;
    } else { // (Sl >= Slplus1)
      if(p > 0) { // i > imin
        // delete position / column p and
        //   move back to previous position / column p-1:
        js_red[right] = js_red[left --];
      } else {
        right ++; // delete position / column 0
      }
      m --;
    }
  }

  for(int r=(left+1); r < m; ++r) {
    js_red[r] = js_red[right++];
  }

  js_red.resize(m);
  return;
}

inline void fill_even_positions
  (int imin, int imax, int istep, int q,
   const std::vector<size_t> & js,
   std::vector< std::vector<ldouble> > & S,
   std::vector< std::vector<size_t> > & J,
   const std::vector<ldouble> & sum_x,
   const std::vector<ldouble> & sum_x_sq,
   const std::vector<ldouble> & sum_w,
   const std::vector<ldouble> & sum_w_sq,
   const enum DISSIMILARITY criterion)
{
  // Derive j for even rows (0-based)
  size_t n = (js.size());
  int istepx2 = (istep << 1);
  size_t jl = (js[0]);
  for(int i=(imin), r(0); i<=imax; i+=istepx2) {
    while(js[r] < jl) {
      // Increase r until it points to a value of at least jmin
      r ++;
    }
    // Initialize S[q][i] and J[q][i]
    S[q][i] = S[q-1][js[r]-1] +
      dissimilarity(criterion, js[r], i, sum_x, sum_x_sq, sum_w, sum_w_sq);
    J[q][i] = js[r]; // rmin

    // Look for minimum S upto jmax within js
    int jh = (i + istep <= imax)
      ? J[q][i + istep] : js[n-1];

    int jmax = std::min((int)jh, (int)i);

    ldouble sjimin(
        dissimilarity(criterion, jmax, i, sum_x, sum_x_sq, sum_w, sum_w_sq)
      );

    for(++ r; r < n && js[r]<=jmax; r++) {

      const size_t & jabs = js[r];

      if(jabs > i) break;

      if(jabs < J[q-1][i]) continue;

      ldouble s =
        dissimilarity(criterion, jabs, i, sum_x, sum_x_sq, sum_w, sum_w_sq);
      ldouble Sj = (S[q-1][jabs-1] + s);

      if(Sj <= S[q][i]) {
        S[q][i] = Sj;
        J[q][i] = js[r];
      } else if(S[q-1][jabs-1] + sjimin > S[q][i]) {
        break;
      } /*else if(S[q-1][js[rmin]-1] + s > S[q][i]) {
 break;
      } */
    }
    r --;
    jl = jh;
  }
}

inline void find_min_from_candidates
  (int imin, int imax, int istep, int q,
   const std::vector<size_t> & js,
   std::vector< std::vector<ldouble> > & S,
   std::vector< std::vector<size_t> > & J,
   const std::vector<ldouble> & sum_x,
   const std::vector<ldouble> & sum_x_sq,
   const std::vector<ldouble> & sum_w,
   const std::vector<ldouble> & sum_w_sq,
   const enum DISSIMILARITY criterion)
{
  size_t rmin_prev = (0);

  for(int i=(imin); i<=imax; i+=istep) {

    size_t rmin = (rmin_prev);

    // Initialization of S[q][i] and J[q][i]
    S[q][i] = S[q-1][js[rmin] - 1] +
      dissimilarity(criterion, js[rmin], i, sum_x, sum_x_sq, sum_w, sum_w_sq);
      // ssq(js[rmin], i, sum_x, sum_x_sq, sum_w);
    J[q][i] = js[rmin];

    for(size_t r = (rmin+1); r<js.size(); ++r) {

      const size_t & j_abs = (js[r]);

      if(j_abs < J[q-1][i]) continue;
      if(j_abs > i) break;

      ldouble Sj = (S[q-1][j_abs - 1] +
        dissimilarity(criterion, j_abs, i, sum_x, sum_x_sq, sum_w, sum_w_sq));
        // ssq(j_abs, i, sum_x, sum_x_sq, sum_w));
      if(Sj <= S[q][i]) {
        S[q][i] = Sj;
        J[q][i] = js[r];
        rmin_prev = r;
      }
    }
  }
}

void SMAWK
  (int imin, int imax, int istep, int q,
   const std::vector<size_t> & js,
   std::vector< std::vector<ldouble> > & S,
   std::vector< std::vector<size_t> > & J,
   const std::vector<ldouble> & sum_x,
   const std::vector<ldouble> & sum_x_sq,
   const std::vector<ldouble> & sum_w,
   const std::vector<ldouble> & sum_w_sq,
   const enum DISSIMILARITY criterion)
{
#ifdef DEBUG_REDUCE
  Rcpp::Rcout << "i:" << '[' << imin << ',' << imax << ']' << '+' << istep
            << std::endl;
#endif

  if(imax - imin <= 0 * istep) { // base case only one element left

    find_min_from_candidates(
      imin, imax, istep, q, js, S, J, sum_x, sum_x_sq, sum_w,
      sum_w_sq, criterion
    );

  } else {

    // REDUCE

#ifdef DEBUG //_REDUCE
    Rcpp::Rcout << "js:";
    for (size_t l=0; l < js.size(); ++l) {
      Rcpp::Rcout << js[l] << ",";
    }
    Rcpp::Rcout << std::endl;
    Rcpp::Rcout << std::endl;
#endif

    std::vector<size_t> js_odd;

    reduce_in_place(imin, imax, istep, q, js, js_odd,
                    S, J, sum_x, sum_x_sq, sum_w,
                    sum_w_sq, criterion);

    int istepx2 = (istep << 1);
    int imin_odd = (imin + istep);
    int imax_odd = (imin_odd + (imax - imin_odd) / istepx2 * istepx2);

    // Recursion on odd rows (0-based):
    SMAWK(imin_odd, imax_odd, istepx2,
          q, js_odd, S, J, sum_x, sum_x_sq, sum_w,
          sum_w_sq, criterion);

#ifdef DEBUG // _REDUCE
    Rcpp::Rcout << "js_odd (reduced):";
    for (size_t l=0; l<js_odd.size(); ++l) {
      Rcpp::Rcout << js_odd[l] << ",";
    }
    Rcpp::Rcout << std::endl << std::endl;

    Rcpp::Rcout << "even pos:";
    for (int i=imin; i<imax; i+=istepx2) {
      Rcpp::Rcout << i << ",";
    }
    Rcpp::Rcout << std::endl << std::endl;

#endif

    fill_even_positions(imin, imax, istep, q, js,
                        S, J, sum_x, sum_x_sq, sum_w,
                        sum_w_sq, criterion);
  }
}

void fill_row_q_SMAWK(int imin, int imax, int q,
                      std::vector< std::vector<ldouble> > & S,
                      std::vector< std::vector<size_t> > & J,
                      const std::vector<ldouble> & sum_x,
                      const std::vector<ldouble> & sum_x_sq,
                      const std::vector<ldouble> & sum_w,
                      const std::vector<ldouble> & sum_w_sq,
                      const enum DISSIMILARITY criterion)
{
  // Assumption: each cluster must have at least one point.

  std::vector<size_t> js(imax-q+1);
  int abs = (q);
  std::generate(js.begin(), js.end(), [&] { return abs++; } );

  SMAWK(imin, imax, 1, q, js, S, J, sum_x, sum_x_sq, sum_w, sum_w_sq, criterion);
}
