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


/* dynamic_prog.cpp  --- dynamic programming function fill_dp_matrix and
 *   various versions of backtrack.
 *
 * Copyright (C) 2010-2016 Mingzhou Song and Haizhou Wang
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
 *
 * Joe Song
 * Created: August 20, 2016. Extracted from Ckmeans.1d.dp.cpp
 */

/*
 * This program was modified by Evan Greene <egreene@fredhutch.org> on March 15, 2018.
 * The modification significantly reduced the functionality of this program, extracting the 
 * parts of the library used for 1-dimensional k-medians. The R-interface was removed,
 * since the reduction in functionality was to integrate this routine into the faust package.
 * I belive that makes it a derived work according to GPL-3.
 * For the full functionality, download the R package "Ckmeans.1d.dp", by Joe Song and Haizhou Wang from CRAN.
 * Their comments are left in full at the head of the file.
 */




#include <string>
#include <iostream>
#include <cassert>
#include <cmath>
#include "kMedDP.h"
void fill_dp_matrix(const std::vector<double> & x, // data
                    const std::vector<double> & w, // weight
                    std::vector< std::vector< ldouble > > & S,
                    std::vector< std::vector< size_t > > & J,
                    const std::string & method,
                    const enum DISSIMILARITY criterion)
  /*
   x: One dimension vector to be clustered, must be sorted (in any order).
   S: K x N matrix. S[q][i] is the sum of squares of the distance from
   each x[i] to its cluster mean when there are exactly x[i] is the
   last point in cluster q
   J: K x N backtrack matrix

   NOTE: All vector indices in this program start at position 0
   */
{
  const int K = (int) S.size();
  const int N = (int) S[0].size();

  std::vector<ldouble> sum_x(N), sum_x_sq(N);
  std::vector<ldouble> sum_w(w.size()), sum_w_sq(w.size());

  std::vector<int> jseq;

  ldouble shift = x[N/2]; // median. used to shift the values of x to
  //  improve numerical stability

  if(w.empty()) { // equally weighted
    sum_x[0] = x[0] - shift;
    sum_x_sq[0] = (x[0] - shift) * (x[0] - shift);
  } else { // unequally weighted
    sum_x[0] = w[0] * (x[0] - shift);
    sum_x_sq[0] = w[0] * (x[0] - shift) * (x[0] - shift);
    sum_w[0] = w[0];
    sum_w_sq[0] = w[0] * w[0];
  }

  S[0][0] = 0;
  J[0][0] = 0;

  for(int i = 1; i < N; ++i) {

    if(w.empty()) { // equally weighted
      sum_x[i] = sum_x[i-1] + x[i] - shift;
      sum_x_sq[i] = sum_x_sq[i-1] + (x[i] - shift) * (x[i] - shift);
    } else { // unequally weighted
      sum_x[i] = sum_x[i-1] + w[i] * (x[i] - shift);
      sum_x_sq[i] = sum_x_sq[i-1] + w[i] * (x[i] - shift) * (x[i] - shift);
      sum_w[i] = sum_w[i-1] + w[i];
      sum_w_sq[i] = sum_w_sq[i-1] + w[i]*w[i];
    }

    // Initialize for q = 0
    S[0][i] = dissimilarity(criterion, 0, i, sum_x, sum_x_sq, sum_w, sum_w_sq); 
    J[0][i] = 0;
  }
  for(int q = 1; q < K; ++q) {
    int imin;
    if(q < K - 1) {
      imin = std::max(1, q);
    } else {
      // No need to compute S[K-1][0] ... S[K-1][N-2]
      imin = N-1;
    }
    fill_row_q_SMAWK(imin, N-1, q, S, J, sum_x, sum_x_sq, sum_w, sum_w_sq, criterion);
  }
}


void backtrack(const std::vector<double> & x,
               const std::vector< std::vector< size_t > > & J,
               std::vector<int>& cluster,
	       std::vector<double>& centers,
	       std::vector<double>& withinss,
               std::vector<double>& count) //int* count)
{
  const size_t K = J.size();
  const size_t N = J[0].size();
  size_t cluster_right = N-1;
  size_t cluster_left;

  // Backtrack the clusters from the dynamic programming matrix
  for(int q = ((int)K)-1; q >= 0; --q) {
    cluster_left = J[q][cluster_right];

    for(size_t i = cluster_left; i <= cluster_right; ++i)
      cluster[i] = q;

    double sum = 0.0;

    for(size_t i = cluster_left; i <= cluster_right; ++i)
      sum += x[i];

    centers[q] = sum / (cluster_right-cluster_left+1);

    for(size_t i = cluster_left; i <= cluster_right; ++i)
      withinss[q] += (x[i] - centers[q]) * (x[i] - centers[q]);

    count[q] = (int) (cluster_right - cluster_left + 1);

    if(q > 0) {
      cluster_right = cluster_left - 1;
    }
  }
}


void backtrack_L1(const std::vector<double> & x,
		  const std::vector< std::vector< size_t > > & J,
		  std::vector<int> & cluster,
		  std::vector<double> & centers,
		  std::vector<double> & withinss,
		  std::vector<double> & count) 
{
  const size_t K = J.size();
  const size_t N = J[0].size();
  size_t cluster_right = N-1;
  size_t cluster_left;

  // Backtrack the clusters from the dynamic programming matrix
  for(int q = ((int)K)-1; q >= 0; --q) {
    cluster_left = J[q][cluster_right];
    
    for(size_t i = cluster_left; i <= cluster_right; ++i)
      cluster[i] = q;

    centers[q] = x[(cluster_right+cluster_left) >> 1];

    for(size_t i = cluster_left; i <= cluster_right; ++i)
      withinss[q] += std::fabs(x[i] - centers[q]);
    
    count[q] = (int) (cluster_right - cluster_left + 1);
    
    if(q > 0) {
      cluster_right = cluster_left - 1;
    }
  }
}
