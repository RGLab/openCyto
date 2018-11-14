/* within_cluster.h---dissimilarity within a cluster
 *
 * Copyright (C) 2017 Mingzhou Song
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
// Created: March 2, 2017.
//   1. ssq() function (L2 norm) was extracted from Ckmeans.1d.dp.h
//   2. added sabs() function using L1 norm
//   3. added a general function within_cluster_dissimilarity

#ifndef KMEDIAN_WITHINC_H
#define KMEDIAN_WITHINC_H

#include <vector>
#include <string>

// typedef long double ldouble;
typedef double ldouble;

enum DISSIMILARITY {
  L1, L2, L2Y
};

inline ldouble ssq(
    const size_t j, const size_t i,
    const std::vector<ldouble> & sum_x, // running sum of xi
    const std::vector<ldouble> & sum_x_sq, // running sum of xi^2
    const std::vector<ldouble> & sum_w=std::vector<ldouble>(0)  // running sum of weights
)
{
  ldouble sji(0.0);

  if(sum_w.empty()) { // equally weighted version
    if(j >= i) {
      sji = 0.0;
    } else if(j > 0) {
      ldouble muji = (sum_x[i] - sum_x[j-1]) / (i - j + 1);
      sji = sum_x_sq[i] - sum_x_sq[j-1] - (i - j + 1) * muji * muji;
    } else {
      sji = sum_x_sq[i] - sum_x[i] * sum_x[i] / (i+1);
    }
  } else { // unequally weighted version
    if(sum_w[j] >= sum_w[i]) {
      sji = 0.0;
    } else if(j > 0) {
      ldouble muji = (sum_x[i] - sum_x[j-1]) / (sum_w[i] - sum_w[j-1]);
      sji = sum_x_sq[i] - sum_x_sq[j-1] - (sum_w[i] - sum_w[j-1]) * muji * muji;
    } else {
      sji = sum_x_sq[i] - sum_x[i] * sum_x[i] / sum_w[i];
    }
  }

  sji = (sji < 0) ? 0 : sji;
  return sji;
}

inline ldouble sabs(
    const size_t j, const size_t i,
    const std::vector<ldouble> & sum_x, // running sum of xi
    const std::vector<ldouble> & sum_w  // running sum of weights
)
{
  ldouble sji(0.0);

  if(sum_w.empty()) { // equally weighted version
    if(j >= i) {
      sji = 0.0;
    } else if(j > 0) {
      size_t l = (i + j) >> 1; // l is the index to the median of the cluster

      if(((i-j+1) % 2) == 1) {
        // If i-j+1 is odd, we have
        //   sum (x_l - x_m) over m = j .. l-1
        //   sum (x_m - x_l) over m = l+1 .. i
        sji = - sum_x[l-1] + sum_x[j-1] + sum_x[i] - sum_x[l];
      } else {
        // If i-j+1 is even, we have
        //   sum (x_l - x_m) over m = j .. l
        //   sum (x_m - x_l) over m = l+1 .. i
        sji = - sum_x[l] + sum_x[j-1] + sum_x[i] - sum_x[l];
      }
    } else { // j==0
      size_t l = i >> 1; // l is the index to the median of the cluster

      if(((i+1) % 2) == 1) {
        // If i-j+1 is odd, we have
        //   sum (x_m - x_l) over m = 0 .. l-1
        //   sum (x_l - x_m) over m = l+1 .. i
        sji = - sum_x[l-1] + sum_x[i] - sum_x[l];
      } else {
        // If i-j+1 is even, we have
        //   sum (x_m - x_l) over m = 0 .. l
        //   sum (x_l - x_m) over m = l+1 .. i
        sji = - sum_x[l] + sum_x[i] - sum_x[l];
      }
    }
  } else { // unequally weighted version
    // no exact solutions are known.
  }

  sji = (sji < 0) ? 0 : sji;
  return sji;
}

inline ldouble dissimilarity(
    const enum DISSIMILARITY dis,
    const size_t j, const size_t i,
    const std::vector<ldouble> & sum_x, // running sum of xi
    const std::vector<ldouble> & sum_x_sq, // running sum of xi^2
    const std::vector<ldouble> & sum_w,  // running sum of weights
    const std::vector<ldouble> & sum_w_sq  // running sum of square of weights
)
{
  ldouble d=0;

  switch(dis) {
  case L1:
    d = sabs(j, i, sum_x, sum_w);
    break;
  case L2:
    d = ssq(j, i, sum_x, sum_x_sq, sum_w);
    break;
  case L2Y:
    d = ssq(j, i, sum_w, sum_w_sq);
    break;
  }
  return d;
}

#endif
