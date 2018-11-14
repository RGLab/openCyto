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



/* Ckmeans_1d_dp.cpp -- Performs 1-D k-means by a dynamic programming
 * approach that is guaranteed to be optimal.
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
 */

/*
 Joe Song
 Computer Science Department
 New Mexico State University
 joemsong@cs.nmsu.edu

 Haizhou Wang
 Computer Science Department
 New Mexico State University
 hwang@cs.nmsu.edu

 Created: May 19, 2007
 Updated: September 3, 2009
 Updated: September 6, 2009.  Handle special cases when K is too big or array
 contains identical elements.  Added number of clusters selection by the
 MCLUST package.
 Updated: Oct 6, 2009 by Haizhou Wang Convert this program from R code into
 C++ code.
 Updated: Oct 20, 2009 by Haizhou Wang, add interface to this function, so that
 it could be called directly in R

 Updated: March 20, 2014. Joe Song.
 1.  When the number of clusters is not uniquely specified, the code now
 automatically select an optimal number of clusters by Bayesian
 information criterion.
 2.  Reduced unnecessary sorting performed on the input string in
 kmeans_1d_dp().

 Updated: February 8, 2015.
 1.  Cleaned up the code by removing commented sections (Haizhou Wang)
 2.  Speed up the code slightly as suggested by a user (Joe Song)
 3.  Throw exceptions for fatal errors (Joe Song)

 Updated: May 3, 2016
 1. Changed from 1-based to 0-based C implementation (MS)
 2. Optimized the code by reducing overhead. See 22% reduction in runtime to
 cluster one million points into two clusters. (MS)
 3. Removed the result class ClusterResult

 Updated: May 7, 2016
 1. Substantial runtime reduction. Added code to check for an upper bound
 for the sum of within cluster square distances. This reduced the runtime
 by half when clustering 100000 points (from standard normal distribution)
 into 10 clusters.
 2. Eliminated the unnecessary calculation of (n-1) elements in the dynamic
 programming matrix that are not needed for the final result. This
 resulted in enormous reduction in run time when the number of cluster
 is 2: assigning one million points into two clusters took half a
 a second on iMac with 2.93 GHz Intel Core i7 processor.

 Updated: May 9, 2016
 1. Added an alternative way to fill in the dynamic programming
 matix by i (columns in the matrix) to achieve speedup

 Updated: May 21, 2016
 1. Moved the weighted solutions to new files. They will be updated separately

 Updated: May 26, 2016
 1. Implemented log-linear algorithm and tested successfully

 Updated: May 29, 2016
 1. Implemented linear algorithm but have bugs

 Updated: May 30, 2016
 1. Debugged the linear algorithm and tested successfully on all examples
 2. Added new test cases in the C++ testing project

 Updated: July 19, 2016.
 1. If the input array is already sorted, not sorting
 is performed.

 Updated: August 20, 2016
 1. The weighted univariate k-means now runs in O(kn), down from O(kn^2).
 This is a result of integrating weighted and unweighted k-means
 clustering into a unified dynamic programming function without sacrificing
 performance.

 */

/*
 * This program was modified by Evan Greene <egreene@fredhutch.org> on August 1, 2017/ March 15, 2018
 * The modification significantly reduced the functionality of this program, extracting the 
 * parts of the library used for 1-dimensional k-medians. The R-interface was removed,
 * since the reduction in functionality was to integrate this routine into the faust package.
 *
 * For the full functionality, download the R package "Ckmeans.1d.dp", by Joe Song and Haizhou Wang from CRAN.
 * Their comments are left in full at the head of the file.
 */


#include "kMedDP.h"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <string>
#include <vector>
#include <cassert>
#include <cstring>

std::vector<int> kMedDP(const std::vector<double> & x,
			size_t K_estimate)
{
  // Input:
  //  x -- an array of double precision numbers, sorted by calling function tsGates
  //  K_estimate -- the number of clusters
  // NOTE: All vectors in this program is considered starting at position 0.
  const size_t N = x.size();
  std::vector<int> cluster(N,0);
  std::vector<double> centers(K_estimate,0.0);
  std::vector<double> withinss(K_estimate,0.0);
  std::vector<double> size(K_estimate,0.0);

  std::vector<size_t> order(N);
  //Number generation using lambda function, not supported by all g++:
  std::size_t n(0);
  std::generate(order.begin(), order.end(), [&]{ return n++; });
  const enum DISSIMILARITY criterion = L1;
  std::vector<double> x_sorted(x);
  std::vector<double> y_sorted; 
  std::vector< std::vector< ldouble > > S( K_estimate, std::vector<ldouble>(N) );
  std::vector< std::vector< size_t > > J( K_estimate, std::vector<size_t>(N) );
  fill_dp_matrix(x_sorted, y_sorted, S, J, "linear", criterion);
  std::vector<int> cluster_sorted(N,0);
  backtrack_L1(x_sorted, J, cluster_sorted, centers, withinss, size);
  
  for(size_t i = 0; i < N; ++i) {
      // Obtain clustering on data in the original order
    cluster[order[i]] = cluster_sorted[i];
  }
  
  return cluster;
}
