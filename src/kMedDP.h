/* Ckmeans_1d_dp.h --- Head file for Ckmeans.1d.dp
 *  Declare wrap function "kmeans_1d_dp()"
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
 Haizhou Wang
 Computer Science Department
 New Mexico State University
 hwang@cs.nmsu.edu

 Created: Oct 10, 2010
 */

#ifndef KMEDIAN_1D_LIBRARY_H
#define KMEDIAN_1D_LIBRARY_H

#include <cstddef> // For size_t
#include <vector>
#include <string>
#include "within_cluster.h"

std::vector<int> kMedDP(const std::vector<double> & x,
			size_t K_estimate);
std::vector<int> kMeanDP(const std::vector<double> & x,
			size_t K_estimate);

void fill_dp_matrix(
    const std::vector<double> & x,
    const std::vector<double> & w,
    std::vector< std::vector< ldouble > > & S,
    std::vector< std::vector< size_t > > & J,
    const std::string & method,
    const enum DISSIMILARITY criterion);

void backtrack(const std::vector<double> & x,
               const std::vector<std::vector<size_t>> & J,
               std::vector<int> & cluster,
	       std::vector<double> & centers,
	       std::vector<double> & withinss,
               std::vector<double> & count);

/*
void backtrack(
    const std::vector<double> & x,
    const std::vector< std::vector< size_t > > & J,
    int* cluster, double* centers, double* withinss,
    double* count); //int* count);
*/
void backtrack_L1(
    const std::vector<double> & x,
    const std::vector< std::vector< size_t > > & J,
    std::vector<int> & cluster,
    std::vector<double> & centers,
    std::vector<double> & withinss,
    std::vector<double> & count);
/*
void backtrack_L1(
    const std::vector<double> & x,
    const std::vector< std::vector< size_t > > & J,
    int* cluster,
    double* centers,
    double* withinss,
    double* count);
*/
/*
void backtrack(
    const std::vector<double> & x,
    const std::vector< std::vector< size_t > > & J,
    std::vector<size_t> & count);

void backtrack_L2Y(
    const std::vector<double> & x, const std::vector<double> & y,
    const std::vector< std::vector< size_t > > & J,
    int* cluster, double* centers, double* withinss,
    double* count); //int* count);
*/
void fill_row_q_SMAWK(
    int imin, int imax, int q,
    std::vector< std::vector<ldouble> > & S,
    std::vector< std::vector<size_t> > & J,
    const std::vector<ldouble> & sum_x,
    const std::vector<ldouble> & sum_x_sq,
    const std::vector<ldouble> & sum_w,
    const std::vector<ldouble> & sum_w_sq,
    const enum DISSIMILARITY criterion);

/*
void fill_row_q(
    int imin, int imax, int q,
    std::vector< std::vector<ldouble> > & S,
    std::vector< std::vector<size_t> > & J,
    const std::vector<ldouble> & sum_x,
    const std::vector<ldouble> & sum_x_sq,
    const std::vector<ldouble> & sum_w,
    const std::vector<ldouble> & sum_w_sq,
    const enum DISSIMILARITY criterion);


void fill_row_q_log_linear(
    int imin, int imax, int q, int jmin, int jmax,
    std::vector< std::vector<ldouble> > & S,
    std::vector< std::vector<size_t> > & J,
    const std::vector<ldouble> & sum_x,
    const std::vector<ldouble> & sum_x_sq,
    const std::vector<ldouble> & sum_w,
    const std::vector<ldouble> & sum_w_sq,
    const enum DISSIMILARITY criterion
);
*/
/* One-dimensional cluster algorithm implemented in C++ */
/* x is input one-dimensional vector and
 Kmin and Kmax stand for the range for the number of clusters*/
/*
void kmeans_1d_dp(
    const double *x, const size_t N,
    const double * y,
    size_t Kmin, size_t Kmax,
    int* cluster, double* centers,
    double* withinss, double * size, // int* size,
    double* BIC,
    const std::string & estimate_k,
    const std::string & method,
    const enum DISSIMILARITY criterion);
*/
/*
void backtrack(
    const std::vector<double> & x,
    const std::vector< std::vector< size_t > > & J,
    std::vector<size_t> & counts, const int K);
*/
/*
size_t select_levels(
    const std::vector<double> & x,
    const std::vector< std::vector< size_t > > & J,
    size_t Kmin, size_t Kmax, double *BIC);

size_t select_levels_3_4_12(
    const std::vector<double> & x,
    const std::vector< std::vector< size_t > > & J,
    size_t Kmin, size_t Kmax, double *BIC);

void fill_weighted_dp_matrix(
    const std::vector<double> & x,
    const std::vector<double> & y,
    std::vector< std::vector< ldouble > > & S,
    std::vector< std::vector< size_t > > & J);

void backtrack_weighted(
    const std::vector<double> & x, const std::vector<double> & y,
    const std::vector< std::vector< size_t > > & J,
    std::vector<size_t> & counts, std::vector<double> & weights,
    const int K);

void backtrack_weighted(
    const std::vector<double> & x, const std::vector<double> & y,
    const std::vector< std::vector< size_t > > & J,
    int* cluster, double* centers, double* withinss,
    double* weights); //int* weights);

size_t select_levels_weighted(
    const std::vector<double> & x, const std::vector<double> & y,
    const std::vector< std::vector< size_t > > & J,
    size_t Kmin, size_t Kmax, double *BIC);

size_t select_levels_weighted_3_4_12(
    const std::vector<double> & x, const std::vector<double> & y,
    const std::vector< std::vector< size_t > > & J,
    size_t Kmin, size_t Kmax, double *BIC);


void range_of_variance(
    const std::vector<double> & x,
    double & variance_min, double & variance_max);
*/

#endif
