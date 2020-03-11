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


#include <iostream>
#include <vector>
#include <cstdlib>
#include <cmath>
#include <algorithm>
#include <numeric>
#include "faust.h"

/*
 * The values for kuiper diff bounds taken from pmden.R in ftnonpar library.
 * In the source file pmden.R the stored values have 15 decimal digits. 
 * However, these values appear to be coerced to double when the library is loaded for use in an R sesssion.
 * To duplicate library performance, we use the truncated values, and simply store as double here.
 */
std::vector<double> kuip01 = {0.30750198,0.22031415,0.15749822,0.10121863,0.07178070,0.05080985,0.03215992};
std::vector<double> kuip02 = {0.26097572,0.18867970,0.13705738,0.08810506,0.06335109,0.04492857,0.02858946};
std::vector<double> kuip03 = {0.14897866,0.11094044,0.08092217,0.05257830,0.03748763,0.02696284,0.01707126};
std::vector<double> kuip04 = {0.12438262,0.09529240,0.06950244,0.04524029,0.03286277,0.02348616,0.01502098};
std::vector<double> kuip05 = {0.09319795,0.07219891,0.05339150,0.03538433,0.02598238,0.01851747,0.01174251};
std::vector<double> kuip06 = {0.08018187,0.06280868,0.04765458,0.03164710,0.02301769,0.01648256,0.01055136};
std::vector<double> kuip07 = {0.066165059,0.052643314,0.040499260,0.027003915,0.019671301,0.014140314,0.009111927};
std::vector<double> kuip08 = {0.05768921,0.04729112,0.03630651,0.02443004,0.01784079,0.01289449,0.00831476};
std::vector<double> kuip09 = {0.049727442,0.041269518,0.031941519,0.021631778,0.015981589,0.011540437,0.007469195};
std::vector<double> kuip10 = {0.044093021,0.037163268,0.029176417,0.019971113,0.014785146,0.010818339,0.006926491};
std::vector<double> kuip11 = {0.038637596,0.033532556,0.026499289,0.018379105,0.013643814,0.009856368,0.006386310};
std::vector<double> kuip12 = {0.035037285,0.030666602,0.024677670,0.017132430,0.012725232,0.009228679,0.005995483};
std::vector<double> kuip13 = {0.031642378,0.028084749,0.022794216,0.015823040,0.011929053,0.008645996,0.005649372};
std::vector<double> kuip14 = {0.029033131,0.025996755,0.021296584,0.014971161,0.011244505,0.008220232,0.005350641};
std::vector<double> kuip15 = {0.025883091,0.024057281,0.020032581,0.014186956,0.010545471,0.007767515,0.005061438};
std::vector<double> kuip16 = {0.023257645,0.022536446,0.018824915,0.013486222,0.010061760,0.007427072,0.004817209};
std::vector<double> kuip17 = {0.020675023,0.020836156,0.017644365,0.012805928,0.009588957,0.007102941,0.004595116};
std::vector<double> kuip18 = {0.019136349,0.019431774,0.016733216,0.012146658,0.009200609,0.006771573,0.004431340};
std::vector<double> kuip19 = {0.017891379,0.018348230,0.015944010,0.011620645,0.008794318,0.006483762,0.004265070};

/*
 * Store the individual kuiper vectors in a common lookup-vector 
 * (Don't have matrix row structure for cpPmeden, contrary to R implementation)
 */
std::vector<std::vector<double>> kuipdiffbounds = {kuip01,kuip02,kuip03,kuip04,kuip05,kuip06,kuip07,kuip08,kuip09,kuip10,
						   kuip11,kuip12,kuip13,kuip14,kuip15,kuip16,kuip17,kuip18,kuip19};

std::vector<double> kuipdiffbounds_x = {50.0, 100.0, 200.0, 500.0, 1000.0, 2000.0, 5000.0};

stringInfo cpPmden(const std::vector<double>& xIn) {
  //sorting is current turned off because cpPmden can only be called by tsGates, which is responsible for the sorting.
  //in the event cpPmden is to be used on a stand-alone basis, re-implement the sorting.
  //std::sort(std::begin(x), std::end(x)); //begin by sorting the data vector in place.
  std::vector<double> x = xIn;
  double locsq_factor = 0.9;
  
  //bool extrema_mean =true; //parameter passed to tautstring C routine. 
  int extrema_mean = 1; //parameter passed to tautstring C routine. 
  int maxkuipnr = 19;
  unsigned long nsamp = x.size();
  std::vector<double> currbounds(maxkuipnr,0.0);  //in this port, set maxkuipnr to 19. Do not allow user choice of mod.
	      
  if (nsamp > 5000)  {
    double kuiper_multiplier = (std::sqrt(5000.0)/std::sqrt(double(nsamp)));
    for (int i = 0; i < maxkuipnr; i++) 
      currbounds[i] = ((kuipdiffbounds[i])[6]) * kuiper_multiplier;
  }
  else {
    std::vector<double> tmp_nsamp(1,nsamp);
    for (int i = 0; i < maxkuipnr; i++) 
      currbounds[i] = (cppApprox(kuipdiffbounds_x, kuipdiffbounds[i], tmp_nsamp)[0]);
  }
  std::vector<double> dataemp(nsamp,(1/double(nsamp-1)));
  dataemp[0] = 0.0;
  double acc_min = 1.0;
  double max_diff = (x[(nsamp-1)]-x[0]);
  double current_ratio = 1.0;
  for (auto i = 0; i < (nsamp-1); i++) {
    current_ratio = (x[(i+1)]-x[i])/max_diff;
    if (current_ratio < acc_min)
      acc_min = current_ratio;
  }
  if (acc_min < 1e-14) {
    double accuracy = ((medianAbsoluteDeviation(xIn))/(1000.0)); 
    std::vector<double> newx;
    for (auto val : x) 
      newx.push_back(val/accuracy);
    decltype(nsamp) j = 0, k, cur_ind;
    std::vector<decltype(nsamp)> ind;
    double rhs, newaccx;
    while((j+1) < nsamp) {
      ind.clear();
      rhs = std::floor(newx[j]) + 1;
      for (auto i = j; i != nsamp; i++) 
	if (newx[i] < rhs)
	  ind.push_back(i);
      
      if (ind.size() == 0)
	k = j;
      else 
	k = *std::max_element(ind.begin(),ind.end());
      cur_ind = 1;
      for (auto i = j; i <= k; i++) {
	newaccx = std::floor(newx[i]) + double(cur_ind)/double(k-j+2);
	cur_ind++;
      }
      j++;
    }
    for (auto &i : newx)
      i = i * accuracy;
    x = newx;
  }

  /*initialize the uniform mesh*/
  std::vector<double> fdist_y(nsamp,0.0);
  std::iota(fdist_y.begin(),fdist_y.end(),0.0);
  for (auto i : fdist_y) {
    fdist_y[i] = fdist_y[i]/double((nsamp-1));
  }
  std::vector<double> eps(nsamp,0.5);
  eps[0] = 0.0;
  eps[(nsamp-1)]=0.0;

  std::vector<double> x_string((nsamp-1),0.0);
  std::vector<int> x_string_compare((nsamp-2),0);
  long first_jump = -1, last_jump = -1;
  std::vector<double> lower(nsamp,0.0);
  std::vector<double> upper(nsamp,0.0);
  stringInfo fts;
  std::vector<double> lastunif(nsamp,0.0);
  std::vector<double> differences(nsamp,0.0);
  std::vector<double> empirical_sum(nsamp,0.0);
  std::vector<double> currkkuip(maxkuipnr,0.0);
  std::vector<int> kuipinds(maxkuipnr,0);
  int first_kuiper_ind = -1;
  double cur_kuiper_val = 0.0;
  double tmp_eps = 0.0;
  long irmax, icomax;
  double prp, currsum;
  std::vector<int> kni(nsamp,0);
  long kni_sum = 0, num_kni_lookups = 0, cur_kni_lookup = 0;
  std::vector<long> kni_lookups(1,0);
  bool stillRepeating = true;
  bool anyJumpDetected = false;
  while (stillRepeating) {
    //first, update the tube around the empirical df
    for (auto i = 0; i < nsamp; i++) {
      lower[i] = fdist_y[i] - eps[i]; 
      upper[i] = fdist_y[i] + eps[i];
    }
    fts = tautString(fdist_y,x,lower,upper,(upper[0]),lower[(nsamp-1)],nsamp,extrema_mean);
    x_string = fts.string;
    anyJumpDetected = false;
    for (auto i=0; i < (nsamp-2); i++) {
      if (x_string[(i+1)] != x_string[i]) { 
	//capture the jump points on the first past, to remove double-checking later
	if (first_jump == -1) {
	  first_jump = i;
	}
	last_jump = i;
	anyJumpDetected = true;
      }
      else {
	x_string_compare[i] = 0;
      }
    }
    //comment out previous test
    //if (first_jump > 0) { //check for points where the taut string jumps, correct nmax count
    if (anyJumpDetected) {
      if (x_string[first_jump] > x_string[(first_jump + 1)]) {
	fts.nmax = fts.nmax + 1;
      }
      if (x_string[last_jump] < x_string[(last_jump + 1)]) {
	fts.nmax = fts.nmax + 1;
      }
    }
    //reset flags for next iteration
    first_jump = -1;
    last_jump = -1;
    //next, linear interploation bewteen knots of the taut string
    lastunif = cppApprox(fts.knotst, fts.knotsy, x);
    std::partial_sum(dataemp.begin(),dataemp.end(),empirical_sum.begin());
    for (auto i = 0; i < nsamp; i++) {
      differences[i] = empirical_sum[i] - lastunif[i];
    }
    //update kuiper bounds and check points where they exceed the current kuiper bounds
    currkkuip = kkuiper(differences,nsamp,maxkuipnr);
    if (currkkuip[0] > (currbounds[0] + 1e-08)) {
      kuipinds[0] = 1;
      first_kuiper_ind = 0;
    }
    else {
      kuipinds[0] = 0;
    }
    for (int i = 1; i < maxkuipnr; i++) {
      if ((currkkuip[i]-currkkuip[(i-1)])> (currbounds[i] + 1e-08)) {
	kuipinds[i] = 1;
	if (first_kuiper_ind == -1) {
	  first_kuiper_ind = i;
	}
      }
      else {
	kuipinds[i] = 0;
      }
    }
    if (first_kuiper_ind != -1) {
      cur_kuiper_val = (currbounds[first_kuiper_ind]/2);
      first_kuiper_ind = -1;//reset the flag for next iteration
      for (auto i = 0; i < nsamp; i++) {
	tmp_eps = eps[i];
	if (tmp_eps > 0)
	  eps[i] = cur_kuiper_val;
      }
    }
    else{
      first_kuiper_ind = -1;//reset the flag for next iteration
      irmax = std::floor(std::log(nsamp)) + 3;
      icomax = 1;
      prp = std::exp(-1);
      currsum = 2 * exp(-1);
      while (std::log(currsum) < (std::log(0.95)/(double(nsamp)))) {
	icomax = icomax + 1;
	prp = prp/(double(icomax));
	currsum = currsum + prp;
      }
      std::fill(kni.begin(),kni.end(),0.0);
      /*this appears to be spurious logic from pmden. since kni is set to zero every iteration of the while loop,
       *the diff checks can never be true. so, they are omitted from the port. we include the R logic for future reference.
       * ind1 <- diff(kni, lag = 1) > 0.04/(nsamp^2)
       * if (sum(ind1) > 0) 
       *   kni[c(FALSE, ind1) | c(ind1, FALSE)] <- 1
       * if (nsamp >= 3) {
       *   ind2 <- diff(kni, lag = 2) > 0.25/(nsamp^1.5)
       *   if (sum(ind2) > 0) 
       *     kni[c(FALSE, FALSE, ind2) | c(FALSE, ind2, FALSE) | c(ind2, FALSE, FALSE)] <- 1 
       * }
       */
      local_density(lastunif,kni,nsamp,icomax,irmax);
      kni_sum = 0; //set to zero every time through the loop
      kni_lookups.clear(); //empty lookups every time.
      for (auto i = 0; i < nsamp; i++) {
	if (kni[i] == 1) {
	  kni_sum += 1;
	  kni_lookups.push_back(i);
	}
      }
      if (kni_sum  == 0) {
	stillRepeating = false;
      }
      else {
	num_kni_lookups = kni_lookups.size();
	for (auto i = 0; i < num_kni_lookups; i++) {
	  cur_kni_lookup =  kni_lookups[i];
	  tmp_eps = eps[cur_kni_lookup] * locsq_factor;
	  eps[cur_kni_lookup] = tmp_eps;
	}
      }
    }
  }
  return fts;
}
    
