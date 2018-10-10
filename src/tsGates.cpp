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
#include <Rcpp.h>
#include "faust.h"
#include "kMedDP.h"

bool isLocalMinTS(double &x_cand,double &x_left,double &x_right) {
  //inequality strict because the string jumps to discrete values
  if ((x_cand < x_left) && (x_cand < x_right))  return true;
  else return false;
}

std::vector<double> findKmedGates(const std::vector<double> &dataVector,
			      const std::vector<int> &classVector,
			      const int numClasses)
{
  //assumes classVector.size() == dataVector.size()
  std::vector<double> gates, tmpData;
  gates.push_back(*std::min_element(dataVector.begin(),dataVector.end()));
  for (int i = 0; i != numClasses; ++i) {
    tmpData.clear();
    for (auto j = 0; j != classVector.size(); j++)
      if (classVector[j] == i)
	tmpData.push_back(dataVector[j]);
    gates.push_back(*std::max_element(tmpData.begin(),tmpData.end()));
  }
  gates.push_back(*std::max_element(dataVector.begin(),dataVector.end()));
  return gates;
}

// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::export]]
std::vector<double> tsGates(const std::vector<double> &xVec, int modePrior) {
  //Given sorted univarite data vector xVec, returns an estimate of cut points to classify data
  //Assumes, at minimum, the DIP test has indicated multimodality in xData
    
  //If modePrior > 0, this indicates multi-dipping has been used to estimate number of modes
  //In this event, we use the estimate via k-meidan cut, kMedDP

  //In the rare event that the taut string estimates 1 mode while the DIP determines multimodality
  //Assume that there are two clusters, and use KMEANS to split data into two clusters.
  stringInfo tautString = cpPmden(xVec); 
  int numModes = tautString.nmax-(tautString.nmax/2); //always check to see if the taut string sees multimodality.
                                                      //Note we are deliberately taking avantage of integer truncation
                                                      //No need to call std::floor.
  int numCenters;
  std::vector<double> gates; 
  if ((numModes <= 1) || (modePrior != 0)) {
    if (modePrior == 0) numCenters = 2; 
    else numCenters = modePrior;
    std::vector<int> kClusters =  kMedDP(xVec,numCenters);
    int kEstimate = *std::max_element(kClusters.begin(),kClusters.end());
    std::vector<double> newGates = findKmedGates(xVec,kClusters,kEstimate);
    return newGates;
  }
  else {
    std::vector<double> ys = tautString.string;
    std::vector<double> yvals = ys; //copy string into yvals
    auto it = std::unique(yvals.begin(),yvals.end()); //collapse to unique elements
    yvals.resize(std::distance(yvals.begin(),it)); //resize vector
    std::vector<bool> localMins(yvals.size(),false);
    for (auto i = 1; i != (localMins.size()-1); ++i) {
      localMins[i] = isLocalMinTS(yvals[i],yvals[(i-1)],yvals[(i+1)]);
    }
    
    std::vector<double> cutValues;
    for (auto i = 0; i != yvals.size(); ++i) 
      if (localMins[i] == true)
	cutValues.push_back(yvals[i]);

    std::vector<double> newGates;
    std::vector<decltype(ys.size())> gateIndices;
    std::vector<decltype(ys.size())> yInds;
    decltype(ys.size()) cutIndex;
    long long indexSum = 0;
    double cCut; 
    for (auto j = 0; j != cutValues.size(); ++j) {
      yInds.clear();
      cCut = cutValues[j];
      for (auto i = 0; i != ys.size(); ++i) 
	if (ys[i] == cCut)
	  yInds.push_back(i);
      indexSum = 0;
      for (auto indexVal : yInds)
	indexSum += indexVal;
      cutIndex = indexSum/yInds.size();
      gateIndices.push_back(cutIndex);
    }
    newGates.push_back(*std::min_element(xVec.begin(),xVec.end()));
    for (auto ci : gateIndices)
      newGates.push_back(xVec[ci]);
    newGates.push_back(*std::max_element(xVec.begin(),xVec.end()));
    std::sort(newGates.begin(),newGates.end()); //ensure newGates are returned in ascending order
    return newGates;
  }
}
