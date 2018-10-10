#include <iostream>
#include <vector>
#include <cstdlib>
#include <cmath>
#include <algorithm>
#include <numeric>
#include <Rcpp.h>
#include "faust.h"


struct stringPairs {
  int startX;
  int endX;
  double stringVal;
};

// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::export]]
Rcpp::List getTautStringApprox(std::vector<double> dataVec) {
  //get the taut-string approximation, and reduce it to pairs of points
  //for compressed plotting.
  std::sort(dataVec.begin(),dataVec.end());
  stringInfo tautString = cpPmden(dataVec); 
  std::vector<double> fullString = tautString.string;
  std::vector<stringPairs> plotSpecs;
  stringPairs tmpSP;
  int currentSP = 0; // staring point of string
  int currentEP; //ending point
  double currentSV = fullString[0]; //string value
  bool accountForEnd = true;
  int epchk = fullString.size() - 1;
  for (int i = 1; i != fullString.size(); ++i) {
    if (fullString[i] != currentSV) {
      currentEP = i;
      if (i == epchk) {
	accountForEnd = false;
      }
      tmpSP.startX = currentSP;
      tmpSP.endX = currentEP;
      tmpSP.stringVal = currentSV;
      plotSpecs.push_back(tmpSP);
      currentSV = fullString[i];
      currentSP = currentEP;
    }
    if ((i == epchk) && (accountForEnd)) {
      currentEP = i;
      tmpSP.startX = currentSP;
      tmpSP.endX=currentEP;
      tmpSP.stringVal=currentSV;
      plotSpecs.push_back(tmpSP);
    }
  }
  Rcpp::Rcout << plotSpecs.size() << std::endl;
  //unpack the specs to pass to R
  std::vector<double> outX, outY;
  for (int i  = 0; i != plotSpecs.size(); ++i) {
    tmpSP = plotSpecs[i];
    currentSP = tmpSP.startX;
    currentEP = tmpSP.endX;
    currentSV = tmpSP.stringVal;
    outX.push_back(dataVec[currentSP]);
    outY.push_back(currentSV);
    outX.push_back(dataVec[currentEP]);
    outY.push_back(currentSV);
  }
  Rcpp::List outList = Rcpp::List::create(Rcpp::Named("stringX")=outX,
					  Rcpp::Named("stringY")=outY);
  return outList;
}

  
  
  
