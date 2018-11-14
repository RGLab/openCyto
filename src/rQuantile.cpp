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

#include "faust.h"
#include <iostream>
#include <vector>
#include <algorithm>
#include <cmath>


std::vector<double> rQuantile(const std::vector<double>& dataVector, std::vector<double> probs) {
  //compute sample quantiles of probs in dataVector.
  //this is a restricted port of the R quantile function -- only bringing over type == 7
  std::vector<double> dVec = dataVector;
  std::sort(dVec.begin(),dVec.end());
  auto n = dVec.size();
  std::vector<double> index;
  for (auto p : probs)
    index.push_back(((n-1.0)*p));
  std::vector<double> lowInd = index;
  std::vector<double> highInd = index;
  
  for (auto i = 0; i != lowInd.size(); ++i) {
    lowInd[i] = std::floor(lowInd[i]);
    highInd[i] = std::ceil(highInd[i]);
  }
  std::vector<double> quantiles(lowInd.size(),0.0);
  std::vector<bool> offSet(index.size(),false);
  for (auto i = 0; i != lowInd.size(); ++i) {
    quantiles[i] = dVec[lowInd[i]];
    if ((index[i] > lowInd[i]) && (dVec[highInd[i]] != quantiles[i])) 
      offSet[i] = true;
  }
  double hRem;
  for (auto i = 0; i != offSet.size(); ++i) 
    if (offSet[i]) {
      hRem = index[i]-lowInd[i];
      quantiles[i] = ((1-hRem) * quantiles[i]) + (hRem * dVec[highInd[i]]);
    }
  return quantiles;
}
