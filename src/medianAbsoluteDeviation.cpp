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
#include <algorithm>
#include <cmath>
#include "faust.h"

double medianAbsoluteDeviation(const std::vector<double>& dataVector) {
  std::vector<double> pv = {0.5};
  double median = (rQuantile(dataVector,pv)[0]);
  double madBase;
  double magicNumber = 1.4826; // approx 1/qnorm(3/4), ensures consistency of E[mad(X_1,...X_n)] = sigma  (see R documentation)
  std::vector<double> devs(dataVector.size(),0.0);
  for (auto i = 0; i != devs.size(); ++i) 
    devs[i] = std::abs(dataVector[i]-median);
  madBase = (rQuantile(devs,pv)[0]);
  return (madBase * magicNumber);
}
