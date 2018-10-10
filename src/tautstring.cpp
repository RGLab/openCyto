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

/*
 *This file is a simple port of the file strings.c from the ftnonpar library, by Laurie Davies and Arne Kovac.
 *I belive that makes it a derived work according to GPL-3.
 *The orginal work is found in version 0.1-88 of that package.
 *The maintainer of that library, Arne Kovac <A.Kovac at bristol.ac.uk>, is not affiliated with this package.
 *Minor changes have for use within faust. 
 *--Evan Greene <egreene@fredhutch.org>, 9/19/2017, 3/15/2018
 */
#include <cstdlib>
#include <cmath>
#include <Rcpp.h>
#include <vector>
#include "faust.h"

/*
 *
 *
 *tautstring.cpp
 *
 *
 */

stringInfo tautString(const std::vector<double>& fdist, 
		      const std::vector<double>& t, 
		      const std::vector<double>& lower, 
		      const std::vector<double>& upper, 
		      double y1,
		      double yn,
		      long n,
		      int EXTRMEAN) 
{

  std::vector<double> i_taut_string((n-1),0.0);
  std::vector<int> knotsind(n,0);
  std::vector<double> knotst(n,0);
  std::vector<double> knotsy(n,0);
  int nknots=1;
  int nmax=0;

  double newmaxderiv,newminderiv,maxderiv=1e+38, minderiv=-maxderiv;
  int actind=2, maxind=1, minind=1, lastbound=0;
  std::vector<int> knotssign(n,0);

  lastbound=0;
  knotssign[0]=lastbound;
  knotsind[0]=1;
  knotsy[0]=y1; 
  knotst[0]=t[0]; 

  while(actind<=n) {
    if(actind<n) {
      newmaxderiv = (upper[(actind-1)]-knotsy[(nknots-1)])/(t[(actind-1)]-(knotst[(nknots-1)]));
      newminderiv = (lower[(actind-1)]-knotsy[(nknots-1)])/(t[(actind-1)]-(knotst[(nknots-1)]));
    }
    else {
      newmaxderiv = (yn-knotsy[(nknots-1)])/(t[(actind-1)]-(knotst[(nknots-1)]));
      newminderiv = (yn-knotsy[(nknots-1)])/(t[(actind-1)]-(knotst[(nknots-1)]));
    }
    if (newminderiv > maxderiv) {
      if(lastbound==-1) nmax++;
      knotssign[nknots]=1;
      knotsind[nknots] = maxind;
      knotsy[nknots] = upper[(maxind-1)];
      knotst[nknots]=t[(maxind-1)]; 
      nknots++;
      actind = maxind;
      lastbound=1;
      maxderiv = 1e+38;
      minderiv = -1e+38;
    }
    else if(newmaxderiv < minderiv) {
      if(lastbound==1) nmax++;
      knotsind[nknots] = minind;
      knotssign[nknots]=-1;
      knotsy[nknots] = lower[(minind-1)];
      knotst[nknots]=t[(minind-1)]; 
      nknots++;
      maxderiv = 1e+38;
      minderiv = -1e+38;
      lastbound=-1;
      actind = minind;
    }
    else {
      if(newmaxderiv < maxderiv) {
	maxderiv = newmaxderiv;
	maxind = actind;
      }
      if(newminderiv > minderiv) {
	minderiv = newminderiv;
	minind = actind;
      }
      if(actind==n) {
	if(lastbound!=0) {
	  knotsind[nknots] = actind;
	  knotsy[nknots] = yn;
	  knotst[nknots]=t[(actind-1)]; 
	  knotssign[nknots]=0; 
	  nknots++;
	}
	else {
	  lastbound=-1;
	  knotsind[0]=1;
	  knotssign[0]=-1;
	  knotsy[0]=y1; 
	  knotst[0]=t[0]; 
	  actind=1;
	  minind=1;
	  maxind=1;
	}
      }
    }
    actind++;
  }
 
  for(auto i=0;i<(nknots-1);i++) {
    for(auto j=knotsind[i];j<knotsind[i+1];j++) {
      if(knotssign[i]==knotssign[i+1]) {
	i_taut_string[j-1]=(knotsy[(i+1)]-knotsy[i])/(knotst[(i+1)]-knotst[i]);
      }
      else {
	if(EXTRMEAN) {
	  i_taut_string[j-1]=(fdist[knotsind[(i+1)]-1]-fdist[knotsind[i]-1])/(knotst[(i+1)]-knotst[i]);
	}
	else {
	  i_taut_string[j-1]=(knotsy[(i+1)]-knotsy[i])/(knotst[(i+1)]-knotst[i]);
	}
      }
    }
  }
  
  knotst.erase(std::remove(knotst.begin(), knotst.end(), 0), knotst.end()); //remove trailling zeros form knot locations
  std::vector<double> knotsy_out(knotst.size(),0.0);
  for (auto i = 0; i < knotst.size(); i++) {
    knotsy_out[i] = knotsy[i];
  }
  stringInfo results = {i_taut_string, knotsind, knotst, knotsy_out, nknots, nmax};
  return results;
}
