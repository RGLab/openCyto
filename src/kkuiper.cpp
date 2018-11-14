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
 *This file is a simple port of the file kkuip.c from the ftnonpar library, by Laurie Davies and Arne Kovac.
 *I belive that makes it a derived work according to GPL-3.
 *The orginal work is found in version 0.1-88 of that package.
 *The maintainer of that library, Arne Kovac <A.Kovac at bristol.ac.uk>, is not affiliated with this package.
 *Minor changes have for use within faust. 
 *--Evan Greene <egreene@fredhutch.org>, 9/19/2017, 3/15,2018
 */


#include <cstdlib>
#include <cmath>
#include <Rcpp.h>
#include <vector>
#include "faust.h"

/*
 *
 *
 *
 *kkuiper.cpp
 *
 *
 *
 */
void easymax(std::vector<double>& x, long n,long left,long right,long *erga,long *ergb,double *ergmax)
{
  double min,max;
  long i,mini,maxi;

  min=1e+38;
  max=-1e+38;
  for(i=left;i<=right;i++) {
    if(x[i]>max) { 
      max=x[i];
      maxi=i;
    }
    if(x[i]<min) {
      min=x[i];
      mini=i;
    }
  }
  if(mini<maxi) {
    *erga=mini;
    *ergb=maxi;
  }
  else {
    *ergb=mini;
    *erga=maxi;
  }
  *ergmax=max-min;
}

void difficultmax(std::vector<double>& x, 
		  long n, 
		  long left, 
		  long right, 
		  long *erga, 
		  long *ergb,
		  double *ergmax)
{
  double min, max;
  std::vector<double> mins(n,0.0);
  std::vector<double> maxs(n,0.0);
  long i, maxi;
  std::vector<int> minis(n,0);
  std::vector<int> maxis(n,0);
  
  min=1e+38;
  max=-1e+38;

  if(x[left]<x[right])
    {
      maxs[left]=x[left];
      maxis[left]=left;
      for(i=left+1;i<=right;i++)
        if(x[i]>maxs[i-1])
          {
	    maxs[i]=x[i];
	    maxis[i]=i;
          }
        else
          {
	    maxs[i]=maxs[i-1];
	    maxis[i]=maxis[i-1];
          }
      mins[right]=x[right];
      minis[right]=right;
      for(i=right-1;i>=left;i--)
        if(x[i]<mins[i+1])
          {
	    mins[i]=x[i];
	    minis[i]=i;
          }
        else
          {
	    mins[i]=mins[i+1];
	    minis[i]=minis[i+1];
          }
      max=-1e+38;
      for(i=left;i<=right;i++)
        if(maxs[i]-mins[i]>max)
          {
	    max=maxs[i]-mins[i];
            maxi=i;
          }
      *erga=maxis[maxi];
      *ergb=minis[maxi];
      *ergmax=std::fabs(x[*erga]-x[left])+std::fabs(x[right]-x[*ergb])-std::fabs(x[right]-x[left]);
    }
  else
    {
      mins[left]=x[left];
      minis[left]=left;
      for(i=left+1;i<=right;i++)
        if(x[i]<mins[i-1])
          {
	    mins[i]=x[i];
	    minis[i]=i;
          }
        else
          {
	    mins[i]=mins[i-1];
	    minis[i]=minis[i-1];
          }
      maxs[right]=x[right];
      maxis[right]=right;
      for(i=right-1;i>=left;i--)
        if(x[i]>maxs[i+1])
          {
	    maxs[i]=x[i];
	    maxis[i]=i;
          }
        else
          {
	    maxs[i]=maxs[i+1];
	    maxis[i]=maxis[i+1];
          }
      max=-1e+38;
      for(i=left;i<=right;i++)
        if(maxs[i]-mins[i]>max)
          {
	    max=maxs[i]-mins[i];
            maxi=i;
          }
      *erga=minis[maxi];
      *ergb=maxis[maxi];
      *ergmax=std::fabs(x[*erga]-x[left])+std::fabs(x[right]-x[*ergb])-std::fabs(x[right]-x[left]);
    }
}

std::vector<double> kkuiper(std::vector<double>& x, //vector with difference between unif ecdf and taut-string approx
			    long n, //length of x
			    int kmax) //the max kuiper number
{
  std::vector<double> norm(kmax,0.0);
  std::vector<int> a(kmax,0);
  std::vector<int> b(kmax,0);
  double min,max,abmax,newmax;
  //long i,k,mini,maxi,maxa,maxb,newa,newb,einordnen ;
  long i,mini,maxi,maxa,maxb,newa,newb,einordnen ;
  min=1e+38;
  max=-1e+38;
  for(i=0; i < n; i++) {
    if(x[i] > max) {
      max=x[i];
      maxi=i;
    }
    if(x[i] < min){
      min=x[i];
      mini=i;
    }
  }
  
  if(mini < maxi) {
    a[0]=mini;
    b[0]=maxi;
  }

  else {
    b[0]=mini;
    a[0]=maxi;
  }
  norm[0]=std::fabs(x[b[0]]-x[a[0]]);

  for(auto k=1; k < kmax; k++) {
    abmax=-1e+38;
    einordnen =-1;
    if(a[0]>0){
      easymax(x,n,0,a[0],&maxa,&maxb,&abmax);
      einordnen =0;
    }
    for(i=0;i<k;i++) {
      difficultmax(x,n,a[i],b[i],&newa,&newb,&newmax);
      if(newmax>abmax) {
	maxa=newa;
	maxb=newb;
	abmax=newmax;
	einordnen =i;
      }
      if(i<k-1) {
	easymax(x,n,b[i],a[i+1],&newa,&newb,&newmax); 
	if(newmax>abmax) {
	  maxa=newa;
	  maxb=newb;
	  abmax=newmax;
	  einordnen =i+1;
	}
      }
      else {
	easymax(x,n,b[i],(n-1),&newa,&newb,&newmax);
	if(newmax>abmax) {
	  maxa=newa;
	  maxb=newb;
	  abmax=newmax;
	  einordnen =i+1;
	}
      }
    }
    if(einordnen <0) {
	break;
    }
    else {
      if(einordnen == k) {
	a[k]=maxa;
	b[k]=maxb;
	norm[k]=norm[k-1]+std::fabs(x[maxb]-x[maxa]);
      }
      else {
	for(i=k;i>=einordnen+1;i--) {
	  a[i]=a[i-1];
	  b[i]=b[i-1];
	}
	if(maxa<a[einordnen]) {
	  a[einordnen]=maxa;
	  b[einordnen]=maxb;
	  norm[k]=norm[k-1]+std::fabs(x[maxb]-x[maxa]);
	}
	else {
	  b[einordnen]=maxa;
	  a[einordnen+1]=maxb;
	  norm[k]=norm[k-1]-std::fabs(x[a[einordnen]]-x[b[einordnen+1]])+std::fabs(x[a[einordnen]]-x[b[einordnen]])+std::fabs(x[a[einordnen+1]]-x[b[einordnen+1]]); 
	}
      }
    }
  }
  return(norm);
}
