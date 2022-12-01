#include <Python.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

//////////////////////////////////////////////////////////////
// ROUTINES
//////////////////////////////////////////////////////////////
/*
  Test cpixx
 */
void test_cpixx()
{
  printf("This is CPIXX in Pryngles!\n");
}

/*
*----------------------------------------------------------------------------
*     Spline interpolation routine from Press et al. (1986, p.88). 
*
*     Given arrays x and y of length n containing a tabulated function,
*     i.e. y(i)=f(x(i)), with x(1)<x(2)<...<x(n), and given values yp1 
*     and ypn for the first derivative of the interpolating function at
*     points 1 and n respectively, this routine returns an array y2 of  
*     length n which contains the second derivatives of the interpola- 
*     ting function at the tabulated points x(i).                     
*
*     If yp1 and/or yp2 are equal to 1x10^30 or larger, the routine is 
*     signalled to set the corresponding boundary condition for a natu- 
*     ral spline, with zero second derivative on that boundary.        
*
*     n is the number of elements in x and y
*----------------------------------------------------------------------------
*/
int spline(double x[],double y[],int n,
	   double y2[])
{
  int i,k;
  double u[1000];
  double sig,p,qn,un;
  
  y2[0] = 0;
  u[0] =  0;

  for(i=1;i<n-1;i++){
    sig= (x[i]-x[i-1])/(x[i+1]-x[i-1]);
    p= sig*y2[i-1]+2;
    y2[i]= (sig-1)/p;
    u[i]= (6*((y[i+1]-y[i])/(x[i+1]-x[i]) - (y[i]-y[i-1])/
	      (x[i]-x[i-1]))/(x[i+1]-x[i-1]) - sig*u[i-1])/p;
  }

  qn= 0;
  un= 0;
  y2[n-1]= (un-qn*u[n-2])/(qn*y2[n-2]+1);

  for(k=n-2;k>=0;k--){
    y2[k]= y2[k]*y2[k+1]+u[k];
  }

  return 0;
}
