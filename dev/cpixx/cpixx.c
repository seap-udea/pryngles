#include <Python.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

//////////////////////////////////////////////////////////////
//MACROS
//////////////////////////////////////////////////////////////
#define MAX_MAT 3
#define MAX_MUS 22
#define MAX_FOU 350
#define MAX_STRING 1000
#define MAX_PIX 1000

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

/*
*----------------------------------------------------------------------------
*     Spline interpolation routine from Press et al. (1986, p.88).   
*                                                                 
*     Given the arrays xa and ya of length n, which tabulate a function
*     (with the xa(i)'s in order), and given the array y2a, which is  
*     the output from SPLINE above, and given a value of x, this     
*     routine returns a cubic-spline interpolated value y.       
*----------------------------------------------------------------------------
*/
double splint(double xa[],double ya[],double y2a[],int n,double x)
{
  int k;
  int klo=0;
  int khi=n-1;
  double h,a,b,y;
  
  while((khi-klo)>1){
    k=floor((klo+khi)/2);
    if(xa[k]>x)
      khi=k;
    else
      klo=k;
  }
  
  h=xa[khi]-xa[klo];

  if(fabs(h)<1e-10)
    fprintf(stderr,"Error in Splint: bad xa input.\n");

  a=(xa[khi]-x)/h;
  b=(x-xa[klo])/h;

  y=a*ya[klo]+b*ya[khi]+((a*a*a-a)*y2a[klo]+(b*b*b-b)*y2a[khi])*(h*h)/6;

  return y;
}

/*
*----------------------------------------------------------------------------
Read Fourier coefficients
*----------------------------------------------------------------------------
*/
int read_fourier(char filename[],
		 int shape[],
		 double xmu[],
		 double rfou[][MAX_MUS][MAX_FOU],
		 double rtra[][MAX_MUS][MAX_FOU]
		 )
{
  FILE *f=fopen(filename,"r");
  char line[MAX_STRING];
  int i,j,k,m,n;
  int nline=0;
  int nmat,nmugs,nfou;
  double fmu;
  int ifou,ibase;
  char *token;
  double coefficient;

  i=0;
  while(fgets(line,sizeof(line),f)) {
    nline++;
    //Avoid comments
    if(strchr(line,'#')!=NULL) continue;

    //Read size of matrix
    if(strlen(line)<10){
      if(!nmat)
	sscanf(line,"%d\n",&nmat);
      else
	sscanf(line,"%d\n",&nmugs);
      continue;
    }

    //Read mus
    sscanf(line,"%lf %lf\n",&xmu[i++],&fmu);
    if(i==nmugs) break;
  }

  //Initialize the fourier cube
  for(i=0;i<nmat*nmugs;i++)
    for(j=0;j<nmugs;j++)
      for(m=0;m<=MAX_FOU;m++){
	rfou[i][j][m]=0.0;
	rtra[i][j][m]=0.0;
      }

  //Read fourier coefficients
  ifou=0;
  while(1){
    if(feof(f)) break;
    for(i=0;i<nmugs;i++){
      ibase=i*nmat;
      for(j=0;j<nmugs;j++){

	//Stop if we get the end of the file
	if(fgets(line,sizeof(line),f)==NULL)
	  break;
	nline++;

	//Read the lines
	token=strtok(line," ");
	n=0;
	k=0;

	do{
	  if(n<3){
	    sscanf(token,"%d",&m);
	  }else{
	    sscanf(token,"%lf",&coefficient);
	    if(k<nmat){
	      rfou[ibase+k][j][ifou]=coefficient;
	    }else{
	      rtra[ibase+k-nmat][j][ifou]=coefficient;
	    }
	    k++;
	  }
	  token=strtok(NULL," ");
	  n++;
	}while(token!=NULL);
      }
    }
    ifou++;
  }
  nfou=ifou-1;
  
  //Store properties of matrix
  shape[0]=nmat;
  shape[1]=nmugs;
  shape[2]=nfou;

  fclose(f);
  return 0;
}

/*
*----------------------------------------------------------------------------
*     Read a Fourier coefficients file and calculate the Stokes vector
*     for a given geometry.
*
*     It is assumed that the Stokes vector of the incoming sunlight 
*     is [1,0,0,0], with the flux measured perpendicular to the direction
*     of incidence equal to pi.
*     
*     phi and beta are assumed to be in radian
*     theta0 and theta are given as cos(theta0) and cos(theta) respectively
*
*     Author: Daphne M. Stam
*     Date: September 2022
*----------------------------------------------------------------------------
 */
int reflection(int npix,
	       double phi[],double beta[],double theta0[],double theta[],
	       int shape[],
	       double xmu[],
	       double rfou[][MAX_MUS][MAX_FOU],
	       double apix[],
	       double Sarr[][MAX_MAT+1])
{
  int i,j,k,m,n;
  double rf[MAX_MAT][MAX_MUS][MAX_MUS],rfsec[MAX_MAT][MAX_MUS][MAX_MUS];
  double rftemp[MAX_MUS],rfmu0[MAX_MAT][MAX_MUS];
  double rfsecmu0[MAX_MUS];
  double Bplus[4],SvR[MAX_MAT],RM[MAX_PIX][MAX_MAT],rf3save[MAX_MAT];
  double slice[MAX_MUS],slicep[MAX_MUS];
  double be,rf3,SvR2,SvR3,P;
  double mu,mu0,muold=1,mu0old=1;
  int ki;
  double fac;

  //Sizes
  int nmat=shape[0];
  int nmugs=shape[1];
  int nfou=shape[2];

  //Initialize the storage matrix and values
  for(i=0;i<npix;i++)
    for(k=0;k<nmat;k++)
      RM[i][k]=0.0;

  //Loop over the Fourier coefficients:
  for(m=0;m<nfou;m++){

    fac=1.0;
    if(m==0) fac=0.5;

    //Initialize the interpolation matrix for the current fourier coefficient:
    for(j=0;j<nmugs;j++){

      for(k=0;k<nmat;k++){
	ki=j*nmat+k;
	
	for(n=0;n<nmugs;n++){
	  rf[k][j][n]=rfou[ki][n][m];
	}
	
	//Slice rf(k,j,:)
	for(n=0;n<nmugs;n++){
	  slice[n]=rf[k][j][n];
	}
	spline(xmu,slice,nmugs,rftemp);

	for(n=0;n<nmugs;n++){
	  rfsec[k][j][n]=rftemp[n];
	}
	
      }//End k

    }//End j
    
    /*
     *----------------------------------------------------------------------------
     *     Loop over the pixels:
     *       If the input angles are (very) similar to a previously calculated case
     *       use those values.
     *       To obtain obtain the fourier coefficient at (mu,mu0) spline has to be
     *       called a second time.
     *----------------------------------------------------------------------------
     */
    for(i=0;i<npix;i++){
      mu = theta[i];
      mu0 = theta0[i];
      Bplus[0]= cos(m*phi[i]);
      Bplus[1]= cos(m*phi[i]);
      Bplus[2]= sin(m*phi[i]);
      Bplus[3]= sin(m*phi[i]);

      if((i>0)&&((fabs(mu-muold)<1e-6))&&(fabs(mu0-mu0old)<1e-6)){
	for(k=0;k<nmat;k++)
	  RM[i][k]=RM[i][k]+2*Bplus[k]*fac*rf3save[k];
      }else{
	for(j=0;j<nmugs;j++){
	  for(k=0;k<nmat;k++){
	    //Slice rf(k,j,:)
	    for(n=0;n<nmugs;n++){
	      slice[n]=rf[k][j][n];
	    }
	    //Slice rfsec(k,j,:)
	    for(n=0;n<nmugs;n++){
	      slicep[n]=rfsec[k][j][n];
	    }
	    rf3=splint(xmu,slice,slicep,nmugs,mu0);
	    rfmu0[k][j]=rf3;
	  }
	}

	for(k=0;k<nmat;k++){
	  //Slice rfmu0(k,:)
	  for(n=0;n<nmugs;n++){
	    slice[n]=rfmu0[k][n];
	  }
	  spline(xmu,slice,nmugs,rfsecmu0);
	  rf3=splint(xmu,slice,rfsecmu0,nmugs,mu);
	  rf3save[k] = rf3;
	  muold = mu;
	  mu0old = mu0;
	  RM[i][k] = RM[i][k] + 2*Bplus[k]*fac*rf3;
	}
      }//End else
    }//End i (pix)	
  }//End loop fourier coefficients

  //Loop again over the pixels to rotate Stokes vector:
  for(i=0;i<npix;i++){
    mu = theta[i];
    mu0 = theta0[i];

    //Calculate the locally reflected Stokes vector:
    for(k=0;k<nmat;k++)
      SvR[k]= mu0*RM[i][k];
    
    //Rotate Stokes elements Q and U to the actual reference plane:
    be= 2*beta[i];
    SvR2= cos(be)*SvR[1] + sin(be)*SvR[2];
    SvR3=-sin(be)*SvR[1] + cos(be)*SvR[2];
    SvR[1]= SvR2;
    SvR[2]= SvR3;
    
    //Compute the local degree of polarisation P:
    if(fabs(SvR[0])<1e-6)
      P=0;
    else if(fabs(SvR[2])<1e-6)
      P=-SvR[1]/SvR[0];
    else
      P= sqrt(SvR[1]*SvR[1]+SvR[2]*SvR[2])/SvR[0];
    if(fabs(P)<1e-6) P=0;
    
    /*
     *----------------------------------------------------------------------------
     *       Add the Stokes elements of the pixel to an array:
     *       Multiply with mu and the actual pixel area to obtain stokes elements
     *----------------------------------------------------------------------------
     */
    for(k=0;k<nmat;k++){
      Sarr[i][k] = SvR[k]*mu*apix[i];
    }

    //The value of the degree of polarization
    Sarr[i][nmat] = P;

  }//End i (pix)
  
  return 0;
}
