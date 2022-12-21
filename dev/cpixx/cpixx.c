//////////////////////////////////////////////////////////////
// DEPENDENCIES
//////////////////////////////////////////////////////////////
#include <Python.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

//////////////////////////////////////////////////////////////
// MACROS
//////////////////////////////////////////////////////////////
#define MAX_FOU 1000
#define MAX_STRING 1000
#define VERBOSITY 0

//////////////////////////////////////////////////////////////
// ROUTINES
//////////////////////////////////////////////////////////////
struct FourierCoefficients{
  int nmat,nmugs,nfou;
  double *xmu;
  double ***rfou;
  double ***rtra;
};

//////////////////////////////////////////////////////////////
// TEST ROUTINES
//////////////////////////////////////////////////////////////
double sum_structure(struct FourierCoefficients F,int n,int m,int p)
{
  double sum=0;
  for(int i=0;i<n;i++){
    for(int j=0;j<m;j++){
      for(int k=0;k<p;k++){
	sum+=F.rfou[i][j][k];
	printf("%lf ",F.rfou[i][j][k]);
	F.rfou[i][j][k]*=2;
      }
      printf("\n");
    }
    printf("\n");
  }
  return sum;
}

double sum_cube(double ***C,int n,int m,int p)
{
  double sum=0;
  printf("%d\n",n);
  for(int i=0;i<n;i++){
    for(int j=0;j<m;j++){
      for(int k=0;k<p;k++){
	sum+=C[i][j][k];
	printf("%lf\n",C[i][j][k]);
	C[i][j][k]*=2;
      }
    }
  }
  return sum;
}

double sum_matrix(double **M,int n,int m)
{
  double sum=0;
  printf("%d\n",n);
  for(int i=0;i<n;i++){
    for(int j=0;j<m;j++){
      sum+=M[i][j];
      printf("%lf\n",M[i][j]);
      M[i][j]*=2;
    }
  }
  return sum;
}

double sum_vector(double *v,int n)
{
  double sum=0;
  printf("%d\n",n);
  for(int i=0;i<n;i++){
    sum+=v[i];
    v[i]*=2;
  }
  return sum;
}

int reflection_test(struct FourierCoefficients F,int qreflection,int npix)
{
  int i,j,m;

  printf("qreflection = %d, npix = %d\n",qreflection,npix);

  //Check sum
  printf("Size = %d\n",F.nmugs);

  //Check xmu
  for(i=0;i<F.nmugs;i++){
    printf("xmu [%d] = %.16e\n",i,F.xmu[i]);
  }
  printf("%.16e\n",F.rfou[0][0][0]);
  exit(0);

  double checksum=0;
  for(i=0;i<F.nmugs*F.nmat;i++){
    for(j=0;j<F.nmugs;j++){
      for(m=0;m<F.nfou;m++){
	checksum+=F.rfou[i][j][m]+F.rtra[i][j][m];
      }
    }
  }
  printf("Checksum = %.16e\n",checksum);
  return 0;
}

//////////////////////////////////////////////////////////////
// ROUTINES
//////////////////////////////////////////////////////////////
double* zeros_vector(int n)
{
  double *v;
  v=(double*)calloc(n,sizeof(double));
  for(int i=0;i<n;i++) v[i]=0.0;
  return v;
}

double** zeros_matrix(int n,int m)
{
  double **M;
  M=(double**)calloc(n,sizeof(double*));
  for(int i=0;i<n;i++){
    M[i]=zeros_vector(m);
    for(int j=0;j<m;j++) M[i][j]=0.0;
  }
  return M;
}

double*** zeros_cube(int n,int m,int p)
{
  double ***C;
  C=(double***)calloc(n,sizeof(double**));
  for(int i=0;i<n;i++)
    C[i]=zeros_matrix(m,p);
  return C;
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
double spline(double x[],double y[],int n,double y2[])
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
  double sum=0.0;
  for(i=0;i<n;i++){
    sum+=y2[i];
    printf("%lf\n",y2[i]);
  }

  printf("%lf\n",sum);
  return sum;
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
int reflection(struct FourierCoefficients F,int qreflection,
	       int npix,
	       double* phi,double* beta,double* theta0,double* theta,double *apix,
	       double **Sarr)
{
  //Declarations
  int i,j,k,m,n;
  double be,rf3,SvR2,SvR3,P;
  double mu,mu0,muold=1,mu0old=1;
  int ki;
  double fac;
  
  //Read sizes
  int nmat=F.nmat;
  int nmugs=F.nmugs;
  int nfou=F.nfou;

  double ***rf,***rfsec;
  double **RM,**rfmu0;
  double *rftemp,*slice,*slicep,*rfsecmu0;
  double *SvR,*rf3save;
  double Bplus[4];
  
  //Allocate dynamically temporal matrices
  //Cubes
  rf=zeros_cube(nmat,nmugs,nmugs);
  rfsec=zeros_cube(nmat,nmugs,nmugs);
  //Matrices
  rfmu0=zeros_matrix(nmat,nmugs);
  RM=zeros_matrix(npix,nmat);
  //Vectos nmugs
  rftemp=zeros_vector(nmugs);
  rfsecmu0=zeros_vector(nmugs);
  slice=zeros_vector(nmugs);
  slicep=zeros_vector(nmugs);
  //Vectos nmats
  rf3save=zeros_vector(nmat);
  SvR=zeros_vector(nmat);

  //Loop over the Fourier coefficients:
  for(m=0;m<nfou;m++){

    fac=1.0;
    if(m==0) fac=0.5;

    //Initialize the interpolation matrix for the current fourier coefficient:
    for(j=0;j<nmugs;j++){

      for(k=0;k<nmat;k++){
	ki=j*nmat+k;
	
	for(n=0;n<nmugs;n++){
	  rf[k][j][n]=qreflection?F.rfou[ki][n][m]:F.rtra[ki][n][m];
	}
	
	//Slice rf(k,j,:)
	for(n=0;n<nmugs;n++)
	  slice[n]=rf[k][j][n];
	spline(F.xmu,slice,nmugs,rftemp);

	for(n=0;n<nmugs;n++)
	  rfsec[k][j][n]=rftemp[n];
	
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
	    for(n=0;n<nmugs;n++)
	      slice[n]=rf[k][j][n];
	    
	    //Slice rfsec(k,j,:)
	    for(n=0;n<nmugs;n++)
	      slicep[n]=rfsec[k][j][n];

	    rf3=splint(F.xmu,slice,slicep,nmugs,mu0);
	    rfmu0[k][j]=rf3;
	  }
	}

	for(k=0;k<nmat;k++){
	  //Slice rfmu0(k,:)
	  for(n=0;n<nmugs;n++)
	    slice[n]=rfmu0[k][n];

	  spline(F.xmu,slice,nmugs,rfsecmu0);
	  rf3=splint(F.xmu,slice,rfsecmu0,nmugs,mu);
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

  //Free temporal arrays
  free(rf);
  free(rfsec);
  //Matrices
  free(rfmu0);
  free(RM);
  //Vectos nmugs
  free(rftemp);
  free(rfsecmu0);
  free(slice);
  free(slicep);
  //Vectos nmats
  free(rf3save);
  free(SvR);
  
  return 0;
}

/*
*----------------------------------------------------------------------------
Read Fourier coefficients
*----------------------------------------------------------------------------
*/
double read_fourier(char filename[],struct FourierCoefficients* F)
{
  int i,j,k,m,n;
  char line[MAX_STRING];
  int nline=0;
  int nmat,nmugs,nfou;
  double fmu;
  int ifou,ibase;
  char *token;
  double coefficient;
  double checksum;

  //Open file
  FILE *f=fopen(filename,"r");
  
  //Read first part of the file
  i=0;

  while(fgets(line,sizeof(line),f)) {
    nline++;
    //Avoid comments
    if(strchr(line,'#')!=NULL) continue;

    //Read size of matrix
    if(strlen(line)<10){
      if(!nmat){
	sscanf(line,"%d\n",&nmat);
	F->nmat=nmat;
      }
      else{
	sscanf(line,"%d\n",&nmugs);
	F->nmugs=nmugs;
	F->xmu=zeros_vector(F->nmugs);
      }
      continue;
    }

    //Read mus
    sscanf(line,"%lf %lf\n",&(F->xmu[i++]),&fmu);
    if(i==nmugs)
      break;
  }
  
  //Initialize the fourier cube
  double*** rfou_read=zeros_cube(F->nmat*F->nmugs,F->nmugs,MAX_FOU);
  double*** rtra_read=zeros_cube(F->nmat*F->nmugs,F->nmugs,MAX_FOU);
  
  //Read fourier coefficients
  ifou=0;
  while(1){
    if(feof(f)) break;
    for(i=0;i<F->nmugs;i++){
      ibase=i*F->nmat;
      for(j=0;j<F->nmugs;j++){

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
	      rfou_read[ibase+k][j][ifou]=coefficient;
	    }else{
	      rtra_read[ibase+k-nmat][j][ifou]=coefficient;
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
  F->nfou=ifou-1;
  
  //Save properties of matrix
  printf("Read fourier: nmat = %d, numgs = %d, nfou = %d\n",F->nmat,F->nmugs,F->nfou);

  F->rfou=zeros_cube(F->nmugs*F->nmat,F->nmugs,F->nfou);
  F->rtra=zeros_cube(F->nmugs*F->nmat,F->nmugs,F->nfou);
  
  //Checksum
  checksum=0;
  for(i=0;i<F->nmugs*F->nmat;i++){
    for(j=0;j<F->nmugs;j++){
      for(m=0;m<F->nfou;m++){
	F->rfou[i][j][m]=rfou_read[i][j][m];
	F->rtra[i][j][m]=rtra_read[i][j][m];
	checksum+=F->rfou[i][j][m]+F->rtra[i][j][m];
      }
    }
  }

  free(rfou_read);
  free(rtra_read);
  
  return checksum;
}
