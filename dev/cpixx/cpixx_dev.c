#include <Python.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

//////////////////////////////////////////////////////////////
// ROUTINES
//////////////////////////////////////////////////////////////
void test_cpixx()
{
  printf("This is Pryngles!\n");
}

double sum_matrix(double **M,int n,int m)
{
  int i,j;
  double suma=0;
  for(i=0;i<n;i++)
    for(j=0;j<m;j++)
      suma+=M[i][j];
  return suma;
}


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

  //fprintf(stderr,"klo = %d, khi = %d, h = %lf\n",klo,khi,h);

  if(fabs(h)<1e-10)
    fprintf(stderr,"Error in Splint: bad xa input.\n");

  a=(xa[khi]-x)/h;
  b=(x-xa[klo])/h;

  //printf("x = %.17lf, xa = %.17lf, xb = %.17lf, a = %.17lf, b = %.17lf\n",x,xa[klo],xa[khi],a,b);
  
  y=a*ya[klo]+b*ya[khi]+((a*a*a-a)*y2a[klo]+(b*b*b-b)*y2a[khi])*(h*h)/6;

  return y;
}

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

#define MAX_MAT 3
#define MAX_MUS 22
#define MAX_FOU 350
#define MAX_STRING 1000
#define MAX_PIX 1000
#define VERBOSITY 0
int read_fourier(char filename[],
		 long int shape[],
		 double xmu[],
		 double rfou[][MAX_MUS][MAX_FOU],
		 double rtra[][MAX_MUS][MAX_FOU]
		 )
{
  int i,j,k,m,n;
  
  char line[MAX_STRING];
  FILE *f=fopen(filename,"r");
  int nline=0;

  int nmat,nmugs,nfou;
  
  //Read first part of the file
  double fmu;

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
    if(VERBOSITY) printf("%d %.16e %.16e\n",i-1,xmu[i-1],fmu);
    if(i==nmugs)
      break;
  }

  //Initialize the fourier cube
  for(i=0;i<nmat*nmugs;i++)
    for(j=0;j<nmugs;j++)
      for(m=0;m<=MAX_FOU;m++){
	rfou[i][j][m]=0.0;
	rtra[i][j][m]=0.0;
      }

  //Read fourier coefficients
  int ifou,ibase;
  char *token;
  double coefficient;

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

	if(VERBOSITY) printf("Line: %s",line);

	//Read the lines
	token=strtok(line," ");
	n=0;
	k=0;
	//printf("\nToken: %s\n",token);
	do{
	  if(n<3){
	    sscanf(token,"%d",&m);
	    ////if(n==2)
	    //if(VERBOSITY) printf("ifou = %d, i = %d, nmu = %d\n",ifou,i,m);
	  }else{
	    sscanf(token,"%lf",&coefficient);
	    //if(VERBOSITY) printf("k = %d, coefficient = %e\n",k,coefficient);
	    if(k<nmat){
	      rfou[ibase+k][j][ifou]=coefficient;
	      if(VERBOSITY) printf("\trfou[%d][%d][%d] = %e\n",ibase+k,j,ifou,rfou[ibase+k][j][ifou]);
	    }else{
	      rtra[ibase+k-nmat][j][ifou]=coefficient;
	      if(VERBOSITY) printf("\trtra[%d][%d][%d] = %e\n",ibase+k,j,ifou,rtra[ibase+k-nmat][j][ifou]);
	    }
	    k++;
	  }
	  token=strtok(NULL," ");
	  //printf("%s\n",token);
	  n++;
	}while(token!=NULL);
	/*
	  sscanf(line,"%d %d %d",&m,&i1,&i2);
	  sscanf(line,"%lf",&fmu);
	  printf("%d %d %d\n",m,i1,i2);
	*/
	//printf("%lf\n",fmu);
	//return 0;
	//break;
      }
      //if(i>=0) break;
    }
    //if(ifou>=0) break;
    ifou++;
  }
  nfou=ifou-1;
  
  //Read size of matrix
  
  //Read coefficients
  //printf("%s",line);
  //printf("%lu\n",strlen(line));
  if(VERBOSITY){
    printf("Number of lines : %d\n",nline);
    printf("nmugs = %d, nmat = %d, nfou = %d\n",nmugs,nmat,nfou);
    printf("Size of matrix = %d\n",nmugs*nmat*nmugs*nfou);
  }

  //Save properties of matrix
  shape[0]=nmat;
  shape[1]=nmugs;
  shape[2]=nfou;
  printf("Shape: %ld,%ld,%ld\n",shape[0],shape[1],shape[2]);
  
  //Checksum
  double checksum=0;
  for(i=0;i<nmugs*nmat;i++){
    for(j=0;j<nmugs;j++){
      for(m=0;m<nfou;m++){
	checksum+=rfou[i][j][m];
      }
    }
  }
  if(VERBOSITY||1) printf("[C] Checksum rfou: %.16e\n",checksum);

  checksum=0;
  for(i=0;i<nmugs*nmat;i++){
    for(j=0;j<nmugs;j++){
      for(m=0;m<nfou;m++){
	checksum+=rtra[i][j][m];
      }
    }
  }
  if(VERBOSITY||1) printf("[C] Checksum rtra: %.16e\n",checksum);
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
	       long int shape[],
	       double xmu[],
	       double rfou[][MAX_MUS][MAX_FOU],
	       double apix[],
	       double Sarr[][MAX_MAT+1])
{
  long int nmat=shape[0];
  long int nmugs=shape[1];
  long int nfou=shape[2];
  printf("C shape: %ld %ld %ld\n",nmat,nmugs,nfou);
  
  int i,j,k,m,n;
  double rf[MAX_MAT][MAX_MUS][MAX_MUS],rfsec[MAX_MAT][MAX_MUS][MAX_MUS];
  double rftemp[MAX_MUS],rfmu0[MAX_MAT][MAX_MUS];
  double rfsecmu0[MAX_MUS];
  double Bplus[4],SvR[MAX_MAT],RM[MAX_PIX][MAX_MAT],rf3save[MAX_MAT];
  double slice[MAX_MUS],slicep[MAX_MUS];
  double be,rf3,SvR2,SvR3,P;
  
  double mu,mu0,muold=1,mu0old=1;

  //Initialize the storage matrix and values
  for(i=0;i<npix;i++)
    for(k=0;k<nmat;k++)
      RM[i][k]=0.0;

  int ki;
  double fac;
  
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
	  ////printf("%d %d %d %lf\n",ki,n,m,rfou[ki][n][m]);
	}
	
	//Slice rf(k,j,:)
	for(n=0;n<nmugs;n++){
	  slice[n]=rf[k][j][n];
	}
	spline(xmu,slice,nmugs,rftemp);

	for(n=0;n<nmugs;n++){
	  rfsec[k][j][n]=rftemp[n];
	}
	
	/*
	for(n=0;n<nmugs;n++){
	  printf("%d %d %d %e %lf %lf\n",k,j,n,xmu[n],slice[n],rftemp[n]);
	}
	if(k>0) return 0;
	*/
	
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
