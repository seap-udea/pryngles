#include <Python.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

//////////////////////////////////////////////////////////////
// ROUTINES
//////////////////////////////////////////////////////////////
void test_cpixx()
{
  printf("This is Pryngles!\n");
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

  double sum=0.0;
  for(i=0;i<n;i++){
    sum+=y2[i];
    printf("%lf\n",y2[i]);
  }

  printf("%lf\n",sum);
  return sum;
}

#define MAX_MAT 5
#define MAX_MUS 30
#define MAX_FOU 5
int reflection(int npix,
	       double phi,double beta,double theta0,double theta,
	       int nmugs,int nmat,int nfou,
	       double xmu[],double rfou[][MAX_FOU],double apix[],
	       double Sarr[][MAX_MAT])
{
  
  return 0;
}


/*
      SUBROUTINE reflection(npix,phi,beta,theta0,theta,
     .                      nmugs,nmat,nfou,xmu,rfou,
     .                      apix,Sarr)
      
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
      IMPLICIT NONE

      DOUBLE PRECISION pi,radfac
      PARAMETER (pi=3.141592653589793D0,radfac=pi/180.D0)

      INTEGER i,j,nmat,nmugs,nfou,npix,ki,m,k,n

      DOUBLE PRECISION fac,mu,mu0,mu0old,muold,
     .                 be,rf3,SVR2,SvR3

      DOUBLE PRECISION xmu(nmugs),
     .                 phi(npix),
     .                 beta(npix),
     .                 theta0(npix),
     .                 theta(npix),
     .                 apix(npix),
     .                 rfou(nmat*nmugs,nmugs,0:nfou),
     .                 rf(nmat,nmugs,nmugs),rfsec(nmat,nmugs,nmugs),
     .                 rftemp(nmugs),rfmu0(nmat,nmugs),
     .                 rfsecmu0(nmugs),     
     .                 Bplus(4),SvR(nmat),RM(npix,nmat),rf3save(nmat),
     .                 P,Sarr(npix,nmat+1)

Cf2py intent(in) npix, phi, beta, theta0, theta, apix, trans
Cf2py intent(in) nmugs, nmat, nfou, xmu, rfou
Cf2py intent(out) Sarr
Cf2py depend(npix) phi, beta, theta0, theta, apix, Sarr
Cf2py depend(nmat) rfou, Sarr
Cf2py depend(nfou) rfou
Cf2py depend(nmugs) rfou, xmu

*----------------------------------------------------------------------------
*     Initialize the storage matrix and values:
*----------------------------------------------------------------------------
      DO i=1,npix
        DO k=1,nmat
            RM(i,k) = 0.D0
        ENDDO
      ENDDO
      muold=1.D0
      mu0old=1.D0
*----------------------------------------------------------------------------
*     Loop over the Fourier coefficients:
*----------------------------------------------------------------------------
      DO m=0,nfou
      
        fac=1.D0
        IF (m.EQ.0) fac=0.5D0
          
*------------------------------------------------------------------------------
*     Initialize the interpolation matrix for the current fourier coefficient:
*------------------------------------------------------------------------------
        DO j=1,nmugs
            DO k=1,nmat
                ki = (j-1)*nmat + k
                DO n=1,nmugs
                    rf(k,j,n) = rfou(ki,n,m)
                ENDDO
                CALL spline(xmu,rf(k,j,:),nmugs,rftemp)
                rfsec(k,j,:) = rftemp
            ENDDO
        ENDDO
        
*----------------------------------------------------------------------------
*     Loop over the pixels:
*       If the input angles are (very) similar to a previously calculated case
*       use those values.
*       To obtain obtain the fourier coefficient at (mu,mu0) spline has to be
*       called a second time.
*----------------------------------------------------------------------------      
        DO i=1,npix
            mu = theta(i)
            mu0 = theta0(i)
            Bplus(1)= DCOS(m*phi(i))
            Bplus(2)= DCOS(m*phi(i))
            Bplus(3)= DSIN(m*phi(i))
            Bplus(4)= DSIN(m*phi(i))
            
            IF ((i.GT.1).AND.
     .          (ABS(mu-muold).LT.1.D-6).AND.
     .          (ABS(mu0-mu0old).LT.1.D-6)) THEN
                DO k=1,nmat 
                    RM(i,k)=RM(i,k)+ 2.D0*Bplus(k)*fac*rf3save(k)
                ENDDO
            ELSE
                DO j=1,nmugs
                    DO k=1,nmat
                        CALL splint(xmu,rf(k,j,:),
     .                              rfsec(k,j,:),
     .                              nmugs,mu0,rf3)
                        rfmu0(k,j) = rf3
                    ENDDO
                ENDDO                    
                DO k=1,nmat
                    CALL spline(xmu,rfmu0(k,:),nmugs,rfsecmu0)
                    CALL splint(xmu,rfmu0(k,:),rfsecmu0,nmugs,mu,rf3)
                    rf3save(k) = rf3
                    muold = mu
                    mu0old = mu0
                    RM(i,k) = RM(i,k) + 2.D0*Bplus(k)*fac*rf3
                ENDDO
            ENDIF
        ENDDO
      ENDDO
*----------------------------------------------------------------------------
*     Loop again over the pixels to rotate Stokes vector:
*----------------------------------------------------------------------------    
      DO i=1,npix
        mu = theta(i)
        mu0 = theta0(i)
        
*----------------------------------------------------------------------------
*       Calculate the locally reflected Stokes vector:
*----------------------------------------------------------------------------
        DO k=1,nmat
           SvR(k)= mu0*RM(i,k) 
        ENDDO

*----------------------------------------------------------------------------
*       Rotate Stokes elements Q and U to the actual reference plane:
*----------------------------------------------------------------------------
        be= 2.D0*beta(i)
        SvR2= DCOS(be)*SvR(2) + DSIN(be)*SvR(3) 
        SvR3=-DSIN(be)*SvR(2) + DCOS(be)*SvR(3) 
        SvR(2)= SvR2
        SvR(3)= SvR3

*----------------------------------------------------------------------------
*       Compute the local degree of polarisation P:
*----------------------------------------------------------------------------
        IF (DABS(Svr(1)).LT.1.D-6) THEN
           P=0.D0
        ELSEIF (DABS(SvR(3)).LT.1.D-6) THEN
           P= -SvR(2)/SvR(1)
        ELSE
           P= DSQRT(SvR(2)*SvR(2)+SvR(3)*SvR(3))/SvR(1)
        ENDIF
        IF (DABS(P).LT.1.D-6) P=0.D0

*----------------------------------------------------------------------------
*       Add the Stokes elements of the pixel to an array:
*       Multiply with mu and the actual pixel area to obtain stokes elements
*----------------------------------------------------------------------------
        DO k=1,nmat
            Sarr(i,k) = SvR(k)*mu * apix(i)
        ENDDO
        
        Sarr(i,nmat+1) = P

*----------------------------------------------------------------------------
*       Next pixel:
*----------------------------------------------------------------------------

      ENDDO
      
*-------------------------------------------------------------------------------
      RETURN
      END
 */


/*
      SUBROUTINE bracks(theta,npix,xmu,nmugs,j1)

*-------------------------------------------------------------------------------
* PURPOSE:
* Find a bracket around "el" in the first "n" elements of
* "array".
* 
* INPUT:
* 	el	: element to be bracketed
* 	array   : array in which bracket must be found
* 	n       : number of elements in array to be considered
* 	ND      : dimension of array
*
* OUTPUT:
*	i1	: index of left bracket element
* 	i2	: index of right bracket element
*
*
* COMMENTS:
* If "el" is outside the range of "array", "i1" and "i2" are 
* set to the appropriate extreme index value.    
* If "el" is equal to one of the elements of "array", "i2" is
* set to the appropriate index value.
*-------------------------------------------------------------------------------
      IMPLICIT DOUBLE PRECISION (a-h,o-z)

      DOUBLE PRECISION eps
      PARAMETER (eps=1.D-10)

      DOUBLE PRECISION pi,radfac
      PARAMETER (pi=3.141592653589793D0,radfac=pi/180.D0)

      INTEGER ni,i,nmugs,npix,i1

      INTEGER j1(npix)

      DOUBLE PRECISION mu,xmu(nmugs),theta(npix)
      
Cf2py intent(in) theta, npix, xmu, nmugs
Cf2py intent(out) j1
Cf2py depend(npix) theta, j1
Cf2py depend(nmugs) xmu
*-------------------------------------------------------------------------------

      DO ni=1,npix
        i1= -1
        j1(ni)= i1
        
        mu= theta(ni)
        
        IF (DACOS(mu)*radfac.LT.(0.D0)) GOTO 33
        
        IF (mu.LE.(xmu(1)+eps)) THEN
           i1= 0
        ELSEIF (mu.GT.(xmu(nmugs)-eps)) THEN
           i1= nmugs
        ELSE
           DO i=1,nmugs-2
               IF ((mu.GT.(xmu(i)+eps)).AND.(mu.LE.(xmu(i+1)+eps))) THEN
                  i1= i
               ENDIF
           ENDDO
           IF ((mu.GT.(xmu(nmugs-1)+eps)).AND.
     .           (mu.LE.(xmu(nmugs)-eps))) THEN
              i1= nmugs-1
           ENDIF
        ENDIF
        j1(ni)= i1

33      CONTINUE

      ENDDO

*-------------------------------------------------------------------------------
      RETURN
*/
