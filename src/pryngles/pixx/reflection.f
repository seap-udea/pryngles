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

      INCLUDE 'max_incl'

      INTEGER i,j,nmat,nmugs,nfou,npix,ki,m,k,n

      DOUBLE PRECISION fac,mu,mu0,mu0old,muold,
     .                 be,rf3,apix,SVR2,SvR3

      DOUBLE PRECISION xmu(nmugs),
     .                 phi(npix),
     .                 beta(npix),
     .                 theta0(npix),
     .                 theta(npix),
     .                 rfou(nmat*nmugs,nmugs,0:nfou),
     .                 rf(nmat,nmugs,nmugs),rfsec(nmat,nmugs,nmugs),
     .                 rftemp(nmugs),rfmu0(nmat,nmugs),
     .                 rfsecmu0(nmugs),     
     .                 Bplus(4),SvR(nmat),RM(npix,nmat),rf3save(nmat),
     .                 P,Sarr(npix,nmat+1)

Cf2py intent(in) npix, phi, beta, theta0, theta, apix, trans
Cf2py intent(in) nmugs, nmat, nfou, xmu, rfou
Cf2py intent(out) Sarr
Cf2py depend(npix) phi, beta, theta0, theta, Sarr
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
            Sarr(i,k) = SvR(k)*mu * apix
        ENDDO
        
        Sarr(i,nmat+1) = P

*----------------------------------------------------------------------------
*       Next pixel:
*----------------------------------------------------------------------------

      ENDDO
      
*-------------------------------------------------------------------------------
      RETURN
      END
