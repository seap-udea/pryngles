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

      INTEGER i,l,nmat,nmugs,nfou,npix,ki,m,k,n
      INTEGER j1(npix)

      DOUBLE PRECISION fac,mu,mu0,be,rf3,SVR2,SvR3

      DOUBLE PRECISION xmu(nmugs),
     .                 phi(npix),
     .                 beta(npix),
     .                 theta0(npix),
     .                 theta(npix),
     .                 apix(npix),
     .                 rfou(nmat*nmugs,nmugs,0:nfou),
     .                 rfm(nmat),Bplus(4),
     .                 rfj(nmugs,nmat),rf2(nmugs),rf(nmugs),
     .                 SvR(nmat),RM(nmat),
     .                 P,
     .                 Sarr(npix,nmat+1)

Cf2py intent(in) npix, phi, beta, theta0, theta, apix, trans
Cf2py intent(in) nmugs, nmat, nfou, xmu, rfou
Cf2py intent(out) Sarr
Cf2py depend(npix) phi, beta, theta0, theta, Sarr, apix
Cf2py depend(nmat) rfou, Sarr
Cf2py depend(nfou) rfou
Cf2py depend(nmugs) rfou, xmu

*----------------------------------------------------------------------------
*     Find the locations in array xmu where the mu-values are:
*----------------------------------------------------------------------------
      CALL bracks(theta,npix,xmu,nmugs,j1)

*----------------------------------------------------------------------------
*     Loop over the pixels:
*----------------------------------------------------------------------------
      DO i=1,npix

*----------------------------------------------------------------------------
*     Initialize the reflected light vector for this pixel:
*----------------------------------------------------------------------------
        DO k=1,nmat
           RM(k)= 0.D0
           SvR(k)= 0.D0
        ENDDO

*       Get the cosines of the angles:
        mu= theta(i)
        mu0=theta0(i)

*----------------------------------------------------------------------------
*       Loop over the Fourier coefficients:
*----------------------------------------------------------------------------
        DO m=0,nfou

           fac=1.D0
           IF (m.EQ.0) fac=0.5D0

           IF (j1(i).NE.0) THEN
              DO k=1,nmat
                 ki= (j1(i)-1)*nmat+k 
                 DO n=1,nmugs
                    rf(n)= rfou(ki,n,m)
                 ENDDO
                 CALL spline(xmu,rf,nmugs,rf2)
                 CALL splint(xmu,rf,rf2,nmugs,mu0,rf3)
                 rfm(k)= rf3
              ENDDO
           ELSE
              DO n=1,nmugs
                 DO k=1,nmat
                    ki= (n-1)*nmat+k 
                    DO l=1,nmugs
                       rf(l)= rfou(ki,l,m)
                    ENDDO
                    CALL spline(xmu,rf,nmugs,rf2)
                    CALL splint(xmu,rf,rf2,nmugs,mu0,rf3)
                    rfj(n,k)= rf3
                 ENDDO
              ENDDO
              DO k=1,nmat
                 DO n=1,nmugs
                    rf(n)= rfj(n,k)
                 ENDDO
                 CALL spline(xmu,rf,nmugs,rf2)
                 CALL splint(xmu,rf,rf2,nmugs,mu,rf3)
                 rfm(k)= rf3
              ENDDO
           ENDIF

*----------------------------------------------------------------------------
*          Calculate the 1st column of the local reflection matrix: 
*----------------------------------------------------------------------------
           Bplus(1)= DCOS(m*phi(i))
           Bplus(2)= DCOS(m*phi(i))
           Bplus(3)= DSIN(m*phi(i))
           Bplus(4)= DSIN(m*phi(i))

           DO k=1,nmat
              RM(k)= RM(k) + 2.D0*Bplus(k)*fac*rfm(k)
           ENDDO

        ENDDO

*----------------------------------------------------------------------------
*       Calculate the locally reflected Stokes vector:
*       apix is the surface area of the pixel (the same for all pixels)
*       apix is initially assumed to be a piece of the illuminated circle 
*       with radius 1
*----------------------------------------------------------------------------
        DO k=1,nmat
           SvR(k)= mu0*RM(k) 
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
*        IF (DABS(SvR(2)).LT.1.D-6) SvR(2)=0.D0
*        IF (DABS(SvR(3)).LT.1.D-6) SvR(3)=0.D0
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
