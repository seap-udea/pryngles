*----------------------------------------------------------------------------
*     Main program
*----------------------------------------------------------------------------

      PROGRAM test

*======================================================================
*     Test spline
*======================================================================
c$$$      PARAMETER (n=10)
c$$$      DOUBLE PRECISION x(n),y(n),y2(n)
c$$$
c$$$      write(*,*) "Test spline"
c$$$      DO i=1,10
c$$$         x(i)=i
c$$$         y(i)=x(i)**2
c$$$      ENDDO
c$$$
c$$$      CALL spline(x,y,n,y2)
c$$$      
c$$$      DO i=1,10
c$$$         write(*,*) i,x(i),y(i),y2(i)
c$$$      ENDDO
      
*======================================================================
*     Test reflection
*======================================================================
      PARAMETER (NPIX_=1,NMUGS_=21,NMAT_=3,NFOU_=3)

      DOUBLE PRECISION phi(NPIX_),beta(NPIX_),
     .     theta0(NPIX_),theta(NPIX_),apix(NPIX_),
     .     Sarr(NPIX_,NMAT_+1)

      CHARACTER (len=100) :: filename

      DOUBLE PRECISION xmu(NMUGS_),rfou(NMAT_*NMUGS_,NMUGS_,0:NFOU_)
      
      INTEGER :: npix=NPIX_,nmat=NMAT_,nmugs=NMUGS_,nfou=NFOU_

      filename='fou_gasplanet.dat'

      CALL rdfous_planet(filename,nfou,nmat,nmugs,xmu,rfou)

      phi(1)=1.0448451569439827
      beta(1)=3.069394277348945
      theta0(1)=0.04990329026929557
      theta(1)=0.02509670973070432
      apix(1)=9.432328787795567e-05
      
      CALL reflection(
     .     npix,phi,beta,theta0,theta,
     .     nmugs,nmat,nfou,xmu,rfou,
     .     apix,Sarr)

      DO i=1,1
         DO j=1,NMAT_+1
            write(*,*) Sarr(i,j)
         ENDDO
      ENDDO

* Expected: 4.59785742e-07, 2.25122997e-07, 1.83440080e-09, 4.89642131e-01
      
      END
