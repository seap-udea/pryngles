*----------------------------------------------------------------------------
*     Main program
*----------------------------------------------------------------------------

      PROGRAM test

*======================================================================
*     Test spline
*======================================================================
      PARAMETER (n=10)
      DOUBLE PRECISION x(n),y(n),y2(n)

      write(*,*) "Test spline"
      DO i=1,10
         x(i)=i
         y(i)=x(i)**2
      ENDDO

      CALL spline(x,y,n,y2)
      
      DO i=1,10
         write(*,*) i,x(i),y(i),y2(i)
      ENDDO
      
*======================================================================
*     Test reflection
*======================================================================
c$$$      PARAMETER (NPIX_=1,NMUGS_=21,NMAT_=3,NFOU_=3)
c$$$
c$$$      DOUBLE PRECISION phi(NPIX_),beta(NPIX_),
c$$$     .     theta0(NPIX_),theta(NPIX_),apix(NPIX_),
c$$$     .     Sarr(NPIX_,NMAT_+1)
c$$$
c$$$      CHARACTER (len=100) :: filename
c$$$
c$$$      DOUBLE PRECISION xmu(NMUGS_),rfou(NMAT_*NMUGS_,NMUGS_,0:NFOU_)
c$$$      
c$$$      INTEGER :: npix=NPIX_,nmat=NMAT_,nmugs=NMUGS_,nfou=NFOU_
c$$$
c$$$      filename='fou_gasplanet.dat'
c$$$
c$$$      CALL rdfous_planet(filename,nfou,nmat,nmugs,xmu,rfou)
c$$$
c$$$      phi(1)=1.0448451569439827
c$$$      beta(1)=3.069394277348945
c$$$      theta0(1)=0.04990329026929557
c$$$      theta(1)=0.02509670973070432
c$$$      apix(1)=9.432328787795567e-05
c$$$      
c$$$      CALL reflection(
c$$$     .     npix,phi,beta,theta0,theta,
c$$$     .     nmugs,nmat,nfou,xmu,rfou,
c$$$     .     apix,Sarr)
c$$$
c$$$      DO i=1,1
c$$$         DO j=1,NMAT_+1
c$$$            write(*,*) Sarr(i,j)
c$$$         ENDDO
c$$$      ENDDO
c$$$
c$$$* Expected: 4.59785742e-07, 2.25122997e-07, 1.83440080e-09, 4.89642131e-01
      
      END
