*----------------------------------------------------------------------------
*     Main program
*----------------------------------------------------------------------------

      PROGRAM test

*======================================================================
*     Test bracks
*======================================================================
c$$$      DOUBLE PRECISION xmu(20)
c$$$      DOUBLE PRECISION x
c$$$      
c$$$      OPEN(unit=1,file='xmu.mat',status='old',err=999)
c$$$
c$$$      DO i=1,20
c$$$         READ(1,*) xmu(i)
c$$$         write (*,*) xmu(i)
c$$$      ENDDO
c$$$
c$$$      DO i=1,10
c$$$         x=i/10
c$$$      ENDDO
c$$$      
c$$$      STOP
c$$$999   WRITE(*,*) 'Error opening file!'
      
*======================================================================
*     Test spline
*======================================================================
c$$$      PARAMETER (n=10)
c$$$      DOUBLE PRECISION x(n),y(n),y2(n)
c$$$      DOUBLE PRECISION xv,yv
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
c$$$
c$$$      xv = 4.3
c$$$
c$$$      DO i=1,9
c$$$         xv = i+3.D-1
c$$$         CALL splint(x,y,y2,n,xv,yv)
c$$$         write(*,*) xv,yv
c$$$      ENDDO
      
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
c$$$c$$$      DO i=1,nmugs*nmat
c$$$c$$$         DO j=1,nmugs
c$$$c$$$            DO m=0,nfou-1
c$$$c$$$               write(*,*) i,j,m,rfou(i,j,m)
c$$$c$$$            ENDDO
c$$$c$$$         ENDDO
c$$$c$$$      ENDDO
c$$$c$$$      STOP
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
c$$$      WRITE(*,*) "What we obtain:"
c$$$      DO i=1,1
c$$$         DO j=1,NMAT_+1
c$$$            write(*,*) Sarr(i,j)
c$$$         ENDDO
c$$$      ENDDO
c$$$
c$$$      WRITE(*,*) "What is expected:"
c$$$      WRITE(*,*) 4.59785742e-07,2.25122997e-07,1.83440080e-09,
c$$$     .          4.89642131e-01
c$$$      
c$$$c$$$* Expected: 4.59785742e-07, 2.25122997e-07, 1.83440080e-09, 4.89642131e-01

*======================================================================
*     Test reflection ring
*======================================================================
      PARAMETER (NPIX_=1,NMUGS_=20,NMAT_=3,NFOU_=344)

      DOUBLE PRECISION phi(NPIX_),beta(NPIX_),
     .     theta0(NPIX_),theta(NPIX_),apix(NPIX_),
     .     Sarr(NPIX_,NMAT_+1)

      CHARACTER (len=100) :: filename

      DOUBLE PRECISION xmu(NMUGS_),rfou(NMAT_*NMUGS_,NMUGS_,0:NFOU_)
      
      INTEGER :: npix=NPIX_,nmat=NMAT_,nmugs=NMUGS_,nfou=NFOU_

      filename='fou_ring_0_4_0_8.dat'
      
c$$$      CALL rdfous_ring(filename,.FALSE.,nfou,nmat,nmugs,xmu,rfou)
c$$$      DO i=1,nmugs*nmat
c$$$         DO j=1,nmugs
c$$$            DO m=0,nfou-1
c$$$               write(*,*) i,j,m,rfou(i,j,m)
c$$$            ENDDO
c$$$         ENDDO
c$$$      ENDDO
c$$$      STOP
      
      phi(1)=1.490116119384765625e-08
      beta(1)=0.000000000000000000e+00
      theta0(1)=4.999999999999998890e-01
      theta(1)=5.000000000000000000e-01
      apix(1)=1.163314390931110409e-04
      
      CALL reflection(
     .     npix,phi,beta,theta0,theta,
     .     nmugs,nmat,nfou,xmu,rfou,
     .     apix,Sarr)

      WRITE(*,*) "What we obtain:"
      DO i=1,1
         DO j=1,NMAT_+1
            write(*,*) Sarr(i,j)
         ENDDO
      ENDDO

      WRITE(*,*) "What is expected:"
      WRITE(*,*) 6.193646058775441171e-06,-3.046132218287162406e-07,
     .          -9.759550180589223642e-15,4.918156751904256135e-02
      
c$$$* Expected: 4.59785742e-07, 2.25122997e-07, 1.83440080e-09, 4.89642131e-01

*======================================================================
*     Test transmission ring
*======================================================================
c$$$      PARAMETER (NPIX_=1,NMUGS_=20,NMAT_=3,NFOU_=344)
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
c$$$      filename='fou_ring_0_4_0_8.dat'
c$$$      
c$$$      CALL rdfous_ring(filename,.TRUE.,nfou,nmat,nmugs,xmu,rfou)
c$$$
c$$$c$$$      DO i=1,nmugs*nmat
c$$$c$$$         DO j=1,nmugs
c$$$c$$$            DO m=0,nfou-1
c$$$c$$$               write(*,*) i,j,m,rfou(i,j,m)
c$$$c$$$            ENDDO
c$$$c$$$         ENDDO
c$$$c$$$      ENDDO
c$$$c$$$      STOP
c$$$      
c$$$      phi(1)=1.601029385538801364e+00
c$$$      beta(1)=1.601029385538801364e+00
c$$$      theta0(1)=1.744974835125044643e-02
c$$$      theta(1)=5.000000000000000000e-01
c$$$      apix(1)=1.163314390931110409e-04
c$$$      
c$$$      CALL reflection(
c$$$     .     npix,phi,beta,theta0,theta,
c$$$     .     nmugs,nmat,nfou,xmu,rfou,
c$$$     .     apix,Sarr)
c$$$
c$$$      WRITE(*,*) "What we obtain:"
c$$$      DO i=1,1
c$$$         DO j=1,NMAT_+1
c$$$            write(*,*) Sarr(i,j)
c$$$         ENDDO
c$$$      ENDDO
c$$$
c$$$      WRITE(*,*) "What is expected:"
c$$$      WRITE(*,*) 1.688771016436214060e-07,-1.485273114913191718e-08,
c$$$     .          2.674513729889046925e-09,8.936444506313186154e-02
c$$$
      END
