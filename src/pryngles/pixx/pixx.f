      SUBROUTINE spline(x,y,n,y2)

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
      IMPLICIT DOUBLE PRECISION (a-h,o-z)

      INTEGER n

      DOUBLE PRECISION x(n),y(n),y2(n),u(n)

Cf2py intent(in) x, y, n
Cf2py intent(out) y2
Cf2py depend(n) x, y, y2
*----------------------------------------------------------------------------
      y2(1)= 0.D0
      u(1)=  0.D0

      DO i=2,n-1
	 sig= (x(i)-x(i-1))/(x(i+1)-x(i-1))
	 p= sig*y2(i-1)+2.D0
	 y2(i)= (sig-1.D0)/p
	 u(i)= (6.D0*((y(i+1)-y(i))/(x(i+1)-x(i))-(y(i)-y(i-1))
     .       /(x(i)-x(i-1)))/(x(i+1)-x(i-1))-sig*u(i-1))/p
      ENDDO

      qn= 0.D0
      un= 0.D0
      y2(n)= (un-qn*u(n-1))/(qn*y2(n-1)+1.D0)

      DO k=n-1,1,-1
	 y2(k)= y2(k)*y2(k+1)+u(k)
      ENDDO

*----------------------------------------------------------------------------
      RETURN
      END

*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      SUBROUTINE splint(xa,ya,y2a,n,x,y)

*----------------------------------------------------------------------------
*     Spline interpolation routine from Press et al. (1986, p.88).   
*                                                                 
*     Given the arrays xa and ya of length n, which tabulate a function
*     (with the xa(i)'s in order), and given the array y2a, which is  
*     the output from SPLINE above, and given a value of x, this     
*     routine returns a cubic-spline interpolated value y.       
*----------------------------------------------------------------------------
      IMPLICIT DOUBLE PRECISION (a-h,o-z)
      
      INTEGER n
      
      DOUBLE PRECISION x,y
      
      DOUBLE PRECISION xa(n),ya(n),y2a(n)

Cf2py intent(in) xa, ya, y2a, n, x
Cf2py intent(out) y
Cf2py depend(n) xa, ya, y2a
*----------------------------------------------------------------------------
      klo=1.0
      khi=n

    1 IF (khi-klo.GT.1) THEN
	 k= (klo+khi)/2.D0  
         IF (xa(k).GT.x) THEN
	    khi= k
         ELSE
	    klo= k
         ENDIF
         GOTO 1
      ENDIF

      h= xa(khi)-xa(klo)

      IF (DABS(h).LT.1.D-10) WRITE(*,10) 
      a= (xa(khi)-x)/h
      b= (x-xa(klo))/h
      y= a*ya(klo)+b*ya(khi)+
     +  ((a**3-a)*y2a(klo)+(b**3-b)*y2a(khi))*(h**2)/6.D0

*----------------------------------------------------------------------------
10    FORMAT('ERROR in routine SPLINT: Bad XA input.')

      RETURN
      END

*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
      END

*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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

      DOUBLE PRECISION eps
      PARAMETER (eps=1.D-10)

      DOUBLE PRECISION pi,radfac
      PARAMETER (pi=3.141592653589793D0,radfac=pi/180.D0)

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

*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      
      SUBROUTINE rdfous_planet(foufile,nfou,
     .                         nmat,nmugs,xmu,rfou)

*----------------------------------------------------------------------------
*     Open and read the Fourier coefficients file:
*----------------------------------------------------------------------------
      IMPLICIT NONE

      DOUBLE PRECISION eps
      PARAMETER (eps=1.D-10)

      DOUBLE PRECISION pi,radfac
      PARAMETER (pi=3.141592653589793D0,radfac=pi/180.D0)

      INTEGER iunf 
      PARAMETER (iunf=23)

      INTEGER m,i1,i2,i,j,ki,ibase,nfou,nmat,nmugs,ifou

      DOUBLE PRECISION xmu(nmugs),
     .                 rfou(nmat*nmugs,nmugs,0:nfou)

      CHARACTER ch*1,foufile*100      
      
Cf2py intent(in) foufile, trans, nfou, nmat, nmugs
Cf2py intent(out) xmu, rfou
Cf2py depend(nmat) rfou
Cf2py depend(nfou) rfou
Cf2py depend(nmugs) rfou, xmu
*----------------------------------------------------------------------------
*     Initialize the aray:
*----------------------------------------------------------------------------
      DO i=1,nmat*nmugs
         DO j=1,nmugs
            DO m=0,nfou
               rfou(i,j,m)= 0.D0
            ENDDO
         ENDDO
      ENDDO

*----------------------------------------------------------------------------
*     Open the Fourier coefficients file:
*----------------------------------------------------------------------------
      OPEN(unit=iunf,file=foufile,status='old',err=999)

*----------------------------------------------------------------------------
*     Read the header (the dap.in file contents):
*----------------------------------------------------------------------------
2     READ(iunf,'(a1)',err=997,end=998) ch
      IF (ch.EQ.'#') GOTO 2
      BACKSPACE(iunf)

*----------------------------------------------------------------------------
*     Read the accuracy, the matrix size, the number of Gauss points and the 
*     Gaussian integration points:
*----------------------------------------------------------------------------
      READ(iunf,*) nmat
      READ(iunf,*) nmugs
      DO i=1,nmugs
         READ(iunf,*) xmu(i)
      ENDDO

      ifou=0
20    DO i=1,nmugs
         ibase= (i-1)*nmat
         DO j=1,nmugs
            READ(iunf,*,END=21) m,i1,i2, 
     .                         (rfou(ibase+ki,j,ifou),ki=1,nmat)

         ENDDO
      ENDDO
      ifou= ifou+1
      GOTO 20
21    CLOSE(iunf)

*----------------------------------------------------------------------------
      GOTO 1000

*----------------------------------------------------------------------------
999   WRITE(*,*) 'Error opening Fourier file!'
      STOP
998   WRITE(*,*) 'Error: unexpected end of the Fourier file!'
      STOP
997   WRITE(*,*) 'Error reading the Fourier file!'
      STOP

*----------------------------------------------------------------------------
1000  RETURN
      END

*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      
      SUBROUTINE rdfous_ring(foufile,trans,nfou,
     .                       nmat,nmugs,xmu,rfou)

*----------------------------------------------------------------------------
*     Open and read the Fourier coefficients file:
*----------------------------------------------------------------------------
      IMPLICIT NONE

      DOUBLE PRECISION eps
      PARAMETER (eps=1.D-10)

      DOUBLE PRECISION pi,radfac
      PARAMETER (pi=3.141592653589793D0,radfac=pi/180.D0)

      INTEGER iunf 
      PARAMETER (iunf=23)

      INTEGER m,i1,i2,i,j,ki,ibase,nfou,nmat,nmugs,ifou

      DOUBLE PRECISION xmu(nmugs),
     .                 rfou(nmat*nmugs,nmugs,0:nfou),
     .                 dummy(nmat)

      CHARACTER ch*1,foufile*100
      
      LOGICAL trans
      
Cf2py intent(in) foufile, trans, nfou, nmat, nmugs
Cf2py intent(out) xmu, rfou
Cf2py depend(nmat) rfou
Cf2py depend(nfou) rfou
Cf2py depend(nmugs) rfou, xmu
*----------------------------------------------------------------------------
*     Initialize the aray:
*----------------------------------------------------------------------------
      DO i=1,nmat*nmugs
         DO j=1,nmugs
            DO m=0,nfou
               rfou(i,j,m)= 0.D0
            ENDDO
         ENDDO
      ENDDO
      
      DO i=1,nmat
         dummy(i) = 0.D0
      ENDDO

*----------------------------------------------------------------------------
*     Open the Fourier coefficients file:
*----------------------------------------------------------------------------
      OPEN(unit=iunf,file=foufile,status='old',err=999)

*----------------------------------------------------------------------------
*     Read the header (the dap.in file contents):
*----------------------------------------------------------------------------
2     READ(iunf,'(a1)',err=997,end=998) ch
      IF (ch.EQ.'#') GOTO 2
      BACKSPACE(iunf)

*----------------------------------------------------------------------------
*     Read the accuracy, the matrix size, the number of Gauss points and the 
*     Gaussian integration points:
*----------------------------------------------------------------------------
      READ(iunf,*) nmat
      READ(iunf,*) nmugs
      DO i=1,nmugs
         READ(iunf,*) xmu(i)
      ENDDO

      ifou=0
20    DO i=1,nmugs
         ibase= (i-1)*nmat
         DO j=1,nmugs
            IF ( trans ) THEN
                READ(iunf,*,END=21) m,i1,i2, (dummy(ki),ki=1,nmat), 
     .                          (rfou(ibase+ki,j,ifou),ki=1,nmat)
            ELSE
                READ(iunf,*,END=21) m,i1,i2, 
     .                              (rfou(ibase+ki,j,ifou),ki=1,nmat), 
     .                              (dummy(ki),ki=1,nmat)
            ENDIF
         ENDDO
      ENDDO
      ifou= ifou+1
      GOTO 20
21    CLOSE(iunf)

*----------------------------------------------------------------------------
      GOTO 1000

*----------------------------------------------------------------------------
999   WRITE(*,*) 'Error opening Fourier file!'
      STOP
998   WRITE(*,*) 'Error: unexpected end of the Fourier file!'
      STOP
997   WRITE(*,*) 'Error reading the Fourier file!'
      STOP

*----------------------------------------------------------------------------
1000  RETURN
      END
