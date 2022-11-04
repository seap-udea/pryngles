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
