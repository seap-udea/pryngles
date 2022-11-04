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
