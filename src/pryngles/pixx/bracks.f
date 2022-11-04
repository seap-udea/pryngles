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

      INCLUDE 'max_incl'

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
