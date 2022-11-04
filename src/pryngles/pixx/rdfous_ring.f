      SUBROUTINE rdfous_ring(foufile,trans,nfou,
     .                       nmat,nmugs,xmu,rfou)

*----------------------------------------------------------------------------
*     Open and read the Fourier coefficients file:
*----------------------------------------------------------------------------
      IMPLICIT NONE

      INCLUDE 'max_incl'

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
