C*********************************************************************
C subroutine EigenVECtor EXtract *************************************
C            -    ---    --      *************************************
C Steve Gibbons Sat Nov 20 15:54:46 GMT 1999                         C
C____________________________________________________________________C
C                                                                    C
C The routine VMEPS (and all others which call the Arnoldi           C
C eigenvalue routines) return the eigenvalues in an array V of       C
C dimensions ( N2, NSD ) with the (up to) NSD eigenvalues being      C
C stored in the columns of V.                                        C
C                                                                    C
C EVECEX simply takes a number I between 1 and NSD and returns       C
C the eigenvector I in the single dimension array VEC.               C
C                                                                    C
C If N2 and NSD are the same dimensions which VMEPS receives for V   C
C then EVECEX will be helpful in avoiding dimensioning problems when C
C retrieving eigenvectors.                                           C
C                                                                    C
C____________________________________________________________________C
C                                                                    C
C  Integer                                                           C
C  -------                                                           C
C                                                                    C
C     N2        : Leading dimension of V in call to VMEPS            C
C     NSD       : Second dimension of V in call to VMEPS             C
C     I         : Number of vector to be retrieved.                  C
C                                                                    C
C  Double Precision                                                  C
C  ----------------                                                  C
C                                                                    C
C     V         : Dimension ( N2, NSD ). Array containing vectors.   C
C     VEC       : Dimension ( N2 ). Output containing vec. number I. C
C                                                                    C
C                                                                    C
C____________________________________________________________________C
C
C*********************************************************************
      SUBROUTINE EVECEX( N2, NSD, I, V, VEC )
      IMPLICIT NONE
C____________________________________________________________________C
C Variable declarations - Parameters ................................C
      INTEGER N2, NSD, I
      DOUBLE PRECISION V( N2, NSD ), VEC( N2 )
C____________________________________________________________________C
C Variable declarations - Working variables .........................C
C
      INTEGER J
C____________________________________________________________________C
C START OF PROGRAM **************************************************C
C____________________________________________________________________C
C
      IF ( I.LT.1 .OR. I.GT.NSD ) THEN
        PRINT *,' Subroutine EVECEX.'
        PRINT *,' I   = ', I
        PRINT *,' NSD = ', NSD
        PRINT *,' Program aborted.'
        STOP
      ENDIF
C     .
      DO J = 1, N2
        VEC( J ) = V( J, I )
      ENDDO
C     .
      RETURN
      END
C*********************************************************************
