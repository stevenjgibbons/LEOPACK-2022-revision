C*********************************************************************
C Subroutine Vector Function Dot Product *****************************
C            -      -        -   -       *****************************
C Steve Gibbons Fri Sep 24 16:35:41 BST 1999                         C
C____________________________________________________________________C
C To have two input vector function arrays of dimension              C
C ( NPHPTS , NTHPTS, 3 ) of which the scalar products are taken to   C
C give a scalar vector SF with dimension ( NPHPTS, NTHPTS ).         C
C____________________________________________________________________C
C Input variables :-                                                 C
C ===============                                                    C
C  Integer                                                           C
C  -------                                                           C
C     NPHPTS	: Number of phi points.                              C
C     NTHPTS	: Number of theta points.                            C
C                                                                    C
C  Double Precision                                                  C
C  ----------------                                                  C
C     VF1	: First vector function. Dim ( NPHPTS , NTHPTS, 3 )  C
C     VF2	: Second vector function. Dim ( NPHPTS , NTHPTS, 3 ) C
C____________________________________________________________________C
C Output variables :-                                                C
C ================                                                   C
C                                                                    C
C  Double Precision                                                  C
C  ----------------                                                  C
C     SF	: Scalar function. Dim ( NPHPTS, NTHPTS )            C
C____________________________________________________________________C
C
C*********************************************************************
      SUBROUTINE VFDP ( VF1, VF2, SF, NPHPTS, NTHPTS )
      IMPLICIT NONE
C____________________________________________________________________C
C Variable declarations - Parameters ................................C
      INTEGER NPHPTS, NTHPTS
      DOUBLE PRECISION VF1 ( NPHPTS , NTHPTS, 3 ),
     1                 VF2 ( NPHPTS , NTHPTS, 3 ),
     2                 SF ( NPHPTS, NTHPTS )
C____________________________________________________________________C
C Variable declarations - Working variables .........................C
      INTEGER ITHETA, IPHI, ICOMP, IOP
      DOUBLE PRECISION ZERO
      PARAMETER ( ZERO = 0.0d0 )
C____________________________________________________________________C
C START OF PROGRAM **************************************************C
C____________________________________________________________________C
C
C .... let's set SF to zero.
c
      IOP = 0
      CALL MATOP ( SF, ZERO, NPHPTS, NTHPTS, IOP )
c
C .... let's loop around icomp, then nthpts, then nphpts
      DO ICOMP = 1, 3
         DO ITHETA = 1, NTHPTS
            DO IPHI = 1, NPHPTS
               SF( IPHI, ITHETA ) = SF( IPHI, ITHETA ) +
     1       VF1( IPHI, ITHETA, ICOMP )*VF2( IPHI, ITHETA, ICOMP )
            ENDDO
         ENDDO
      ENDDO
C____________________________________________________________________C
      RETURN
      END
C*********************************************************************

