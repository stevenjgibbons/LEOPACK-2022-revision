C*********************************************************************
C subroutine Poloidal Velocity Completeness Coefficients Find ********
C            -        -        -            -            -    ********
C Steve Gibbons Thu Nov  9 10:10:59 WET 2000                         C
C____________________________________________________________________C
C                                                                    C
C Because of the additional boundary conditions on the poloidal      C
C velocity radial function, the values at grid nodes 2 and NR-1 are  C
C not solved for in the solution vector.                             C
C                                                                    C
C This routine extracts from the array SVFDC, two vectors of length  C
C NBN, called PVLC (poloidal velocity left coefficients) and PVRC    C
C (poloidal velocity right coefficients) such that for poloidal      C
C harmonic IH, the vectors can be completed with the following       C
C segments of code:-                                                 C
C                                                                    C
C     IRND = 2                                                       C
C     IND2 = INDFUN( IRND, IH, INARR )                               C
C     TEMP = 0.0d0                                                   C
C     DO I = 1, NBN                                                  C
C       IR   = IRND + I                                              C
C       IND  = INDFUN( IR, IH, INARR )                               C
C       TEMP = TEMP + PVLC( I )*SV( IND )                            C
C     ENDDO                                                          C
C     SV( IND2 ) = TEMP                                              C
C                                                                    C
C     IRND = NR - 1                                                  C
C     IND2 = INDFUN( IRND, IH, INARR )                               C
C     TEMP = 0.0d0                                                   C
C     DO I = 1, NBN                                                  C
C       IR   = IRND - I                                              C
C       IND  = INDFUN( IR, IH, INARR )                               C
C       TEMP = TEMP + PVRC( I )*SV( IND )                            C
C     ENDDO                                                          C
C     SV( IND2 ) = TEMP                                              C
C                                                                    C
C____________________________________________________________________C
C                                                                    C
C Input variables :-                                                 C
C ===============                                                    C
C                                                                    C
C  Integer                                                           C
C  -------                                                           C
C                                                                    C
C     NR        : Number of grid nodes.                              C
C                                                                    C
C     NDCS      : Number of distinct differencing coeff.s            C
C                  represented in SVFDC.                             C
C                                                                    C
C     NBN       : Maximum number of nodes on either side for         C
C                  central differencing.                             C
C                                                                    C
C     NFDCM     : Leading order of the array SVFDC.                  C
C                 This must be atleast (2*NBN + 1)                   C
C                                                                    C
C     NDRVM     : Maximum number of derivatives allowed.             C
C                                                                    C
C     IS        : Number of difference scheme for poloidal velocity. c
C                                                                    C
C  Double Precision                                                  C
C  ----------------                                                  C
C                                                                    C
C     SVFDC     : Finite difference coefficient matrix.              C
C                  Dimension ( NFDCM, NR, NDRVM+1, NDCS ).           C
C                                                                    C
C     PVLC      : Dim (NBN). See above.                              C
C     PVRC      : Dim (NBN). See above.                              C
C____________________________________________________________________C
C
C*********************************************************************
      SUBROUTINE PVCCF( NR, NDCS, NBN, NFDCM, NDRVM, IS, SVFDC,
     1                  PVLC, PVRC )
      IMPLICIT NONE
C____________________________________________________________________C
C Variable declarations - Parameters ................................C
      INTEGER          NR, NDCS, NBN, NFDCM, NDRVM, IS
      DOUBLE PRECISION SVFDC( NFDCM, NR, NDRVM+1, NDCS ),
     1                 PVLC( NBN ), PVRC( NBN )
C____________________________________________________________________C
C Variable declarations - Working variables .........................C
      INTEGER I, IR, IND
C____________________________________________________________________C
C START OF PROGRAM **************************************************C
C____________________________________________________________________C
C
      IR  = 2
      DO I = 1, NBN
        IND       = NBN + 1 + I
        PVLC( I ) = SVFDC( IND, IR, 1, IS )
      ENDDO
C     .
      IR  = NR - 1
      DO I = 1, NBN
        IND       = NBN + 1 - I
        PVRC( I ) = SVFDC( IND, IR, 1, IS )
      ENDDO
C     .
      RETURN
      END
C*********************************************************************
