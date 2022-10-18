C*********************************************************************
C subroutine Arbitrarily Spaced Solution Vector Addition Routine *****
C            -           -      -        -      -        -       *****
C Steve Gibbons Tue May 30 11:39:09 BST 2000                         C
C____________________________________________________________________C
C                                                                    C
C Let VEC1 be a solution vector indexed by INARR1, MHT1, MHL1, MHM1  C
C and with radial spacings given by XARR1.                           C
C                                                                    C
C Let VEC2 be a solution vector indexed by INARR2, MHT2, MHL2, MHM2  C
C and with radial spacings given by XARR2.                           C
C                                                                    C
C ASSVAR will loop around all the harmonics that constitute VEC2,    C
C and if VEC1 and VEC2 have a harmonic in common then we loop around C
C the grid nodes in XARR1 and interpolate (using SVRINT) VEC2 and    C
C perform                                                            C
C                                                                    C
C  VEC1( IR, IH ) :=  FAC1*VEC1( IR, IH ) + FAC2*VEC2( IR, IH )      C
C                                                                    C
C____________________________________________________________________C
C                                                                    C
C Input variables :-                                                 C
C ===============                                                    C
C  Integer                                                           C
C  -------                                                           C
C                                                                    C
C     INARR1    : Integer array dimension ( * ).                     C
C                  Elements may be arbitrary except for              C
C                  INARR1( 1 ) = IFORMF - flag for vector format.    C
C                  INARR1( 2 ) = NR. Number of radial grid nodes.    C
C                  INARR1( 3 ) = NH = total number of radial func.s  C
C                                                                    C
C     MHT1      : Harmonic types for VEC1. Dim (*).                  C
C     MHL1      : Harmonic degree for VEC1. Dim (*).                 C
C     MHM1      : Harmonic order for VEC1. Dim (*).                  C
C                                                                    C
C     INARR2    : Integer array dimension ( * ).                     C
C                  Elements may be arbitrary except for              C
C                  INARR2( 1 ) = IFORMF - flag for vector format.    C
C                  INARR2( 2 ) = NR. Number of radial grid nodes.    C
C                  INARR2( 3 ) = NH = total number of radial func.s  C
C                                                                    C
C     MHT2      : Harmonic types for VEC2. Dim (*).                  C
C     MHL2      : Harmonic degree for VEC2. Dim (*).                 C
C     MHM2      : Harmonic order for VEC2. Dim (*).                  C
C                                                                    C
C     NNDS      : Number of nodes to be used in interpolation.       C
C                 Must be atleast 2.                                 C
C     IWORK     : Dimension ( NNDS ). Work array.                    C
C                                                                    C
C  Double Precision                                                  C
C  ----------------                                                  C
C                                                                    C
C     VEC1      : Solution vector 1. Dim ( * )                       C
C     XARR1     : Radial spacings 1. Dim ( * )                       C
C     VEC2      : Solution vector 2. Dim ( * )                       C
C     XARR2     : Radial spacings 2. Dim ( * )                       C
C                                                                    C
C     WORKA     : Work array. Dim (NNDS)                             C
C     WORKB     : Work array. Dim (NNDS)                             C
C     WORKC     : Work array. Dim (NNDS,NNDS)                        C
C                                                                    C
C     FAC1      : Scaling factor 1.                                  C
C     FAC2      : Scaling factor 2.                                  C
C                                                                    C
C____________________________________________________________________C
C
C*********************************************************************
      SUBROUTINE ASSVAR( INARR1, MHT1, MHL1, MHM1, INARR2, MHT2, MHL2,
     1                   MHM2, NNDS, IWORK, VEC1, XARR1, VEC2, XARR2,
     2                   WORKA, WORKB, WORKC, FAC1, FAC2 )
      IMPLICIT NONE
C____________________________________________________________________C
C Variable declarations - Parameters ................................C
      INTEGER INARR1( * ), MHT1( * ), MHL1( * ), MHM1( * ),
     1        INARR2( * ), MHT2( * ), MHL2( * ), MHM2( * ),
     2        NNDS, IWORK( NNDS )
      DOUBLE PRECISION VEC1( * ), XARR1( * ), VEC2( * ), XARR2( * ),
     1        WORKA( NNDS ), WORKB( NNDS ), WORKC( NNDS, NNDS ),
     2        FAC1, FAC2
C____________________________________________________________________C
C Variable declarations - Working variables .........................C
      INTEGER IR, IH1, IH2, NR, IND, INDFUN, NH1, NH2
      DOUBLE PRECISION RAD
C____________________________________________________________________C
C START OF PROGRAM **************************************************C
C____________________________________________________________________C
C
      NR     = INARR1( 2 )
      NH1    = INARR1( 3 )
      NH2    = INARR2( 3 )
C
      DO IH1 = 1, NH1
        DO IH2 = 1, NH2
          IF (    MHT1( IH1 ).EQ.MHT2( IH2 )   .AND.
     1            MHL1( IH1 ).EQ.MHL2( IH2 )   .AND.
     2            MHM1( IH1 ).EQ.MHM2( IH2 )          ) THEN
C           .
C           . Loop around grid nodes in VEC1
C           .
            DO IR = 1, NR
              RAD = XARR1( IR )
              IND = INDFUN( IR, IH1, INARR1 )
              CALL SVRINT( RAD, VEC2, XARR2, INARR2, IH2, NNDS,
     1                     WORKA, IWORK, WORKB, WORKC )
              VEC1( IND ) =  FAC1*VEC1( IND ) + FAC2*WORKA( 1 )
            ENDDO
C           .
            GOTO 50
          ENDIF
        ENDDO
 50     CONTINUE
      ENDDO
C
      RETURN
      END
C*********************************************************************
