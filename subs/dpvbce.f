C*********************************************************************
C subroutine Double Precision Vector Boundary Condition Enforce ******
C            -      -         -      -        -         -       ******
C Steve Gibbons 22.3.99                                              C
C____________________________________________________________________C
C Treats a vector prior to algebraic solution in order to enforce    C
C the boundary conditions. (The matrix must have been treated with   C
C DPMBCE prior to LU decomposition.)                                 C
C____________________________________________________________________C
C Input variables :-                                                 C
C ===============                                                    C
C  Integer                                                           C
C  -------                                                           C
C     N2        : Length of the vector.                              C
C     NH        : Number of spherical harmonics.                     C
C     NR        : Number of radial grid nodes.                       C
C     MHT       : Dimension ( NH )                                   C
C                                                                    C
C   MHT( I ) = 1 for a poloidal velocity vector                      C
C   MHT( I ) = 2 for a toroidal velocity vector                      C
C   MHT( I ) = 3 for a temperature / codensity term                  C
C                                                                    C
C                                                                    C
C  Double Precision                                                  C
C  ----------------                                                  C
C     RVEC      : Vector dimension ( N2 )                            C
C                                                                    C
C     CIB       : Coefficients for inner boundary - dimension ( * )  C
C     COB       : Coefficients for outer boundary - dimension ( * )  C
C   ( cib and cob are only referred to if OIB is true )              C
C                                                                    C
C  Logical                                                           C
C  -------                                                           C
C     OIB       : Set to .TRUE. for inhomogeneous boundaries         C
C                 Set to .FALSE. for homogeneous boundaries          C
C                                                                    C
C____________________________________________________________________C
C
C*********************************************************************
      SUBROUTINE DPVBCE ( N2, NH, NR, MHT, RVEC, CIB, COB, OIB )
      IMPLICIT NONE
C____________________________________________________________________C
C Variable declarations - Parameters ................................C
      INTEGER N2, NH, NR, MHT( NH )
      DOUBLE PRECISION RVEC( N2 ), CIB( * ), COB( * )
      LOGICAL OIB
C____________________________________________________________________C
C Variable declarations - Working variables .........................C
      INTEGER IRN, IH, IPTT, IND
      DOUBLE PRECISION ZERO, FAC
      PARAMETER ( ZERO = 0.0d0 )
C____________________________________________________________________C
C START OF PROGRAM **************************************************C
C____________________________________________________________________C
C
C Check the length of the vector
C
      IF ( N2.NE.NR*NH ) THEN
         PRINT *,' Subroutine DPVBCE. Bad dimensions '
         PRINT *,' Array length = ',N2
         PRINT *,' NR = ', NR,' NH = ',NH
         STOP
      ENDIF
C
      DO IH = 1, NH
        IPTT = MHT( IH )
        FAC = ZERO
        IF ( OIB ) FAC = CIB( IH )
        IRN = 1
        IND = ( IRN - 1 )*NH + IH
        RVEC( IND ) = FAC
        FAC = ZERO
        IF ( OIB ) FAC = COB( IH )
        IRN = NR
        IND = ( IRN - 1 )*NH + IH
        RVEC( IND ) = FAC
        IF ( IPTT.EQ.1 ) THEN
          FAC = ZERO
          IRN = 2
          IND = ( IRN - 1 )*NH + IH
          RVEC( IND ) = FAC
          IRN = NR - 1
          IND = ( IRN - 1 )*NH + IH
          RVEC( IND ) = FAC
        ENDIF
      ENDDO
C
      RETURN
      END
C*********************************************************************
