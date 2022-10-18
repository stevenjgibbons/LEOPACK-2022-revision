C*********************************************************************
C subroutine Dynamo Benchmark Initial State Vector Form **************
C            -      -         -       -     -      -    **************
C Steve Gibbons Fri Feb  4 08:36:23 GMT 2000                         C
C____________________________________________________________________C
C                                                                    C
C Fills in vector with the initial state for dynamo benchmark        C
C as suggested by U. Christensen. Values are given by the function   C
C DBPICS.                                                            C
C                                                                    C
C____________________________________________________________________C
C                                                                    C
C Input variables :-                                                 C
C ===============                                                    C
C  Integer                                                           C
C  -------                                                           C
C                                                                    C
C     INARR     : Indexes the array VEC. Dim (3). See INDFUN.        C
C                                                                    C
C     MHT       : Dim (*). MHT(ih) = 1 --> poloidal velocity.        C
C                          MHT(ih) = 2 --> toroidal velocity.        C
C                          MHT(ih) = 3 --> temperature.              C
C                          MHT(ih) = 4 --> poloidal magnetic field   C
C                          MHT(ih) = 5 --> toroidal magnetic field   C
C                                                                    C
C     MHL       : Dim (*). MHL(ih) = l, Spherical harmonic degree.   C
C                                                                    C
C     MHM       : Dim (*). MHM(ih) = m,  for cos m phi dependence.   C
C                          MHM(ih) = -m, for sin m phi dependence.   C
C                                                                    C
C     IFLAG     : 1 --> full system including magnetic field.        C
C                 2 --> velocity and temperature only.               C
C                                                                    C
C  Double Precision                                                  C
C  ----------------                                                  C
C                                                                    C
C     VEC       : Dim ( * ). Vector containing initial solution.     C
C                                                                    C
C     XARR      : Dim ( * ). Vector containing grid node radii.      C
C                                                                    C
C____________________________________________________________________C
C
C*********************************************************************
      SUBROUTINE DBISVF( INARR, MHT, MHL, MHM, IFLAG, VEC, XARR )
      IMPLICIT NONE
C____________________________________________________________________C
C Variable declarations - Parameters ................................C
      INTEGER INARR( 3 ), MHT( * ), MHL( * ), MHM( * ), IFLAG
      DOUBLE PRECISION VEC( * ), XARR( * )
C____________________________________________________________________C
C Variable declarations - Working variables .........................C
      INTEGER LH, NTH, NPH, MMAX, NPHPTS
      PARAMETER ( LH = 4, NTH = 6, NPH = 16 )
      DOUBLE PRECISION SHC( LH*( LH + 2) ), ZCOEF, X1, X2,
     1                 SF( NPH, NTH ), GAUW( NTH ), GAUX( NTH )
      DOUBLE PRECISION PA ( ( LH + 1 )*( LH + 2 )/2 , NTH ),
     1                 DPA ( ( LH + 1 )*( LH + 2 )/2 , NTH )
      DOUBLE PRECISION FTF1( 2*NPH ), FTF2( 2*NPH ), FTF3( 2*NPH ),
     1                 VF( 2, NTH, 3)
      DOUBLE PRECISION QST( LH*( LH + 2), 3 ), RI, RO, RAD, LOW,
     1                 THE, PHI, COSTH, DPHI, PI, ZERO, DBPISC, COEF,
     2                 SQRLL1
      PARAMETER ( PI=3.14159265358979312D0, ZERO = 0.0d0,
     1            X1 = -1.0d0, X2 =  1.0d0, LOW = 1.0d-7 )
      INTEGER IR, NR, NH, IOP, ILEN, ITH, IPH, ICOMP, IHARM, NHARM,
     1        L, M, ICS, IH, MM, IND, INDFUN, ILN
C____________________________________________________________________C
C START OF PROGRAM **************************************************C
C____________________________________________________________________C
C
      NR     = INARR( 2 )
      NH     = INARR( 3 )
C
      RI     = XARR(  1 )
      RO     = XARR( NR )
      NHARM  = LH*(LH + 2)
C
      IF ( DABS( RI ).LT.LOW ) THEN
        ILN = 2
      ELSE
        ILN = 1
      ENDIF
C
C Zero vector
C
      IOP  = 0
      ILEN = NR*NH
      CALL VECOP( VEC, ZERO, ILEN, IOP )
C
C Calculate Gauss points, Legendre functions etc.
C
      CALL GAUWTS ( X1, X2, GAUX, GAUW, NTH )
      CALL SCHNLA ( PA, DPA, GAUX, LH, NTH )
C
C Do initial temperature field
C
      ICOMP  = 4
      MMAX   = 4
      NPHPTS = NPH
      DPHI   = 2.0d0*PI/DBLE( NPHPTS )
      DO  IR = ILN, NR
        RAD  = XARR( IR )
C       .
C       . Fill scalar function array.
C       .
        DO ITH = 1, NTH
          COSTH = GAUX( ITH )
          THE   = ACOS( COSTH )
          DO IPH = 1, NPHPTS
            PHI  = DBLE( IPH - 1 )*DPHI
            SF( IPH, ITH ) = DBPISC( ICOMP, RAD, THE, PHI, RI, RO )
          ENDDO
        ENDDO
C       .
C       . SF is filled. Convert to spectral coefficients
C       .
        CALL FORSST ( SHC, SF, GAUW, PA, FTF1, LH, MMAX,
     1                NTH, NPHPTS, ZCOEF )
C       .
        DO IHARM = 1, NHARM
          COEF   = SHC( IHARM )
          IF ( DABS( COEF ).GT.LOW ) THEN
            CALL LMFIND( IHARM, L, M, ICS )
            IF ( ICS.EQ.1 ) MM = M
            IF ( ICS.EQ.2 ) MM = -M
            DO IH = 1, NH
              IF ( MHT( IH ).EQ.3 .AND. MHL( IH ).EQ.L
     1             .AND. MHM( IH ).EQ.MM ) THEN
                IND = INDFUN( IR, IH, INARR )
                VEC( IND ) = COEF
              ENDIF
            ENDDO
          ENDIF
        ENDDO
C       .
      ENDDO
C
      IF ( IFLAG.EQ.2 ) RETURN
C
C Do initial magnetic field.
C
      MMAX   = 0
      NPHPTS = 2
      DO  IR = ILN, NR
        RAD  = XARR( IR )
C       .
C       . Fill vector function array.
C       .
        DO ITH = 1, NTH
          COSTH = GAUX( ITH )
          THE   = ACOS( COSTH )
          DO ICOMP = 1, 3
            VF( 1, ITH, ICOMP ) =
     1                DBPISC( ICOMP, RAD, THE, PHI, RI, RO )
            VF( 2, ITH, ICOMP ) = VF( 1, ITH, ICOMP )
          ENDDO
        ENDDO
C       .
C       . VF is filled. Convert to spectral coefficients
C       .
        CALL VF2QST ( QST, VF, GAUX, GAUW, PA, DPA, FTF1, FTF2,
     1                FTF3, ZCOEF, LH, NTH, NPHPTS, MMAX )
C       .
        DO IHARM = 1, NHARM
C         .
C         . First do poloidal parts
C         .
          COEF   = QST( IHARM, 1 )
          IF ( DABS( COEF ).GT.LOW ) THEN
            CALL LMFIND( IHARM, L, M, ICS )
            IF ( ICS.EQ.1 ) MM = M
            IF ( ICS.EQ.2 ) MM = -M
            DO IH = 1, NH
              IF ( MHT( IH ).EQ.4 .AND. MHL( IH ).EQ.L
     1             .AND. MHM( IH ).EQ.MM ) THEN
                IND = INDFUN( IR, IH, INARR )
                VEC( IND ) = COEF*RAD/DBLE( L*L + L )
              ENDIF
            ENDDO
          ENDIF
C         .
C         . Now do toroidal parts
C         .
          COEF   = QST( IHARM, 3 )
          IF ( DABS( COEF ).GT.LOW ) THEN
            CALL LMFIND( IHARM, L, M, ICS )
            IF ( ICS.EQ.1 ) MM = M
            IF ( ICS.EQ.2 ) MM = -M
            DO IH = 1, NH
              IF ( MHT( IH ).EQ.5 .AND. MHL( IH ).EQ.L
     1             .AND. MHM( IH ).EQ.MM ) THEN
                IND = INDFUN( IR, IH, INARR )
                VEC( IND ) = (-1.0d0)*COEF/SQRLL1( L )
              ENDIF
            ENDDO
          ENDIF
C         .
        ENDDO
C       .
      ENDDO
C
      IF ( IFLAG.EQ.1 ) RETURN
C
      PRINT *,' Subroutine DBPSVF.'
      PRINT *,' IFLAG is neither 1 nor 2.'
      PRINT *,' Program aborted.'
      STOP
      END
C*********************************************************************
