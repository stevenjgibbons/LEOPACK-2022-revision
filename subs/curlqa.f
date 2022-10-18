C*********************************************************************
C subroutine CURL of Qst Array ***************************************
C            ----    -   -     ***************************************
C Steve Gibbons  4. 8.99                                             C
C____________________________________________________________________C
C                                                                    C
C Given a full QST decomposition ( at all radial nodes! ), CURLQA    C
C will put the curl of that vector into a second such array.         C
C                                                                    C
C____________________________________________________________________C
C                                                                    C
C Input variables :-                                                 C
C ===============                                                    C
C  Integer                                                           C
C  -------                                                           C
C                                                                    C
C     LH        : Maximum spherical harmonic degree.                 C
C     NR        : Number of radial grid nodes.                       C
C                                                                    C
C  Double Precision                                                  C
C  ----------------                                                  C
C                                                                    C
C     QSTA1     : QST decomp. Dimension ( LH*(LH+2), 3, NR )         C
C     QSTA2     : QST decomp. Dimension ( LH*(LH+2), 3, NR )         C
C                                                                    C
C  The curl of the vector stored in QSTA1 is returned in QSTA2.      C
C                                                                    C
C     RI        : Radius of inner boundary.                          C
C     RO        : Radius of outer boundary.                          C
C                                                                    C
C  Character                                                         C
C  ---------                                                         C
C     ORD       : *(2). Determines the order of accuracy. Options :  C
C             SS - Strictly second order                             C
C             SF - Strictly fourth order                             C
C             O5 - Optimum accuracy for bandwidth 5; this gives      C
C                  Fourth order accuracy for 1st and 2nd derivatives C
C                  and second order accuracy for 3rd and 4th der.s   C
C             O7 - Optimum accuracy for bandwidth 7; this gives      C
C                  Sixth order accuracy for 1st and 2nd derivatives  C
C                  and fourth order accuracy for 3rd and 4th der.s   C
C                                                                    C
C____________________________________________________________________C
C
C*********************************************************************
      SUBROUTINE CURLQA ( LH, NR, QSTA1, QSTA2, RI, RO, ORD )
      IMPLICIT NONE
C____________________________________________________________________C
C Variable declarations - Parameters ................................C
      INTEGER LH, NR
      DOUBLE PRECISION QSTA1( LH*(LH+2), 3, NR ),
     1                 QSTA2( LH*(LH+2), 3, NR ), RI, RO
      CHARACTER *(2) ORD
C____________________________________________________________________C
C Variable declarations - Working variables .........................C
      INTEGER IH, NH, IRAD, L, M, ICS, IQST
      DOUBLE PRECISION ZERO, LOW, H, TOTAL, SQRLL1, RAD,
     1                 D0F, D1F, D2F, D3F, D4F
      PARAMETER ( ZERO = 0.0d0, LOW = 1.0d-7 )
C____________________________________________________________________C
C START OF PROGRAM **************************************************C
C____________________________________________________________________C
C
      NH = LH*(LH + 2)
C
C Clear the array QSTA2
C
      DO IRAD = 1, NR
        DO IQST = 1, 3
          DO IH = 1, NH
            QSTA2( IH, IQST, IRAD ) = ZERO
          ENDDO
        ENDDO
      ENDDO
C
      H = ( RO - RI )/DBLE( NR - 1 )
C
C Let's loop around the harmonics
C
      DO IH = 1, NH
        CALL LMFIND ( IH, L, M, ICS )
        SQRLL1 = DBLE( L + L*L )
        SQRLL1 = DSQRT( SQRLL1 )
        DO IQST = 1, 3
          TOTAL = ZERO
          DO IRAD = 1, NR
            TOTAL = TOTAL + DABS( QSTA1( IH, IQST, IRAD ) )
          ENDDO
c
c    exit this harmonic if it is simply zero
c
          IF ( ABS(TOTAL).LT.LOW ) GOTO 500
c
          DO IRAD = 1, NR
            RAD = RI + H*DBLE( IRAD - 1 )
C ----------                     dealing with scaloidal
            IF ( IQST.EQ.1 ) THEN
              QSTA2( IH, 3, IRAD ) =
     1          (-1.0d0)*SQRLL1*QSTA1( IH, 1, IRAD )/RAD
            ENDIF
C ----------                     dealing with spheroidal
            IF ( IQST.EQ.2 ) THEN
              CALL QSTADF ( LH, IH, NR, IRAD, IQST, QSTA1, ORD,
     1                    D0F, D1F, D2F, D3F, D4F, RI, RO )
              QSTA2( IH, 3, IRAD ) = QSTA2( IH, 3, IRAD ) +
     1              D1F + D0F/RAD
            ENDIF
C ----------                     dealing with toroidal
            IF ( IQST.EQ.3 ) THEN
              CALL QSTADF ( LH, IH, NR, IRAD, IQST, QSTA1, ORD,
     1                    D0F, D1F, D2F, D3F, D4F, RI, RO )
              QSTA2( IH, 1, IRAD ) = 
     1           (-1.0d0)*SQRLL1*D0F/RAD
              QSTA2( IH, 2, IRAD ) = 
     1           (-1.0d0)*( D1F + D0F/RAD )
            ENDIF
C ---------- 
          ENDDO
c
 500      CONTINUE
        ENDDO
      ENDDO
C
      RETURN
      END
C*********************************************************************

