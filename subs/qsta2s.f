C*********************************************************************
C Subroutine QST Array 2 Solution vector *****************************
C            --- -     - -               *****************************
C Steve Gibbons  4. 8.99                                             C
C____________________________________________________________________C
C                                                                    C
C  Note that this routine adds the result to SV and does not zero    C
C  this array before use.                                            C
C                                                                    C
C____________________________________________________________________C
C                                                                    C
C Input variables :-                                                 C
C ===============                                                    C
C  Integer                                                           C
C  -------                                                           C
C     LH        : Level of harmonics.                                C
C     NDIM      : Length of vector, SV.                              C
C     NR        : Number of radial grid nodes                        C
C     NH        : Number of harmonics (all types)                    C
C     MHT       : Dimension ( NH )                                   C
C     MHL       : Dimension ( NH )                                   C
C     MHM       : Dimension ( NH )                                   C
C     MHC       : Dimension ( NH )                                   C
C                                                                    C
C     IPOLD     : Destination of poloidal harmonics                  C
C     ITORD     : Destination of toroidal harmonics                  C
C     ICURL     : Specifies whether curl is taken before             C
C                  putting into solution vector.                     C
C                  1 = no. 2 = yes                                   C
C                                                                    C
C       The flags ipold and itord determine which components of the  C
C     solution vector the poloidal and toroidal components of the    C
C     qst array should go. The options are                           C
C                                                                    C
C         1 --- Poloidal velocity harmonics                          C
C         2 --- Toroidal velocity harmonics                          C
C         4 --- Poloidal magnetic field harmonics                    C
C         5 --- Toroidal magnetic field harmonics                    C
C                                                                    C
C     For instance, if the QST array contains a decomposition of     C
C ( u \times B ) then set ICURL = 2 and set IPOLD = 4 and ITORD = 5. C
C                                                                    C
C     If the QST array contains a decomposition of curl (u \times B) C
C  then set ICURL = 1 and set IPOLD = 4 and ITORD = 5.               C
C                                                                    C
C     If the QST array contains a decomposition of (u \times curl u) C
C  then set ICURL = 2 and set IPOLD = 2 and ITORD = 1.               C
C  (The poloidal and toroidal parts are reversed because we are      C
C  are considering the vorticity equation.                           C
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
C                                                                    C
C  Double Precision                                                  C
C  ----------------                                                  C
C                                                                    C
C     SV        : Old Solution Vector Dimensions                     C
C                  Dim ( NDIM ) with NDIM = NH * NR                  C
C     QST       : Input array containing scaloidal/spheroidal        C
C                  decomposition of vector c.f. eqn (38).            C
C                  Has dimensions (  LH*(LH+2) , 3, NR ).            C
C     FAC       : Multiplication factor for term.                    C
C     RI        : Radius of inner core.                              C
C     RO        : Radius of outer core.                              C
C____________________________________________________________________C
C
C*********************************************************************
      SUBROUTINE QSTA2S ( LH, NDIM, NR, NH, MHT, MHL, MHM, MHC,
     1                    IPOLD, ITORD, ICURL, SV, QST, FAC, RI, RO,
     2                    ORD )
      IMPLICIT NONE
C____________________________________________________________________C
C Variable declarations - Parameters ................................C
      INTEGER LH, NDIM, NR, NH, MHT( NH ), MHL( NH ),
     1        MHM( NH ), MHC( NH ), IPOLD, ITORD, ICURL
      CHARACTER *(2) ORD
      DOUBLE PRECISION SV( NDIM ), RI, RO, FAC,
     1                 QST( LH*(LH+2), 3, NR )
C____________________________________________________________________C
C Variable declarations - Working variables .........................C
      INTEGER IRN, IH, L, M, ICS, NOHARM, INDSHC, IND
      DOUBLE PRECISION RAD, H, DPLL1, RDPLL1, TOL,
     1                 D0F, D1F, D2F, D3F, D4F
      PARAMETER ( TOL = 1.0d-8 )
C____________________________________________________________________C
C START OF PROGRAM **************************************************C
C____________________________________________________________________C
C
      H = (RO - RI)/DBLE( NR - 1 )
C
C Check value of NDIM against NH and NR etc.
C
      IF ( NDIM.NE.NH*NR ) THEN
        PRINT *,' Subroutine QSTA2S, bad array size'
        PRINT *,' NDIM = ',NDIM
        PRINT *,' NH = ',NH,' NR = ',NR
        STOP
      ENDIF
C
C Early exit for zero factor
C
      IF ( ABS( FAC ).LT.TOL ) RETURN
C
C Check values of IPOLD, ITORD and ICURL
C
      IF ( IPOLD.NE.1 .AND. IPOLD.NE.2 .AND. IPOLD.NE.4 .AND.
     1     IPOLD.NE.5 ) THEN
        PRINT *,' Subroutine QSTA2S, IPOLD = ', IPOLD
        STOP
      ENDIF
C
      IF ( ITORD.NE.1 .AND. ITORD.NE.2 .AND. ITORD.NE.4 .AND.
     1     ITORD.NE.5 ) THEN
        PRINT *,' Subroutine QSTA2S, ITORD = ', ITORD 
        STOP
      ENDIF
C
      IF ( ICURL.NE.1 .AND. ICURL.NE.2 ) THEN
        PRINT *,' Subroutine QSTA2S, ICURL = ', ICURL 
        STOP
      ENDIF
C
      DO IRN = 1, NR
        RAD = RI + H*DBLE( IRN - 1 )
        DO IH = 1, NH
           L   = MHL( IH )
           M   = MHM( IH )
           ICS = MHC( IH )
           NOHARM = INDSHC( L, M, ICS )
           DPLL1 = DBLE( L*L + L )
           RDPLL1 = DSQRT( DPLL1 )
           IND = ( IRN - 1 )*NH + IH
c...........
C                           destination poloidal harmonics
c...........
           IF ( MHT( IH ).EQ.IPOLD ) THEN
C
C First do case where curl does not need to be taken
C
             IF ( ICURL.EQ.1 ) THEN
               SV( IND ) = SV( IND ) +
     1             FAC*RAD*QST( NOHARM, 1, IRN )/DPLL1
             ENDIF
C
C Now do case where curl needs to be taken
C
             IF ( ICURL.EQ.2 ) THEN
               SV( IND ) = SV( IND ) -
     1             FAC*QST( NOHARM, 3, IRN )/RDPLL1
             ENDIF
C
           ENDIF
c...........
C                           destination toroidal harmonics
c...........
           IF ( MHT( IH ).EQ.ITORD ) THEN
C
C First do case where curl does not need to be taken
C
             IF ( ICURL.EQ.1 ) THEN
               SV( IND ) = SV( IND ) -
     1             FAC*QST( NOHARM, 3, IRN )/RDPLL1
             ENDIF
C
C Now do case where curl needs to be taken
C
             IF ( ICURL.EQ.2 ) THEN
               CALL QSTADF ( LH, NOHARM, NR, IRN, 2, QST, ORD,
     1                    D0F, D1F, D2F, D3F, D4F, RI, RO )
               SV( IND ) = SV( IND ) +
     1                     FAC*QST( NOHARM, 1, IRN )/RAD -
     2                     FAC*( D1F+D0F/RAD )/RDPLL1
             ENDIF
C
           ENDIF
c...........
        ENDDO
      ENDDO
C
C
C
      RETURN
      END
C*********************************************************************

