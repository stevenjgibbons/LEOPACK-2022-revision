C*********************************************************************
C subroutine Double Precision Vector Drifting Frame Time derivative **
C            -      -         -      -        -     -               **
C Steve Gibbons 22.3.99                                              C
C____________________________________________________________________C
C                                                                    C
C Adds to RVEC, the time derivative of the vector OLDSV for the case C
C where the system drifts with constant rate C.                      C
C                                                                    C
C Define PHI = phi - c t      and so                                 C
C                                                                    C
C d PHI / dt = - c                                                   C
C                                                                    C
C____________________________________________________________________C
C Input variables :-                                                 C
C ===============                                                    C
C  Integer                                                           C
C  -------                                                           C
C     NDIM      : Length of the solution vector.                     C
C     NR        : Number of radial grid nodes.                       C
C     NH        : Number of harmonics (all types)                    C
C                                                                    C
C     MHT       : Dimension ( NH )                                   C
C     MHL       : Dimension ( NH )                                   C
C                                                                    C
C  mht, mhl, mhm and mhc define the spherical harmonics present      C
C  in the solution vector. For spherical harmonic number I;          C
C                                                                    C
C   MHT( I ) = 1 for a poloidal velocity vector                      C
C   MHT( I ) = 2 for a toroidal velocity vector                      C
C   MHT( I ) = 3 for a temperature / codensity term                  C
C                                                                    C
C   MHL( I ) = spherical harmonic degree, l                          C
C                                                                    C
C   MHM( I ) = spherical harmonic order, m                           C
C                                                                    C
C   MHC( I ) = 1 for a cosine dependence in phi and                  C
C            = 2  "  "  sine     "        "  "                       C
C                                                                    C
C  Double Precision                                                  C
C  ----------------                                                  C
C     RI        : Radius of inner boundary.                          C
C     RO        : Radius of outer boundary.                          C
C     RVEC      : DP vector of dimension ( NDIM )                    C
C     OLDSV     : DP vector of dimension ( NDIM ) - old vector       C
C     FACT      : Multiplication factor for temperature term         C
C     FACV      : Multiplication factor for velocity term            C
C     CVAL      : Value of drift rate.                               C
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
      SUBROUTINE DPVDFT ( NDIM, NR, NH, MHT, MHL, MHM, MHC, RI, RO,
     1                    RVEC, OLDSV, FACT, FACV, ORD, CVAL )
      IMPLICIT NONE
C____________________________________________________________________C
C Variable declarations - Parameters ................................C
      INTEGER NDIM, NR, NH,
     1        MHT( NH ), MHL( NH ), MHM( NH ), MHC( NH )
      DOUBLE PRECISION RI, RO, RVEC( NDIM ), FACT, FACV,
     1                 OLDSV( NDIM ), CVAL
      CHARACTER *(2) ORD
C____________________________________________________________________C
C Variable declarations - Working variables .........................C
      DOUBLE PRECISION H, RAD, DL, FAC, TOL, DPCM,
     1                 D0F, D1F, D2F, D3F, D4F
      INTEGER IIH, INDI, L, IRN, IOH, M, ICS, IT, JH, INDO
      PARAMETER ( TOL = 1.0d-6 )
C____________________________________________________________________C
C START OF PROGRAM **************************************************C
C____________________________________________________________________C
C
C Check the length of the vector
C
      IF ( NDIM.NE.NR*NH ) THEN
         PRINT *,' Subroutine DPVDFT. Bad dimensions '
         PRINT *,' Array length = ',NDIM
         PRINT *,' NR = ', NR,' NH = ',NH
         STOP
      ENDIF
C
C Quick exit for zero CVAL
C
      IF ( ABS( CVAL ).LT.TOL ) RETURN
C
      FAC = 0.0d0
      H = ( RO - RI ) / DBLE( NR - 1 )
C
c     . Loop around the harmonics
c
      DO IIH = 1, NH
        IT  = MHT( IIH )
        L   = MHL( IIH )
        M   = MHM( IIH )
        ICS = MHC( IIH )
        IF ( M.EQ.0 ) GOTO 500
        IF ( IT.EQ.1 .AND. ABS( FACV ).LT.TOL ) GOTO 500
        IF ( IT.EQ.2 .AND. ABS( FACV ).LT.TOL ) GOTO 500
        IF ( IT.EQ.3 .AND. ABS( FACT ).LT.TOL ) GOTO 500
C
C Ok - so it seems this harmonic needs doing ...
C Loop around the harmonics to find it's 'partner'
C
        DO JH = 1, NH
           IF (     MHT( JH ).EQ.IT      .AND.
     1              MHL( JH ).EQ.L       .AND.
     2              MHM( JH ).EQ.M       .AND.
     3              MHC( JH ).NE.ICS    ) THEN
             IOH = JH
             GOTO 499
           ENDIF 
        ENDDO 
 499    CONTINUE
C
C Calculate the additional factor M*C with minus
C sign if necessary
C
        IF ( ICS.EQ.1 ) THEN
          DPCM = DBLE( M )*CVAL
        ENDIF
        IF ( ICS.EQ.2 ) THEN
          DPCM = DBLE( M )*CVAL*(-1.0d0)
        ENDIF
C
C First do poloidal velocity terms
        IF ( MHT( IIH ).EQ.1 ) THEN
c
           DO IRN = 1, NR
             INDI = ( IRN - 1 )*NH + IIH
             INDO = ( IRN - 1 )*NH + IOH
             RAD = RI + H*DBLE( IRN - 1 )
             CALL SVDERF ( NDIM, NH, IIH, NR, IRN, OLDSV, ORD,
     1                    D0F, D1F, D2F, D3F, D4F, RI, RO )
             FAC  = DL( L, RAD, D0F, D1F, D2F )*(-1.0d0)
             RVEC( INDO ) = RVEC( INDO ) + FAC*FACV*DPCM
           ENDDO
C
        ENDIF
C Now do toroidal velocity terms
        IF ( MHT( IIH ).EQ.2 ) THEN
c
           DO IRN = 1, NR
             INDI = ( IRN - 1 )*NH + IIH
             INDO = ( IRN - 1 )*NH + IOH
             RVEC( INDO ) = RVEC( INDO ) +
     1                       FACV*DPCM*OLDSV( INDI )
           ENDDO
C
        ENDIF
C Now do temperature terms
        IF ( MHT( IIH ).EQ.3 ) THEN
c
           DO IRN = 1, NR
             INDI = ( IRN - 1 )*NH + IIH
             INDO = ( IRN - 1 )*NH + IOH
             RVEC( INDO ) = RVEC( INDO ) +
     1                       FACT*DPCM*OLDSV( INDI )
           ENDDO
C
        ENDIF
 500  CONTINUE
      ENDDO
C
      RETURN
      END
C*********************************************************************
