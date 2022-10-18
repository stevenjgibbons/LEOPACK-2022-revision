C*********************************************************************
C subroutine Solution Vector Kumar Roberts Velocity Functions ********
C            -        -      -     -       -        -         ********
C Steve Gibbons Wed Oct  6 10:20:50 BST 1999                         C
C____________________________________________________________________C
C                                                                    C
C Enters the radial functions to the Kumar Roberts velocity into     C
C the appropriate locations of the solution vector, V.               C
C____________________________________________________________________C
C                                                                    C
C Input variables :-                                                 C
C ===============                                                    C
C  Integer                                                           C
C  -------                                                           C
C                                                                    C
C     NR        : Number of radial grid nodes.                       C
C     INARR     : Int. parameter array corresponding to V.           C
C                 Dim ( * ) but the following must be true.          C
C                                                                    C
C                 INARR( 1 ) = IFORMF                                C
C                 INARR( 2 ) = NRR     See INDFUN for details        C
C                 INARR( 3 ) = NH      nrr must equal nr             C
C                                                                    C
C     MHT       : Array length ( * ) - atleast length NH             C
C                                                                    C
C     MHT( ih ) = 1 if 'ih' is poloidal velocity harmonic.           C
C     MHT( ih ) = 2 if 'ih' is toroidal velocity harmonic.           C
C     MHT( ih ) = 3 if 'ih' is temperature harmonic.                 C
C     MHT( ih ) = 4 if 'ih' is poloidal magnetic field harmonic.     C
C     MHT( ih ) = 5 if 'ih' is toroidal magnetic field harmonic.     C
C                                                                    C
C     MHL       : Spherical harmonic degree, l, of harmonic 'ih'.    C
C                                                                    C
C     MHM       : Spherical harmonic degree, m, of harmonic 'ih'     C
C                if harmonic has (cos m phi dependence) - otherwise  C
C                MHM( ih ) = -m                                      C
C                                                                    C
C     ILNR      : Lowest radial node to be used.                     C
C     IRNR      : Highest radial node to be used.                    C
C                                                                    C
C  Double Precision                                                  C
C  ----------------                                                  C
C                                                                    C
C     V         : DP vector of dimension ( * ) containing solution   C
C                 Length must be atleast NH*NR                       C
C                                                                    C
C     DPARR     : Array for use by RADVLF                            C
C                 DPARR( 1 ) MUST contain RI and                     C
C                 DPARR( 2 ) MUST contain RO.                        C
C                                                                    C
C     XARR      : Array ( * ). Only referenced if IFORMF = 3 or 4.   C
C                 In this case dim( xarr ) must be atleast NR and    C
C                 xarr( IR ) will contain the x value at node IR.    C
C                                                                    C
C     EPS0      : Velocity parameter for t_1^0 ( r )                 C
C     EPS1      : Velocity parameter for s_2^0 ( r )                 C
C     EPS2      : Velocity parameter for s_2^2s( r )                 C
C     EPS3      : Velocity parameter for s_2^2c( r )                 C
C     LAMBDA    : Part of t_1^0 ( r ) - generally set to zero.       C
C     PPAR      : Parameter P in the cyclonic velocity functions.    C
C                  Sin / Cos ( PPAR * pi * r ) .....                 C
C                                                                    C
C____________________________________________________________________C
C
C*********************************************************************
      SUBROUTINE SVKRVF( NR, INARR, MHT, MHL, MHM, ILNR, IRNR, V,
     1                   DPARR, XARR, EPS0, EPS1, EPS2, EPS3, LAMBDA,
     2                   PPAR )
      IMPLICIT NONE
C____________________________________________________________________C
C Variable declarations - Parameters ................................C
      INTEGER NR, INARR( * ), MHT( * ), MHL( * ), MHM( * ), ILNR,
     1        IRNR
      DOUBLE PRECISION V( * ), DPARR( * ), XARR( NR ), 
     1                 EPS0, EPS1, EPS2, EPS3, LAMBDA, PPAR
C____________________________________________________________________C
C Variable declarations - Working variables .........................C
      DOUBLE PRECISION R, RAD, RI, RO, FAC, TOL, PI, H, VEL, PPI, RR
      INTEGER ICOUNT, IR, NH, NRR, IH, IND, INDFUN
      LOGICAL SWITC0, SWITC1, SWITC2, SWITC3
      PARAMETER ( TOL = 1.0d-7, PI=3.14159265358979312D0 )
C____________________________________________________________________C
C START OF PROGRAM **************************************************C
C____________________________________________________________________C
C
      PPI = PPAR * PI
C     .
      NRR = INARR( 2 )
      NH  = INARR( 3 )
C     .
      IF ( NRR.NE.NR ) THEN
        PRINT *,' Subroutine SVKRVF.'
        PRINT *,' NR         = ', NR
        PRINT *,' INARR( 2 ) = ', INARR( 2 )
        PRINT *,' Program aborted.'
        STOP
      ENDIF
C     .
      RI = DPARR( 1 )
      RO = DPARR( 2 )
C     .
      IF ( (RO-RI).LT.TOL ) THEN
        PRINT *,' Subroutine SVKRVF.'
        PRINT *,' RI = ', RI
        PRINT *,' RO = ', RO
        PRINT *,' Program aborted.'
        STOP
      ENDIF
C     .
      FAC = 1.0d0/(RO-RI)
C     .
      SWITC0 = .TRUE.
      SWITC1 = .TRUE.
      SWITC2 = .TRUE.
      SWITC3 = .TRUE.
C     .
      ICOUNT = 0
C     .
      DO IH = 1, NH
C       .
C       . First, do toroidal velocity function t_1^0 ( r )
C       . With our treatment of the radial function,
C       . t_1^0 ( r ) = r ( 1 - r*r ) + LAMBDA * r
C       .
        IF (     MHT( IH ).EQ.2 .AND. MHL( IH ).EQ.1 .AND.
     1           MHM( IH ).EQ.0 .AND. SWITC0 ) THEN
          DO IR = ILNR, IRNR
            IND = INDFUN( IR, IH, INARR )
            CALL RADVLF( RAD, IR, INARR, DPARR, XARR, H )
            R = FAC*( RAD - RI )
            RR = R*R
            VEL = R * ( 1.0d0 - RR ) + LAMBDA * R
            V( IND ) = VEL*EPS0
          ENDDO
          SWITC0 = .FALSE.
          ICOUNT = ICOUNT + 1
        ENDIF
C       .
C       . Now, do poloidal velocity function s_2^0 ( r )
C       . With our treatment of the radial function,
C       . s_2^0 ( r ) = r^5 ( 1 - r^2 )^3             
C       .
        IF (     MHT( IH ).EQ.1 .AND. MHL( IH ).EQ.2 .AND.
     1           MHM( IH ).EQ.0 .AND. SWITC1 ) THEN
          DO IR = ILNR, IRNR
            IND = INDFUN( IR, IH, INARR )
            CALL RADVLF( RAD, IR, INARR, DPARR, XARR, H )
            R = FAC*( RAD - RI )
            RR = R*R
            VEL = R*RR*RR*( 1.0d0 - RR )**3
            V( IND ) = VEL*EPS1
          ENDDO
          SWITC1 = .FALSE.
          ICOUNT = ICOUNT + 1
        ENDIF
C       .
C       . Now, do poloidal velocity function s_2^2s ( r )
C       . With our treatment of the radial function,
C       . s_2^2s ( r ) = r^3.(1-r^2)^2.cos( PPI*r )
C       .
        IF (     MHT( IH ).EQ.1 .AND. MHL( IH ).EQ.2 .AND.
     1           MHM( IH ).EQ.(-2) .AND. SWITC2 ) THEN
          DO IR = ILNR, IRNR
            IND = INDFUN( IR, IH, INARR )
            CALL RADVLF( RAD, IR, INARR, DPARR, XARR, H )
            R = FAC*( RAD - RI )
            RR = R*R
            VEL = R*RR*(1.0d0-RR)*(1.0d0-RR)*COS( PPI*R )
            V( IND ) = VEL*EPS2
          ENDDO
          SWITC2 = .FALSE.
          ICOUNT = ICOUNT + 1
        ENDIF
C       .
C       . Now, do poloidal velocity function s_2^2c ( r )
C       . With our treatment of the radial function,
C       . s_2^2c ( r ) = r^3.(1-r^2)^2.sin( PPI*r )
C       .
        IF (     MHT( IH ).EQ.1 .AND. MHL( IH ).EQ.2 .AND.
     1           MHM( IH ).EQ.2 .AND. SWITC3 ) THEN
          DO IR = ILNR, IRNR
            IND = INDFUN( IR, IH, INARR )
            CALL RADVLF( RAD, IR, INARR, DPARR, XARR, H )
            R = FAC*( RAD - RI )
            RR = R*R
            VEL = R*RR*(1.0d0-RR)*(1.0d0-RR)*SIN( PPI*R )
            V( IND ) = VEL*EPS3
          ENDDO
          SWITC3 = .FALSE.
          ICOUNT = ICOUNT + 1
        ENDIF
C       .
      ENDDO
C     .
      IF ( ICOUNT.NE.4 ) THEN
        PRINT *,' Subroutine SVKRVF. ICOUNT = ', ICOUNT
        PRINT *,' Error in input harmonics '
        IF ( SWITC0 ) WRITE ( 6, 900 )
        IF ( SWITC1 ) WRITE ( 6, 901 )
        IF ( SWITC2 ) WRITE ( 6, 902 )
        IF ( SWITC3 ) WRITE ( 6, 903 )
        PRINT *,' Program aborted.'
        STOP
      ENDIF
 900  FORMAT(' Harmonic t_1^0 not found. ')
 901  FORMAT(' Harmonic s_2^0 not found. ')
 902  FORMAT(' Harmonic s_2^2s not found. ')
 903  FORMAT(' Harmonic s_2^2c not found. ')
C     .
      RETURN
      END
C*********************************************************************
