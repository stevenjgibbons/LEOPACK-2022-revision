C*********************************************************************
C subroutine Solution Vector Kumar Roberts Velocity Functions ********
C            -        -      -     -       -        -         ********
C Steve Gibbons Tue Mar 13 12:26:55 MET 2001                         C
C____________________________________________________________________C
C                                                                    C
C Enters the radial functions to the Dudley James velocity into      C
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
C     XARR      : Array ( * ). Only referenced if IFORMF = 3 or 4.   C
C                 In this case dim( xarr ) must be atleast NR and    C
C                 xarr( IR ) will contain the x value at node IR.    C
C                                                                    C
C     EPS       : Velocity parameter for s_2^0 ( r )                 C
C     LAMBDA    : Part of t_1^0 ( r ) - generally set to zero.       C
C                                                                    C
C____________________________________________________________________C
C
C*********************************************************************
      SUBROUTINE SVDJVF( NR, INARR, MHT, MHL, MHM, ILNR, IRNR, V,
     1                   XARR, EPS, LAMBDA )
      IMPLICIT NONE
C____________________________________________________________________C
C Variable declarations - Parameters ................................C
      INTEGER NR, INARR( * ), MHT( * ), MHL( * ), MHM( * ), ILNR,
     1        IRNR
      DOUBLE PRECISION V( * ), XARR( NR ), EPS, LAMBDA
C____________________________________________________________________C
C Variable declarations - Working variables .........................C
      DOUBLE PRECISION R, RAD, RI, RO, FAC, TOL, PI, VEL, RPI
      INTEGER          ICOUNT, IR, NH, NRR, IH, IND, INDFUN
      LOGICAL          SWITC0, SWITC1
      PARAMETER ( TOL = 1.0d-7, PI=3.14159265358979312D0 )
C____________________________________________________________________C
C START OF PROGRAM **************************************************C
C____________________________________________________________________C
C
      NRR = INARR( 2 )
      NH  = INARR( 3 )
C     .
      IF ( NRR.NE.NR ) THEN
        PRINT *,' Subroutine SVDJVF.'
        PRINT *,' NR         = ', NR
        PRINT *,' INARR( 2 ) = ', INARR( 2 )
        PRINT *,' Program aborted.'
        STOP
      ENDIF
C     .
      RI = XARR(  1 )
      RO = XARR( NR )
C     .
      IF ( (RO-RI).LT.TOL ) THEN
        PRINT *,' Subroutine SVDJVF.'
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
C     .
      ICOUNT = 0
C     .
      DO IH = 1, NH
C       .
C       . First, do toroidal velocity function t_1^0 ( r )
C       . With our treatment of the radial function,
C       . t_1^0 ( r ) = sin ( pi*r ) + LAMBDA * r
C       .
        IF (     MHT( IH ).EQ.2 .AND. MHL( IH ).EQ.1 .AND.
     1           MHM( IH ).EQ.0 .AND. SWITC0 ) THEN
          DO IR = ILNR, IRNR
            IND = INDFUN( IR, IH, INARR )
            RAD = XARR( IR )
            R   = FAC*( RAD - RI )
            RPI = PI*R
            VEL = DSIN( RPI ) + LAMBDA*R
            V( IND ) = VEL
          ENDDO
          SWITC0 = .FALSE.
          ICOUNT = ICOUNT + 1
        ENDIF
C       .
C       . Now, do poloidal velocity function s_2^0 ( r )
C       . With our treatment of the radial function,
C       . s_2^0 ( r ) = r*sin( pi*r )
C       .
        IF (     MHT( IH ).EQ.1 .AND. MHL( IH ).EQ.2 .AND.
     1           MHM( IH ).EQ.0 .AND. SWITC1 ) THEN
          DO IR = ILNR, IRNR
            IND = INDFUN( IR, IH, INARR )
            RAD = XARR( IR )
            R   = FAC*( RAD - RI )
            RPI = PI*R
            VEL = R*DSIN( RPI )
            V( IND ) = VEL*EPS
          ENDDO
          SWITC1 = .FALSE.
          ICOUNT = ICOUNT + 1
        ENDIF
C       .
      ENDDO
C     .
      IF ( ICOUNT.NE.2 ) THEN
        PRINT *,' Subroutine SVDJVF. ICOUNT = ', ICOUNT
        PRINT *,' Error in input harmonics '
        IF ( SWITC0 ) WRITE ( 6, 900 )
        IF ( SWITC1 ) WRITE ( 6, 901 )
        PRINT *,' Program aborted.'
        STOP
      ENDIF
 900  FORMAT(' Harmonic t_1^0 not found. ')
 901  FORMAT(' Harmonic s_2^0 not found. ')
C     .
      RETURN
      END
C*********************************************************************
