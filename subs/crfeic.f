C*********************************************************************
C subroutine Chebyshev Radial Function Energy Integral Calculate *****
C            -         -      -        -      -        -         *****
C Steve Gibbons Sat May  6 09:48:22 BST 2000                         C
C____________________________________________________________________C
C                                                                    C
C Returns 0.5d0 * \int_{volume} A.A dV where A is a single vector    C
C harmonic - the type being indicated by MHT( ih ).                  C
C                                                                    C
C This can be either magnetic energy or kinetic energy.              C
C                                                                    C
C The value for A = ( A_r, A_theta, A_phi ) is returned in ENINT.    C
C                                                                    C
C If MHT( ih ) does not correspond to a poloidal or toroidal         C
C vector harmonic, ENINT is returned zero.                           C
C                                                                    C
C____________________________________________________________________C
C Input variables :-                                                 C
C ===============                                                    C
C  Double Precision                                                  C
C  ----------------                                                  C
C     ENINT	: The integral.					     C
C     DLOW	: The accepted difference between successive         C
C                  integrations. Suggested value about 1.0d-7        C
C     CCV  	: Chebyshev coefficient vector. Dim (*)              C
C                  with dimension ( * ).                             C
C  Integer                                                           C
C  -------                                                           C
C                                                                    C
C     INARR     : Array of integers. Dimension ( * )                 C
C                 See CHINDF.                                        C
C                                                                    C
C     IH        : Number of harmonic to be evaluated.                C
C                                                                    C
C     MHT       : Array length ( * ) - atleast length NH             C
C                                                                    C
C MHT defines what each scalar function in a solution vector         C
C represents.                                                        C
C                                                                    C
C         MHT( IH ) = 1 --> harmonic is poloidal velocity            C
C         MHT( IH ) = 2 --> harmonic is toroidal velocity            C
C         MHT( IH ) = 3 --> harmonic is temperature.                 C
C         MHT( IH ) = 4 --> harmonic is poloidal magnetic field      C
C         MHT( IH ) = 5 --> harmonic is toroidal magnetic field      C
C                                                                    C
C     MHL       : Array length ( * ) - atleast length NH             C
C                  Sph. harm. degree, l.                             C
C____________________________________________________________________C
C Calls the trapezoidal routine CRFETR                               C
C____________________________________________________________________C
C
C*********************************************************************
      SUBROUTINE CRFEIC ( ENINT, DLOW, INARR, IH, CCV, MHT, MHL )
      IMPLICIT NONE
C____________________________________________________________________C
C Variable declarations - Parameters ................................C
      DOUBLE PRECISION X1, X2, ENINT, DLOW, CCV( * )
      INTEGER INARR( * ), IH, MHT( * ), MHL( * )
C____________________________________________________________________C
C Variable declarations - Working variables .........................C
      INTEGER J, JMAX, IPARR( 2 ), NH, K, N, INDCCV, CHINDF,
     1        IFLAG, L, ITYPE
      DOUBLE PRECISION ST, OS, OST, FAC, SQRLL1, PI, DK
      PARAMETER ( JMAX = 20, PI = 3.14159265358979312D0 )
C____________________________________________________________________C
C START OF PROGRAM **************************************************C
C____________________________________________________________________C
C Check validity of arguments ....
C
      IF ( DLOW.LT.1.0d-20 .OR. DLOW.GT.1.0d-4 ) THEN
         PRINT *,' Subroutine CRFEIC.'
         PRINT *,' DLOW = ',DLOW,'. Should be approx. 1.0d-7 '
         PRINT *,' Program aborted.'
         STOP
      ENDIF
C
      NH = INARR( 3 )
      IF ( IH.LT.1 .OR. IH.GT.NH ) THEN
        PRINT *,'  Subroutine CRFEIC.'
        PRINT *,'  IH = ', IH,' NH = ', NH
        PRINT *,'  Program aborted.'
        STOP
      ENDIF
C
      L     = MHL( IH )
      ITYPE = MHT( IH )
      IFLAG = 0
      ENINT = 0.0d0
      DK    = SQRLL1( L )
      FAC   = 2.0d0*PI*DK*DK/(2.0d0*DBLE( L ) + 1.0d0)
      IF ( ITYPE.EQ.1 .OR. ITYPE.EQ.4 ) IFLAG = 1
      IF ( ITYPE.EQ.2 .OR. ITYPE.EQ.5 ) IFLAG = 2
      IF ( IFLAG.EQ.0 ) RETURN
C
      IPARR( 1 ) = IFLAG
      IPARR( 2 ) = L
C
      N             = INARR( 2 ) - 2
      K             = N
      INDCCV        = CHINDF( K, IH, INARR )
      X1            = CCV( INDCCV )
C
      K             = N+1
      INDCCV        = CHINDF( K, IH, INARR )
      X2            = CCV( INDCCV )
C
      OST = -1.0d30
      OS = -1.0d30
C
      DO J = 1, JMAX
        CALL CRFETR( X1, X2, ST, J, INARR, IPARR, IH, CCV )
        ENINT = (4.0d0*ST-OST)/3.0d0
        IF ( ABS( ENINT - OS ).LT.DLOW*ABS(OS) ) THEN
          ENINT = ENINT*FAC
          RETURN
        ENDIF
        OS = ENINT
        OST = ST
      ENDDO
C ........ integration hasn't worked; set output to zero.
      ENINT = 0.0d0
C
      RETURN
      END
C*********************************************************************
