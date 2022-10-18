C*********************************************************************
C subroutine Scalar Function 2 Spectral Decomposition add xsv 1 ******
C            -      -        - -        -                     - ******
C Steve Gibbons Tue Feb 13 10:43:55 WET 2001                         C
C____________________________________________________________________C
C                                                                    C
C Pre-calculates the array SF2SA( NH3, NTHP ) which is used by       C
C scalar real --> spectral transform routine  SF2SD2.                C
C                                                                    C
C____________________________________________________________________C
C                                                                    C
C Input variables :-                                                 C
C ===============                                                    C
C                                                                    C
C  Integer                                                           C
C  -------                                                           C
C                                                                    C
C     LH        : Maximum spherical harmonic degree, l.              C
C     NH3       : Number of scalar functions                         C
C     ML3       : Array dim ( NH3 ). Sph. harm degree, L.            C
C     MM3       : Array dim ( NH3 ). Sph. harm order, M, or -M.      C
C     M0        : Smallest non-zero wavenumber in solution.          C
C     NTHP      : Number of theta points.                            C
C     NPHP      : Number of phi points.                              C
C                                                                    C
C  Double Precision                                                  C
C  ----------------                                                  C
C                                                                    C
C     PA        : Schmidt Normalised Legendre Functions              C
C                  Dim { ( LH + 1 )*( LH + 2 )/2 , NTHP }            C
C                  P_l^m(X_j) = PA( L*(L+1)/2 + M + 1 , j )          C
C                                                                    C
C     GAUW      : Gauss weights as calculated by GAUWTS. ( NTHP )    C
C                                                                    C
C     FAC       : Coefficient of function to be added.               C
C                                                                    C
C     SF2SA     : Dim ( NH3, NTHP ). Array for use by SF2SD2.        C
C____________________________________________________________________C
C
C*********************************************************************
      SUBROUTINE SF2SD1( LH, NH3, ML3, MM3, M0, NTHP, NPHP, PA, GAUW,
     1                   FAC, SF2SA )
      IMPLICIT NONE
C____________________________________________________________________C
C Variable declarations - Parameters ................................C
      INTEGER          LH, NH3, ML3( NH3 ), MM3( NH3 ), M0, NTHP, NPHP
      DOUBLE PRECISION PA( ( LH + 1 )*( LH + 2 )/2, NTHP ),
     1                 GAUW( NTHP ), FAC, SF2SA( NH3, NTHP )
C____________________________________________________________________C
C Variable declarations - Working variables .........................C
      INTEGER          ITHETA, IH3, ICS, L, M, IP, MPS
      DOUBLE PRECISION FACL, DALF, WEIGHT
C____________________________________________________________________C
C START OF PROGRAM **************************************************C
C____________________________________________________________________C
C
C ............... begin looping around theta points ..................
      DO ITHETA = 1, NTHP
        WEIGHT = GAUW( ITHETA )
C       .
C       . Loop around harmonics
C       .
        DO IH3 = 1, NH3
C
           L      = ML3( IH3 )
C
C First do monopole term
C
           IF ( L.EQ.0 ) THEN
             SF2SA( IH3, ITHETA ) = 0.5d0*FAC*WEIGHT/DBLE( NPHP )
             GOTO 70
           ENDIF
C
           FACL   = 0.25d0*DBLE( 2*L + 1 )
           M      = MM3( IH3 )
           ICS    = 1
C          .
C          . Modify M and ICS if the harmonic is sin (m phi) dep.
C          .
           IF ( M.LT.0   ) ICS = 2
           IF ( ICS.EQ.2 ) M   = -M
C          .
C          . Store associated Legendre Functions
C          .
           IP     = L*(L+1)/2+M+1
           DALF   = PA( IP, ITHETA )
C          .
C          . Make sure that harmonic is
C          . compatible with the symmetry imposition
C          . MPS (pseudo M) = M/M0
C          .
           MPS = M/M0
           IF ( MPS*M0.NE.M ) THEN
             PRINT *,' Subroutine SF2SD1.'
             PRINT *,' MM3(',IH3,') = ', MM3( IH3 )
             PRINT *,' Program aborted.'
             STOP
           ENDIF
C          .
           SF2SA( IH3, ITHETA ) = FAC*WEIGHT*DALF*FACL
C          .
 70     CONTINUE
        ENDDO
C       .      Ended loop ih3 = 1, nh3
C       .
      ENDDO
C ............. ended looping around theta points ...................
      RETURN
      END
C*********************************************************************
