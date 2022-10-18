C*********************************************************************
C subroutine Vector Operations Banded Matrix Auxiliary Routine *******
C            -      -          -      -      -         -       *******
C Steve Gibbons Mon Nov 27 15:20:35 WET 2000                         C
C____________________________________________________________________C
C                                                                    C
C We wish to form banded matrices with which to perform many         C
C different operations. This subroutine may be called as an          C
C EXTERNAL routine from subroutines such as AMLICA and NMLICA        C
C in order to form these matrices.                                   C
C                                                                    C
C DPARS is not referred to.                                          C
C                                                                    C
C RAD is the radius.                                                 C
C                                                                    C
C IPARS( 1 ) contains the integer IOPT which specifies which         C
C matrix we are building. This may take any of the values listed     C
C below.                                                             C
C                                                                    C
C IPARS( 2 ) contains the spherical harmonic degree, l.              C
C                                                                    C
C IHD is dependent upon IOPT (see below)                             C
C                                                                    C
C Possible values for IOPT are                                       C
C ----------------------------                                       C
C                                                                    C
C  1: multiply p(r) to get Q(r)                                      C
C                                                                    C
C        now Q( r ) = L(L+1)/RAD p( r )                              C
C                                                                    C
C        and so   IHD = 0, CVEC( 1 ) = L(L+1)/RAD                    C
C                 *******  **********************                    C
C                                                                    C
C  2: multiply p(r) to get S(r)                                      C
C                                                                    C
C        now S( r ) = DSQRT( L(L+1) )( p/RAD + dp/dr )               C
C                                                                    C
C        and so   IHD = 1, CVEC( 1 ) = DSQRT( L(L+1) )/RAD           C
C                 *******  CVEC( 2 ) = DSQRT( L(L+1) )               C
C                          ***************************               C
C                                                                    C
C  3: multiply tau(r) to get T(r)                                    C
C                                                                    C
C        now T( r ) = -DSQRT( L(L+1) ) tau( r )                      C
C                                                                    C
C        and so   IHD = 0, CVEC( 1 ) = -DSQRT( L(L+1) )              C
C                 *******  ****************************              C
C                                                                    C
C In the following q, s and t refer to the vectors                   C
C                                                                    C
C  4: curl of scaloidal vector Q(r) q                                C
C                                                                    C
C       curl[ Q(r) q ] = -DSQRT( L(L+1) ) Q( r )/RAD  t              C
C                                                                    C
C        and so   IHD = 0, CVEC( 1 ) = -DSQRT( L(L+1) )/RAD          C
C                 *******  ********************************          C
C                                                                    C
C  5: curl of spheroidal vector S(r) s                               C
C                                                                    C
C       curl[ S(r) s ] = t[ dS/dr + S(r)/RAD ]                       C
C                                                                    C
C        and so   IHD = 1, CVEC( 1 ) = 1/RAD                         C
C                 *******  CVEC( 2 ) = 1.0d0                         C
C                          *****************                         C
C                                                                    C
C  6: curl of toroidal vector T(r) t: (scaloidal component)          C
C                                                                    C
C       curl[ T(r) t ] = - q[ T(r) DSQRT( L(L+1) )/RAD ]             C
C                        - s[ dT/dr + T(r)/RAD ]                     C
C                                                                    C
C        and so   IHD = 0, CVEC( 1 ) = -DSQRT( L(L+1) )/RAD          C
C                 *******  ********************************          C
C                                                                    C
C  7: curl of toroidal vector T(r) t: (spheroidal component)         C
C                                                                    C
C       curl[ T(r) t ] = - q[ T(r) DSQRT( L(L+1) )/RAD ]             C
C                        - s[ dT/dr + T(r)/RAD ]                     C
C                                                                    C
C        and so   IHD = 1, CVEC( 1 ) = - 1.0d0 / RAD                 C
C                 *******  CVEC( 2 ) = - 1.0d0                       C
C                          *************************                 C
C                                                                    C
C  8: Mulitply Q(r) to get p(r)                                      C
C                                                                    C
C       p( r ) = RAD/( L*L + L ) Q( r )                              C
C                                                                    C
C        and so   IHD = 0, CVEC( 1 ) = RAD/( L*L + L )               C
C                 *******  ***************************               C
C                                                                    C
C  9: Mulitply T(r) to get tau(r)                                    C
C                                                                    C
C       tau( r ) = - T( r ) / SQRT( L*L + L )                        C
C                                                                    C
C        and so   IHD = 0, CVEC( 1 ) = - 1 / SQRT( L*L + L )         C
C                 *******  *********************************         C
C                                                                    C
C 10: Calculate a pure first derivative.                             C
C                                                                    C
C        and so   IHD = 1, CVEC( 1 ) = 0                             C
C                 *******  CVEC( 2 ) = 1.0d0                         C
C                          *****************                         C
C                                                                    C
C 11: Curl of scaloidal vector Q(r) q to give tau( r ) radial        C
C      function                                                      C
C                                                                    C
C        curl[ Q(r) q ] = T_c( r ) t                                 C
C                                                                    C
C     where T_c( r ) = -DSQRT( L(L+1) ) Q( r )/RAD                   C
C                                                                    C
C     and so tau( r ) = Q( r )/RAD                                   C
C                                                                    C
C        and so   IHD = 0, CVEC( 1 ) = 1.0d0 / RAD                   C
C                 *******  ***********************                   C
C                                                                    C
C 12: Curl of spheroidal vector S(r) s to give tau( r ) radial       C
C      function                                                      C
C                                                                    C
C        curl[ S(r) s ] = T_c( r ) t                                 C
C                                                                    C
C     where T_c( r ) = dS/dr + S(r)/RAD                              C
C                                                                    C
C       and so tau( r ) =  -( dS/dr + S(r)/RAD )/DSQRT( L(L+1) )     C
C                                                                    C
C        and so   IHD = 1, CVEC( 1 ) = -1.0d0/(RAD*DSQRT( L(L+1) ))  C
C                 *******  CVEC( 2 ) = -1.0d0/DSQRT( L(L+1) )        C
C                          ****************************************  C
C                                                                    C
C 13: Curl of toroidal vector T( r ) t to give p( r )                C
C                                                                    C
C        curl[ T(r) t ] = Q_c( r ) q   +   S_c( r ) s                C
C                                                                    C
C         where  Q_c( r ) = -T(r) DSQRT( L(L+1) )/RAD   and          C
C                S_c( r ) = -[ dT/dr + T(r)/RAD ]                    C
C                                                                    C
C   hence p( r ) = -T(r)/DSQRT( L(L+1) )                             C
C                                                                    C
C        and so   IHD = 0, CVEC( 1 ) = -1.0d0/DSQRT( L(L+1) )        C
C                 *******  **********************************        C
C____________________________________________________________________C
C
C*********************************************************************
      SUBROUTINE VOBMAR( CVEC, RAD, IPARS, DPARS, IHD )
      IMPLICIT NONE
C____________________________________________________________________C
C Variable declarations - Parameters ................................C
      INTEGER IPARS( * ), IHD
      DOUBLE PRECISION CVEC( * ), RAD, DPARS( * )
C____________________________________________________________________C
C Variable declarations - Working variables .........................C
      INTEGER IHDMIN( 13 ), IOPT, L, ND
      DOUBLE PRECISION TOL, DLFAC, DLL1, DSLL1
      PARAMETER ( TOL = 1.0d-8 )
C____________________________________________________________________C
C START OF PROGRAM **************************************************C
C____________________________________________________________________C
C Check on value of RAD
C     .
      IF ( DABS( RAD ).LT.TOL ) THEN
        PRINT *,' Subroutine VOBMAR.'
        PRINT *,' RAD = ', RAD
        PRINT *,' Program aborted.'
        STOP
      ENDIF
C     .
      IOPT = IPARS( 1 )
      L    = IPARS( 2 )
C     .
      IF ( IOPT.LT.1 .OR. IOPT.GT.13 ) THEN
        PRINT *,' Subroutine VOBMAR.'
        PRINT *,' IOPT = ', IOPT
        PRINT *,' Program aborted.'
        STOP
      ENDIF
C     .
      DLFAC = DBLE( L )
      DLL1  = DLFAC*DLFAC + DLFAC
      DSLL1 = DSQRT( DLL1 )
C     .
      IHDMIN(  1 ) = 0
      IHDMIN(  2 ) = 1
      IHDMIN(  3 ) = 0
      IHDMIN(  4 ) = 0
      IHDMIN(  5 ) = 1
      IHDMIN(  6 ) = 0
      IHDMIN(  7 ) = 1
      IHDMIN(  8 ) = 0
      IHDMIN(  9 ) = 0
      IHDMIN( 10 ) = 1
      IHDMIN( 11 ) = 0
      IHDMIN( 12 ) = 1
      IHDMIN( 13 ) = 0
C     .
      IF ( IHD.NE.IHDMIN( IOPT ) ) THEN
        PRINT *,' Subroutine VOBMAR.'
        PRINT *,' IHD = ', IOPT
        PRINT *,' Should be ',IHDMIN( IOPT ),' for option.'
        PRINT *,' Program aborted.'
        STOP
      ENDIF
C     .
      DO ND = 1, IHD + 1
        CVEC( ND ) = 0.0d0
      ENDDO
C     .
      IF ( IOPT.EQ.1 ) THEN
        CVEC( 1 ) = DLL1/RAD
        RETURN
      ENDIF
C     .
      IF ( IOPT.EQ.2 ) THEN
        CVEC( 1 ) = DSLL1/RAD
        CVEC( 2 ) = DSLL1
        RETURN
      ENDIF
C     .
      IF ( IOPT.EQ.3 ) THEN
        CVEC( 1 ) = (-1.0d0)*DSLL1
        RETURN
      ENDIF
C     .
      IF ( IOPT.EQ.4 ) THEN
        CVEC( 1 ) = (-1.0d0)*DSLL1/RAD
        RETURN
      ENDIF
C     .
      IF ( IOPT.EQ.5 ) THEN
        CVEC( 1 ) = 1.0d0/RAD
        CVEC( 2 ) = 1.0d0
        RETURN
      ENDIF
C     .
      IF ( IOPT.EQ.6 ) THEN
        CVEC( 1 ) = (-1.0d0)*DSLL1/RAD
        RETURN
      ENDIF
C     .
      IF ( IOPT.EQ.7 ) THEN
        CVEC( 1 ) = (-1.0d0)/RAD
        CVEC( 2 ) = (-1.0d0)
        RETURN
      ENDIF
C     .
      IF ( IOPT.EQ.8 ) THEN
        CVEC( 1 ) = RAD/DLL1
        RETURN
      ENDIF
C     .
      IF ( IOPT.EQ.9 ) THEN
        CVEC( 1 ) = (-1.0d0)/DSLL1
        RETURN
      ENDIF
C     .
      IF ( IOPT.EQ.10 ) THEN
        CVEC( 1 ) = 0.0d0
        CVEC( 2 ) = 1.0d0
        RETURN
      ENDIF
C     .
      IF ( IOPT.EQ.11 ) THEN
        CVEC( 1 ) = 1.0d0/RAD
        RETURN
      ENDIF
C     .
      IF ( IOPT.EQ.12 ) THEN
        CVEC( 1 ) = -1.0d0/(DSLL1*RAD)
        CVEC( 2 ) = -1.0d0/DSLL1
        RETURN
      ENDIF
C     .
      IF ( IOPT.EQ.13 ) THEN
        CVEC( 1 ) = -1.0d0/DSLL1
        RETURN
      ENDIF
C     .
      END
C*********************************************************************
