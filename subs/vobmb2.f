C*********************************************************************
C subroutine Vector Operations Banded Matrix Building Routine 2 ******
C            -      -          -      -      -        -       - ******
C Steve Gibbons Tue Nov 28 13:43:55 MET 2000                         C
C____________________________________________________________________C
C                                                                    C
C Forms banded matrices to perform operations on a vector of         C
C spherical harmonic radial functions. The operation required is     C
C specified by the integer flag IOPT:                                C
C                                                                    C
C Possible values for IOPT are                                       C
C ----------------------------                                       C
C                                                                    C
C  1: multiply p(r) to get Q(r)                                      C
C                                                                    C
C        now Q( r ) = L(L+1)/RAD p( r )                              C
C                                                                    C
C    IHD( highest derivative required ) = 0.                         C
C                                                                    C
C  2: multiply p(r) to get S(r)                                      C
C                                                                    C
C        now S( r ) = DSQRT( L(L+1) )( p/RAD + dp/dr )               C
C                                                                    C
C    IHD( highest derivative required ) = 1.                         C
C                                                                    C
C  3: multiply tau(r) to get T(r)                                    C
C                                                                    C
C        now T( r ) = -DSQRT( L(L+1) ) tau( r )                      C
C                                                                    C
C    IHD( highest derivative required ) = 0.                         C
C                                                                    C
C In the following q, s and t refer to the vectors                   C
C                                                                    C
C  4: curl of scaloidal vector Q(r) q                                C
C                                                                    C
C       curl[ Q(r) q ] = -DSQRT( L(L+1) ) Q( r )/RAD  t              C
C                                                                    C
C    IHD( highest derivative required ) = 0.                         C
C                                                                    C
C  5: curl of spheroidal vector S(r) s                               C
C                                                                    C
C       curl[ S(r) s ] = t[ dS/dr + S(r)/RAD ]                       C
C                                                                    C
C    IHD( highest derivative required ) = 1.                         C
C                                                                    C
C  6: curl of toroidal vector T(r) t: (scaloidal component)          C
C                                                                    C
C       curl[ T(r) t ] = - q[ T(r) DSQRT( L(L+1) )/RAD ]             C
C                        - s[ dT/dr + T(r)/RAD ]                     C
C                                                                    C
C    IHD( highest derivative required ) = 0.                         C
C                                                                    C
C  7: curl of toroidal vector T(r) t: (spheroidal component)         C
C                                                                    C
C       curl[ T(r) t ] = - q[ T(r) DSQRT( L(L+1) )/RAD ]             C
C                        - s[ dT/dr + T(r)/RAD ]                     C
C                                                                    C
C    IHD( highest derivative required ) = 1.                         C
C                                                                    C
C  8: Mulitply Q(r) to get p(r)                                      C
C                                                                    C
C       p( r ) = RAD/( L*L + L ) Q( r )                              C
C                                                                    C
C    IHD( highest derivative required ) = 0.                         C
C                                                                    C
C  9: Mulitply T(r) to get tau(r)                                    C
C                                                                    C
C       tau( r ) = - T( r ) / SQRT( L*L + L )                        C
C                                                                    C
C    IHD( highest derivative required ) = 0.                         C
C                                                                    C
C 10: Calculate a pure first derivative.                             C
C                                                                    C
C    IHD( highest derivative required ) = 1.                         C
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
C    IHD( highest derivative required ) = 0.                         C
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
C    IHD( highest derivative required ) = 1.                         C
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
C    IHD( highest derivative required ) = 0.                         C
C                                                                    C
C____________________________________________________________________C
C                                                                    C
C The matrix DPOM has dimensions ( N1, N2 ). The vector which is     C
C multiplied by DPOM contains NH radial functions, each with NR      C
C radial grid nodes. The radius at grid node IR is XARR( IR ).       C
C The IRth value of function IH is stored in element                 C
C                                                                    C
C   ( IH - 1 )*NR + IR                                               C
C                                                                    C
C The spherical harmonic degree of the radial function is stored     C
C in MLT( IH ).                                                      C
C                                                                    C
C N2 must be equal to NR*NH                                          C
C N1 must be equal to 1 + 2*K where K is the number of diagonals.    C
C                                                                    C
C If IOPT is such that IHD = 0, then K = 0.                          C
C If IOPT is such that IHD = 1, then K = NBN.                        C
C                                                                    C
C The matrix is zero-ed on entry.                                    C
C The terms are added to the matrix between nodes ILNR and IRNR.     C
C                                                                    C
C The finite difference coefficients are stored in the array FDCM    C
C which has dimensions ( NCFM, NR, 1 ) (NCFM=2*NBN) and this         C
C matrix is formed by the subroutine FDCMBD with the following       C
C parameters ...                                                     C
C                                                                    C
C                     NDRVM = 1                                      C
C                     NLMN  = ILNR                                   C
C                     NRMN  = IRNR                                   C
C                     NLMC  = NLMC                                   C
C                     NRMC  = NRMC                                   C
C                                                                    C
C____________________________________________________________________C
C                                                                    C
C Input variables :-                                                 C
C ===============                                                    C
C  Integer                                                           C
C  -------                                                           C
C                                                                    C
C     N1        : Leading dimension of DPOM.                         C
C     N2        : Second dimension of DPOM.                          C
C     K         : Number of diagonals in banded matrix.              C
C     NR        : Number of radial grid nodes.                       C
C     NH        : Number of radial functions.                        C
C     NBN       : Number of side nodes for radial derivatives.       C
C     NCFM      : Leading coefficient of FDCM.                       C
C     ILNR      : First grid node to operate upon.                   C
C     IRNR      : Last grid node to operate upon.                    C
C     MLT       : Dim (NH). Element IH contains degree L.            C
C     IOPT      : Choice of matrix operation. See above list.        C
C                                                                    C
C  Double Precision                                                  C
C  ----------------                                                  C
C                                                                    C
C     XARR      : Dim (NR). Radial node values.                      C
C     DPOM      : Dim (N1, N2). Double precision operations matrix.  C
C     FDCM      : Dim (NCFM,NR,1). Finite difference coefficients.   C
C                                                                    C
C____________________________________________________________________C
C
C*********************************************************************
      SUBROUTINE VOBMB2( N1, N2, K, NR, NH, NBN, NCFM, ILNR, IRNR,
     1                   MLT, IOPT, XARR, DPOM, FDCM, NLMC, NRMC )
      IMPLICIT NONE
C____________________________________________________________________C
C Variable declarations - Parameters ................................C
      INTEGER          N1, N2, K, NR, NH, NBN, NCFM, ILNR, IRNR,
     1                 MLT( NH ), IOPT, NLMC, NRMC
      DOUBLE PRECISION XARR( NR ), DPOM( N1, N2 ), FDCM(NCFM,NR,1)
C____________________________________________________________________C
C Variable declarations - Working variables .........................C
      INTEGER          KARR( 13 ), IHD, IOP, IPARS( 2 ), IMF,
     1                 INARR( 3 ), IH
      DOUBLE PRECISION CVEC( 2 ), ZERO, DPARS( 1 ), FAC
      PARAMETER        ( IOP = 0, ZERO = 0.0d0 )
      EXTERNAL         VOBMAR
C____________________________________________________________________C
C START OF PROGRAM **************************************************C
C____________________________________________________________________C
C
      KARR(  1 ) = 0
      KARR(  2 ) = NBN
      KARR(  3 ) = 0
      KARR(  4 ) = 0
      KARR(  5 ) = NBN
      KARR(  6 ) = 0
      KARR(  7 ) = NBN
      KARR(  8 ) = 0
      KARR(  9 ) = 0
      KARR( 10 ) = NBN
      KARR( 11 ) = 0
      KARR( 12 ) = NBN
      KARR( 13 ) = 0
C
      IF ( K.NE.KARR( IOPT ) ) THEN
        PRINT *,' Subroutine VOBMB2.'
        PRINT *,' IOPT = ', IOPT,' K = ',K
        PRINT *,' Program aborted.'
        STOP
      ENDIF
C
      IHD = 0
      IF ( K.EQ.NBN ) IHD = 1
C
      IF ( N1.NE.(2*K+1) .OR. N2.NE.(NR*NH) ) THEN
        PRINT *,' Subroutine VOBMB2.'
        PRINT *,' N1 = ', N1,' N2 = ', N2
        PRINT *,' NR = ', NR,' NH = ', NH
        PRINT *,' K  = ', K
        PRINT *,' Program aborted.'
        STOP
      ENDIF
C
C Clear the matrix
C
      CALL MATOP( DPOM, ZERO, N1, N2, IOP )
C
      INARR( 1 ) = 4
      INARR( 2 ) = NR
      INARR( 3 ) = NH
      FAC        = 1.0d0
      IMF        = 1
      IPARS( 1 ) = IOPT
      DO IH = 1, NH
        IPARS( 2 ) = MLT( IH )
        CALL NMLICA( N1, N2, K, K, IOP, IMF, IH, IH, INARR,
     1               IHD, K, ILNR, IRNR, NLMC, NRMC, NCFM,
     2               IMF, VOBMAR, DPOM, FAC, XARR, FDCM, DPARS,
     3               IPARS, NR, CVEC )
      ENDDO
C
      END
C*********************************************************************
