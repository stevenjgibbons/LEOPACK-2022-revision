C*********************************************************************
C subroutine Adapted Vector operation Banded Matrix Building Routine *
C            -       -                       -      -        -       *
C Steve Gibbons Tue Nov 28 14:54:39 MET 2000                         C
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
C Unlike AVBMBR, this routine uses the SVFDC finite derivative       C
C array and so will in general assume a particular kind of boundary  C
C condition, or will require values from points 1 to NR in order to  C
C take derivatives.                                                  C
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
C in MLT( IH ). The finite difference scheme required for this       C
C radial function is given by MPT( IH )                              C
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
C     NCFM      : Leading coefficient of SVFDC.                      C
C     NDRVM     : Maximum number of derivatives stored in SVFDC.     C
C     ILNR      : First grid node to operate upon.                   C
C     IRNR      : Last grid node to operate upon.                    C
C     MLT       : Dim (NH). Element IH contains degree L.            C
C     MPT       : Dim (NH). Element IH contains degree IS.           C
C     IOPT      : Choice of matrix operation. See above list.        C
C     NDCS      : Number of finite difference schemes.               C
C                                                                    C
C  Double Precision                                                  C
C  ----------------                                                  C
C                                                                    C
C     XARR      : Dim (NR). Radial node values.                      C
C     DPOM      : Dim (N1, N2). Double precision operations matrix.  C
C                                                                    C
C     SVFDC     : Finite difference coefficient matrix.              C
C                  Dimension ( NCFM, NR, NDRVM+1, NDCS ).            C
C                                                                    C
C____________________________________________________________________C
C
C*********************************************************************
      SUBROUTINE AVBMBR( N1, N2, K, NR, NH, NBN, NCFM, NDRVM, ILNR,
     1             IRNR, MLT, MPT, IOPT, NDCS, XARR, DPOM, SVFDC )
      IMPLICIT NONE
C____________________________________________________________________C
C Variable declarations - Parameters ................................C
      INTEGER          N1, N2, K, NR, NH, NBN, NCFM, NDRVM, ILNR,
     1                 IRNR, MLT( NH ), MPT( NH ), IOPT, NDCS
      DOUBLE PRECISION XARR( NR ), DPOM( N1, N2 ),
     1                 SVFDC( NCFM, NR, NDRVM+1, NDCS )
C____________________________________________________________________C
C Variable declarations - Working variables .........................C
      INTEGER          IARR( 13 ), IHD, IOP, IPARS( 2 ), IMF,
     1                 INARR( 3 ), IH, IS
      DOUBLE PRECISION CVEC( 2 ), ZERO, DPARS( 1 ), FAC
      PARAMETER        ( IOP = 0, ZERO = 0.0d0 )
      EXTERNAL         VOBMAR
C____________________________________________________________________C
C START OF PROGRAM **************************************************C
C____________________________________________________________________C
C
      IARR(  1 ) = 0
      IARR(  2 ) = 1
      IARR(  3 ) = 0
      IARR(  4 ) = 0
      IARR(  5 ) = 1
      IARR(  6 ) = 0
      IARR(  7 ) = 1
      IARR(  8 ) = 0
      IARR(  9 ) = 0
      IARR( 10 ) = 1
      IARR( 11 ) = 0
      IARR( 12 ) = 1
      IARR( 13 ) = 0
C
      IHD = IARR( IOPT )
C
      IF ( N1.NE.(2*K+1) .OR. N2.NE.(NR*NH) ) THEN
        PRINT *,' Subroutine AVBMBR.'
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
        IS         = MPT( IH )
        CALL AMLICA( N1, N2, K, K, IOP, IMF, IH, IH, INARR,
     1               IHD, NBN, ILNR, IRNR, NCFM, NR, NDCS, IS,
     2               NDRVM, NDRVM, IPARS, VOBMAR, DPOM, FAC,
     3               XARR, CVEC, DPARS, SVFDC )
      ENDDO
C
      END
C*********************************************************************
