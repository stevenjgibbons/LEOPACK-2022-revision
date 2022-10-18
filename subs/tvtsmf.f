C*********************************************************************
C subroutine Toroidal Velocity Time-Step Matrices Form ***************
C            -        -        -    -    -        -    ***************
C Steve Gibbons Fri Nov  3 09:22:18 MET 2000                         C
C____________________________________________________________________C
C                                                                    C
C If we are time-stepping, solving for a toroidal velocity, t^{i+1}, C
C such that                                                          C
C                                                                    C
C           [ t^{i+1} - t^i ]            [ (1-c) \nabla^2 t^{i+1} ]  C
C  c_e curl [ ------------- ] = c_i curl [    +  c \nabla^2 t^i   ]  C
C           [  delta t      ]    + Additional_forcing_terms          C
C                                                                    C
C then TVTSMF builds the two matrices, AM2 and BM2, such that        C
C                                                                    C
C AM2 t^{i+1} = BM2 t^{i} + Forcing terms.                           C
C                                                                    C
C The solution vector must consist only of toroidal velocity         C
C harmonics, of which there are NH2, and the values of these         C
C radial functions, at NR grid nodes, must be given by the index     C
C                                                                    C
C   ind = ( ih - 1 )*nr + ir                                         C
C                                                                    C
C where ih is the number of the radial function and ir is the        C
C number of the radial grid node.                                    C
C                                                                    C
C The radial functions are characterised by the integer arrays       C
C ML2, MM2 (not referred to here) and MP2.                           C
C                                                                    C
C ML2( ih ) gives the spherical harmonic degree, l.                  C
C MM2( ih ) gives the spherical harmonic order, m, when ih has a     C
C (cos m phi) dependence and -m for a (sin m phi) dependence.        C
C                                                                    C
C NBN is the number of bounding nodes.                               C
C Rows are added to the matrices between grid nodes 2 and NR - 1.    C
C The matrix AM2 is LU decomposed and so additional diagonal         C
C elements must be added using the routine AMSDEA. Due to the        C
C LU decomposition, AM2 must have the dimensions                     C
C  ( 3*NBN+1, NH2*NR )                                               C
C whereas BM2 is only required for a matrix-vector multiplication    C
C and so only the dimensions ( 2*NBN+1, NH2*NR ) are needed.         C
C                                                                    C
C It zeroes both matrices on input and assumes a LAPACK format for   C
C the banded matrix ( IMF = 1 in MATIND ).                           C
C                                                                    C
C The boundary conditions imposed in the different finite diff.      C
C schemes are given by the arrays MHIBC and MHOBC. Now MP2(ih) = is  C
C where MHIBC( is ) is one of the following:-                        C
C                                                                    C
C  MHIBC( is ) = 2 --> Function must vanish at the bndry.            C
C  MHIBC( is ) = 6 --> r df/dr - f(r) = 0   at the bndry.            C
C                                                                    C
C Similarly for MHOBC. MHIBC and MHOBC need not be identical.        C
C____________________________________________________________________C
C                                                                    C
C Input variables :-                                                 C
C ===============                                                    C
C  Integer                                                           C
C  -------                                                           C
C                                                                    C
C     NR        : Number of grid nodes.                              C
C     NH2       : Number of toroidal velocity harmonics.             C
C     NDCS      : Number of distinct differencing coeff.s            C
C                 represented in SVFDC.                              C
C                                                                    C
C     NBN       : Maximum number of nodes on either side for         C
C                  central differencing.                             C
C                                                                    C
C     MHIBC     : Dimension ( NDCS ). See above.                     C
C     MHOBC     : Dimension ( NDCS ). See above.                     C
C     NFDCM     : Leading order of the array SVFDC.                  C
C                 This must be atleast (2*NBN + 1)                   C
C     NDRVM     : Maximum number of derivatives in SVFDC.            C
C                                                                    C
C     ML2       : Dim ( NH2 ). See above.                            C
C     MP2       : Dim ( NH2 ). See above.                            C
C                                                                    C
C     IPIV2     : Dim (NH2*NR). Pivotting information for            C
C                 LAPACK routine DGBTRF (LU decomposition).          C
C                                                                    C
C  Double Precision                                                  C
C  ----------------                                                  C
C                                                                    C
C     XARR      : Array of dimension (  NR  ).                       C
C                  XARR( i ) contains the value of x or r at the     C
C                   i^{th} radial grid node.                         C
C                                                                    C
C     SVFDC     : Finite difference coefficient matrix.              C
C                  Dimension ( NFDCM, NR, NDRVM+1, NDCS ).           C
C                                                                    C
C     AM2      : Matrix for solution of p^{i+1}. Has dimensions      C
C                 ( 3*NBN + 1, NH2*NR )                              C
C                                                                    C
C     BM2      : Matrix for multiplication of p^i. Has dimensions    C
C                 ( 2*NBN + 1, NH2*NR )                              C
C                                                                    C
C     CE        : Coefficient of time-derivative.                    C
C     CI        : Coefficient of viscous diffusion.                  C
C     CFAC      : Determines how implicit/explicit integration is.   C
C                 Must be strictly greater than 0.0 and strictly     C
C                 less than 1.0  The higher CFAC is, the more        C
C                 explicit the integration. CFAC = 0.5 --> Crank-    C
C                 Nicolson scheme.                                   C
C     DELTAT    : Time-step size.                                    C
C                                                                    C
C                                                                    C
C____________________________________________________________________C
C
C*********************************************************************
      SUBROUTINE TVTSMF( NR, NH2, NDCS, NBN, MHIBC, MHOBC, NFDCM,
     1                   NDRVM, ML2, MP2, IPIV2, XARR, SVFDC,
     2                   AM2, BM2, CE, CI, CFAC, DELTAT )
      IMPLICIT NONE
C____________________________________________________________________C
C Variable declarations - Parameters ................................C
      INTEGER NR, NH2, NDCS, NBN, MHIBC( NDCS ), MHOBC( NDCS ),
     1        NFDCM, NDRVM, ML2( * ), MP2( * ), IPIV2( * )
      DOUBLE PRECISION XARR( * ), SVFDC( NFDCM, NR, NDRVM+1, NDCS ),
     1                AM2( 3*NBN+1, NH2*NR ), BM2( 2*NBN+1, NH2*NR ),
     2                CE, CI, CFAC, DELTAT 
C____________________________________________________________________C
C Variable declarations - Working variables .........................C
      INTEGER IOP, IMF, INARR( 3 ), ILNR, IRNR, INFO, N1, N2, KL,
     1        IH2, IHD, IPARS( 2 ), IS, KLE
      DOUBLE PRECISION ZERO, ONE, FAC, DPARS( 1 ), WORK( 5 )
      PARAMETER ( ZERO = 0.0d0, ONE = 1.0d0, IOP = 0 )
      EXTERNAL AMDLT
C____________________________________________________________________C
C START OF PROGRAM **************************************************C
C____________________________________________________________________C
C Check input parameters
C
      INARR( 1 ) = 4
      INARR( 2 ) = NR
      INARR( 3 ) = NH2
      DPARS( 1 ) = ZERO
C (actually dpars is not referred to by AMDLT)
C
C Enforce matrix is built in LAPACK format
C
      IMF = 1
      KL  = NBN
C
C Check CFAC
C
      IF ( CFAC.LE.ZERO .OR. CFAC.GE.ONE ) THEN
        PRINT *,' Subroutine TVTSMF.'
        PRINT *,' CFAC = ', CFAC
        PRINT *,' Program aborted.'
        STOP
      ENDIF
C
      IF ( DELTAT.EQ.ZERO .OR. CE.EQ.ZERO .OR. CI.EQ.ZERO ) THEN
        PRINT *,' Subroutine TVTSMF.'
        PRINT *,' CE = ', CE,' CI = ', CI,' DELTAT = ', DELTAT
        PRINT *,' Program aborted.'
        STOP
      ENDIF
C
C Zero the matrix AM2
C
      N2 = NH2*NR
      N1 = 3*NBN + 1
      CALL MATOP( AM2, ZERO, N1, N2, IOP )
C
      KLE  = KL
      ILNR = 2
      IRNR = NR - 1
C
C Add the curl of \nabla^2 parts onto AM2 matrix
C Curl Lap of a toroidal vel. harmonic is the D_l operator
C
      FAC  = CI*DELTAT*(CFAC - 1.0d0)
C
      IHD  = 2
      DO IH2 = 1, NH2
        IS         = MP2( IH2 )
        IPARS( 1 ) = 2
        IPARS( 2 ) = ML2( IH2 )
        CALL AMLICA( N1, N2, KL, KL, KLE, IMF, IH2, IH2, INARR,
     1               IHD, NBN, ILNR, IRNR, NFDCM, NR, NDCS, IS,
     2               NDRVM, NDRVM, IPARS, AMDLT, AM2, FAC, XARR,
     3               WORK, DPARS, SVFDC )
      ENDDO
C
C Add the curl parts of time derivative onto AM2 matrix
C Curl is D_l^0(!) for toroidal harmonics and so we multiply
C CE by +1.0d0 to make FAC
C
      FAC  = CE
      IHD  = 0
      DO IH2 = 1, NH2
        IS         = MP2( IH2 )
        IPARS( 1 ) = 1
        IPARS( 2 ) = ML2( IH2 )
        CALL AMLICA( N1, N2, KL, KL, KLE, IMF, IH2, IH2, INARR,
     1               IHD, NBN, ILNR, IRNR, NFDCM, NR, NDCS, IS,
     2               NDRVM, NDRVM, IPARS, AMDLT, AM2, FAC, XARR,
     3               WORK, DPARS, SVFDC )
      ENDDO
C
C Add the diagonal elements to AM2 matrix.
C
      FAC = 1.0d0
      CALL AMSDEA( AM2, N1, N2, KL, KL, KLE, IMF, INARR,
     1             MP2, MHIBC, 'Inner Boundary', FAC, NDCS )
      CALL AMSDEA( AM2, N1, N2, KL, KL, KLE, IMF, INARR,
     1             MP2, MHOBC, 'Outer Boundary', FAC, NDCS )
C
C Now we must attempt to do an LU decomposition of the
C AM2 matrix ... for this, we need the LAPACK routine DGBTRF
C
      CALL DGBTRF( N2, N2, KL, KL, AM2, N1, IPIV2, INFO )
C
C Check for an error from LU decomp.
C
      IF ( INFO.NE.0 ) THEN
        PRINT *,' Subroutine TVTSMF.'
        PRINT *,' The LAPACK subroutine DGBTRF has been called'
        PRINT *,' and has returned ',INFO,' as a value of INFO. '
        PRINT *,' Program aborted.'
        STOP
      ENDIF
C
C Now zero the matrix BM2 ...
C
      N1 = 2*NBN + 1
      CALL MATOP( BM2, ZERO, N1, N2, IOP )
      KLE = 0
C
C Add the curl of \nabla^2 parts onto BM2 matrix.
C Curl Lap of a toroidal harmonic is the D_l operator
C
      FAC  = CI*DELTAT*CFAC
      IHD  = 2
      DO IH2 = 1, NH2
        IS         = MP2( IH2 )
        IPARS( 1 ) = 2
        IPARS( 2 ) = ML2( IH2 )
        CALL AMLICA( N1, N2, KL, KL, KLE, IMF, IH2, IH2, INARR,
     1               IHD, NBN, ILNR, IRNR, NFDCM, NR, NDCS, IS,
     2               NDRVM, NDRVM, IPARS, AMDLT, BM2, FAC, XARR,
     3               WORK, DPARS, SVFDC )
      ENDDO
C
C Add the curl parts of time derivative onto BM2 matrix
C Curl is D_l^0(!) for toroidal harmonics and so we multiply
C CE by +1.0d0 to make FAC
C
      FAC  = CE
      IHD  = 0
      DO IH2 = 1, NH2
        IS         = MP2( IH2 )
        IPARS( 1 ) = 1
        IPARS( 2 ) = ML2( IH2 )
        CALL AMLICA( N1, N2, KL, KL, KLE, IMF, IH2, IH2, INARR,
     1               IHD, NBN, ILNR, IRNR, NFDCM, NR, NDCS, IS,
     2               NDRVM, NDRVM, IPARS, AMDLT, BM2, FAC, XARR,
     3               WORK, DPARS, SVFDC )
      ENDDO
C
      RETURN
      END
C*********************************************************************
