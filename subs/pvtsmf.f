C*********************************************************************
C subroutine Poloidal Velocity Time-Step Matrices Form ***************
C            -        -        -    -    -        -    ***************
C Steve Gibbons Wed Nov  1 15:04:22 MET 2000                         C
C____________________________________________________________________C
C                                                                    C
C If we are time-stepping, solving for a poloidal velocity, p^{i+1}, C
C such that                                                          C
C                                                                    C
C           [ p^{i+1} - p^i ]            [ (1-c) \nabla^2 p^{i+1} ]  C
C  c_e curl [ ------------- ] = c_i curl [    +  c \nabla^2 p^i   ]  C
C           [  delta t      ]    + Additional_forcing_terms          C
C                                                                    C
C then PVTSMF builds the two matrices, AM1 and BM1, such that        C
C                                                                    C
C AM1 p^{i+1} = BM1 p^{i} + Forcing terms.                           C
C                                                                    C
C The solution vector must consist only of poloidal velocity         C
C harmonics, of which there are NH1, and the values of these         C
C radial functions, at NR grid nodes, must be given by the index     C
C                                                                    C
C   ind = ( ih - 1 )*nr + ir                                         C
C                                                                    C
C where ih is the number of the radial function and ir is the        C
C number of the radial grid node.                                    C
C                                                                    C
C The radial functions are characterised by the integer arrays       C
C ML1, MM1 (not referred to here) and MP1.                           C
C                                                                    C
C ML1( ih ) gives the spherical harmonic degree, l.                  C
C MM1( ih ) gives the spherical harmonic order, m, when ih has a     C
C (cos m phi) dependence and -m for a (sin m phi) dependence.        C
C                                                                    C
C NBN is the number of bounding nodes.                               C
C Rows are added to the matrices between grid nodes 3 and NR - 2.    C
C The matrix AM1 is LU decomposed and so additional diagonal         C
C elements must be added using the routine AMSDEA. Due to the        C
C LU decomposition, AM1 must have the dimensions                     C
C  ( 3*NBN+1, NH1*NR )                                               C
C whereas BM1 is only required for a matrix-vector multiplication    C
C and so only the dimensions ( 2*NBN+1, NH1*NR ) are needed.         C
C                                                                    C
C It zeroes both matrices on input and assumes a LAPACK format for   C
C the banded matrix ( IMF = 1 in MATIND ).                           C
C                                                                    C
C The boundary conditions imposed in the different finite diff.      C
C schemes are given by the arrays MHIBC and MHOBC. Now MP1(ih) = is  C
C where MHIBC( is ) is one of the following:-                        C
C                                                                    C
C  MHIBC( is ) = 4 --> Function must vanish at the bndry AND         C
C                       first derivative must vanish at the bndry.   C
C  MHIBC( is ) = 5 --> Function must vanish at the bndry AND         C
C                       second derivative must vanish at the bndry.  C
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
C     NH1       : Number of poloidal velocity harmonics.             C
C     NDCS       : Number of distinct differencing coeff.s           C
C                  represented in SVFDC.                             C
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
C     ML1       : Dim ( NH1 ). See above.                            C
C     MP1       : Dim ( NH1 ). See above.                            C
C                                                                    C
C     IPIV1     : Dim (NH1*NR). Pivotting information for            C
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
C     AM1      : Matrix for solution of p^{i+1}. Has dimensions      C
C                 ( 3*NBN + 1, NH1*NR )                              C
C                                                                    C
C     BM1      : Matrix for multiplication of p^i. Has dimensions    C
C                 ( 2*NBN + 1, NH1*NR )                              C
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
      SUBROUTINE PVTSMF( NR, NH1, NDCS, NBN, MHIBC, MHOBC, NFDCM,
     1                   NDRVM, ML1, MP1, IPIV1, XARR, SVFDC,
     2                   AM1, BM1, CE, CI, CFAC, DELTAT )
      IMPLICIT NONE
C____________________________________________________________________C
C Variable declarations - Parameters ................................C
      INTEGER NR, NH1, NDCS, NBN, MHIBC( NDCS ), MHOBC( NDCS ),
     1        NFDCM, NDRVM, ML1( * ), MP1( * ), IPIV1( * )
      DOUBLE PRECISION XARR( * ), SVFDC( NFDCM, NR, NDRVM+1, NDCS ),
     1                AM1( 3*NBN+1, NH1*NR ), BM1( 2*NBN+1, NH1*NR ),
     2                CE, CI, CFAC, DELTAT 
C____________________________________________________________________C
C Variable declarations - Working variables .........................C
      INTEGER IOP, IMF, INARR( 3 ), ILNR, IRNR, INFO, N1, N2, KL,
     1        IH1, IHD, IPARS( 2 ), IS, KLE
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
      INARR( 3 ) = NH1
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
        PRINT *,' Subroutine PVTSMF.'
        PRINT *,' CFAC = ', CFAC
        PRINT *,' Program aborted.'
        STOP
      ENDIF
C
      IF ( DELTAT.EQ.ZERO .OR. CE.EQ.ZERO .OR. CI.EQ.ZERO ) THEN
        PRINT *,' Subroutine PVTSMF.'
        PRINT *,' CE = ', CE,' CI = ', CI,' DELTAT = ', DELTAT
        PRINT *,' Program aborted.'
        STOP
      ENDIF
C
C Zero the matrix AM1
C
      N2 = NH1*NR
      N1 = 3*NBN + 1
      CALL MATOP( AM1, ZERO, N1, N2, IOP )
C
      KLE  = KL
      ILNR = 3
      IRNR = NR - 2
C
C Add the curl of \nabla^2 parts onto AM1 matrix
C Curl of a poloidal harmonic is the -D_l^2 operator
C
      FAC  = CI*DELTAT*(1.0d0 - CFAC)
C
      IHD  = 4
      DO IH1 = 1, NH1
        IS         = MP1( IH1 )
        IPARS( 1 ) = 3
        IPARS( 2 ) = ML1( IH1 )
        CALL AMLICA( N1, N2, KL, KL, KLE, IMF, IH1, IH1, INARR,
     1               IHD, NBN, ILNR, IRNR, NFDCM, NR, NDCS, IS,
     2               NDRVM, NDRVM, IPARS, AMDLT, AM1, FAC, XARR,
     3               WORK, DPARS, SVFDC )
      ENDDO
C
C Add the curl parts of time derivative onto AM1 matrix
C Curl is -D_l for poloidal harmonics and so we multiply
C CE by -1.0d0 to make FAC
C
      FAC  = (-1.0d0)*CE
      IHD  = 2
      DO IH1 = 1, NH1
        IS         = MP1( IH1 )
        IPARS( 1 ) = 2
        IPARS( 2 ) = ML1( IH1 )
        CALL AMLICA( N1, N2, KL, KL, KLE, IMF, IH1, IH1, INARR,
     1               IHD, NBN, ILNR, IRNR, NFDCM, NR, NDCS, IS,
     2               NDRVM, NDRVM, IPARS, AMDLT, AM1, FAC, XARR,
     3               WORK, DPARS, SVFDC )
      ENDDO
C
C Add the diagonal elements to AM1 matrix.
C
      FAC = 1.0d0
      CALL AMSDEA( AM1, N1, N2, KL, KL, KLE, IMF, INARR,
     1             MP1, MHIBC, 'Inner Boundary', FAC, NDCS )
      CALL AMSDEA( AM1, N1, N2, KL, KL, KLE, IMF, INARR,
     1             MP1, MHOBC, 'Outer Boundary', FAC, NDCS )
C
C Now we must attempt to do an LU decomposition of the
C DMAT matrix ... for this, we need the LAPACK routine DGBTRF
C
      CALL DGBTRF( N2, N2, KL, KL, AM1, N1, IPIV1, INFO )
C
C Check for an error from LU decomp.
C
      IF ( INFO.NE.0 ) THEN
        PRINT *,' Subroutine PVTSMF.'
        PRINT *,' The LAPACK subroutine DGBTRF has been called'
        PRINT *,' and has returned ',INFO,' as a value of INFO. '
        PRINT *,' Program aborted.'
        STOP
      ENDIF
C
C Now zero the matrix BM1 ...
C
      N1 = 2*NBN + 1
      CALL MATOP( BM1, ZERO, N1, N2, IOP )
      KLE = 0
C
C Add the curl of \nabla^2 parts onto BM1 matrix.
C Curl of a poloidal harmonic is the -D_l^2 operator
C
      FAC  = CI*DELTAT*CFAC*(-1.0d0)
      IHD  = 4
      DO IH1 = 1, NH1
        IS         = MP1( IH1 )
        IPARS( 1 ) = 3
        IPARS( 2 ) = ML1( IH1 )
        CALL AMLICA( N1, N2, KL, KL, KLE, IMF, IH1, IH1, INARR,
     1               IHD, NBN, ILNR, IRNR, NFDCM, NR, NDCS, IS,
     2               NDRVM, NDRVM, IPARS, AMDLT, BM1, FAC, XARR,
     3               WORK, DPARS, SVFDC )
      ENDDO
C
C Add the curl parts of time derivative onto BM1 matrix
C Curl is -D_l for poloidal harmonics and so we multiply
C CE by -1.0d0 to make FAC
C
      FAC  = (-1.0d0)*CE
      IHD  = 2
      DO IH1 = 1, NH1
        IS         = MP1( IH1 )
        IPARS( 1 ) = 2
        IPARS( 2 ) = ML1( IH1 )
        CALL AMLICA( N1, N2, KL, KL, KLE, IMF, IH1, IH1, INARR,
     1               IHD, NBN, ILNR, IRNR, NFDCM, NR, NDCS, IS,
     2               NDRVM, NDRVM, IPARS, AMDLT, BM1, FAC, XARR,
     3               WORK, DPARS, SVFDC )
      ENDDO
C
      RETURN
      END
C*********************************************************************
