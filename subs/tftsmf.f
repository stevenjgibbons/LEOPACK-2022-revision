C*********************************************************************
C subroutine Toroidal Field Time-Step Matrices Form ******************
C            -        -     -    -    -        -    ******************
C Steve Gibbons Fri Nov  3 11:00:51 MET 2000                         C
C____________________________________________________________________C
C                                                                    C
C If we are time-stepping, solving for a toroidal field, t^{i+1},    C
C such that                                                          C
C                                                                    C
C       [ t^{i+1} - t^i ]         [ (1-c) \nabla^2 t^{i+1} ]         C
C  c_k  [ ------------- ] =  c_l  [    +  c \nabla^2 t^i   ]         C
C       [    delta t    ]      +    Additional_forcing_terms         C
C                                                                    C
C then TFTSMF builds the two matrices, AM5 and BM5, such that        C
C                                                                    C
C AM5 t^{i+1} = BM5 t^{i} + Forcing terms.                           C
C                                                                    C
C The solution vector must consist only of toroidal field            C
C harmonics, of which there are NH5, and the values of these         C
C radial functions, at NR grid nodes, must be given by the index     C
C                                                                    C
C   ind = ( ih - 1 )*nr + ir                                         C
C                                                                    C
C where ih is the number of the radial function and ir is the        C
C number of the radial grid node.                                    C
C                                                                    C
C The radial functions are characterised by the integer arrays       C
C ML5, MM5 (not referred to here) and MP5.                           C
C                                                                    C
C ML5( ih ) gives the spherical harmonic degree, l.                  C
C MM5( ih ) gives the spherical harmonic order, m, when ih has a     C
C (cos m phi) dependence and -m for a (sin m phi) dependence.        C
C                                                                    C
C NBN is the number of bounding nodes.                               C
C Rows are added to the matrices between grid nodes 2 and NR - 1.    C
C The matrix AM5 is LU decomposed and so additional diagonal         C
C elements must be added using the routine AMSDEA. Due to the        C
C LU decomposition, AM5 must have the dimensions                     C
C  ( 3*NBN+1, NH5*NR )                                               C
C whereas BM5 is only required for a matrix-vector multiplication    C
C and so only the dimensions ( 2*NBN+1, NH5*NR ) are needed.         C
C                                                                    C
C It zeroes both matrices on input and assumes a LAPACK format for   C
C the banded matrix ( IMF = 1 in MATIND ).                           C
C                                                                    C
C The boundary conditions imposed in the different finite diff.      C
C schemes are given by the arrays MHIBC and MHOBC. Now MP5(ih) = is  C
C where MHIBC( is ) is one of the following:-                        C
C                                                                    C
C  MHIBC( is ) = 2 --> Function must vanish at the bndry.            C
C  MHIBC( is ) = 3 --> First derivative must vanish at the bndry.    C
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
C     NH5       : Number of toroidal velocity harmonics.             C
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
C     ML5       : Dim ( NH5 ). See above.                            C
C     MP5       : Dim ( NH5 ). See above.                            C
C                                                                    C
C     IPIV5     : Dim (NH5*NR). Pivotting information for            C
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
C     AM5      : Matrix for solution of p^{i+1}. Has dimensions      C
C                 ( 3*NBN + 1, NH5*NR )                              C
C                                                                    C
C     BM5      : Matrix for multiplication of p^i. Has dimensions    C
C                 ( 2*NBN + 1, NH5*NR )                              C
C                                                                    C
C     CK        : Coefficient of time-derivative.                    C
C     CL        : Coefficient of magnetic diffusion.                 C
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
      SUBROUTINE TFTSMF( NR, NH5, NDCS, NBN, MHIBC, MHOBC, NFDCM,
     1                   NDRVM, ML5, MP5, IPIV5, XARR, SVFDC,
     2                   AM5, BM5, CK, CL, CFAC, DELTAT )
      IMPLICIT NONE
C____________________________________________________________________C
C Variable declarations - Parameters ................................C
      INTEGER NR, NH5, NDCS, NBN, MHIBC( NDCS ), MHOBC( NDCS ),
     1        NFDCM, NDRVM, ML5( * ), MP5( * ), IPIV5( * )
      DOUBLE PRECISION XARR( * ), SVFDC( NFDCM, NR, NDRVM+1, NDCS ),
     1                AM5( 3*NBN+1, NH5*NR ), BM5( 2*NBN+1, NH5*NR ),
     2                CK, CL, CFAC, DELTAT 
C____________________________________________________________________C
C Variable declarations - Working variables .........................C
      INTEGER IOP, IMF, INARR( 3 ), ILNR, IRNR, INFO, N1, N2, KL,
     1        IH5, IHD, IPARS( 2 ), IS, KLE
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
      INARR( 3 ) = NH5
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
        PRINT *,' Subroutine TFTSMF.'
        PRINT *,' CFAC = ', CFAC
        PRINT *,' Program aborted.'
        STOP
      ENDIF
C
      IF ( DELTAT.EQ.ZERO .OR. CK.EQ.ZERO .OR. CL.EQ.ZERO ) THEN
        PRINT *,' Subroutine TFTSMF.'
        PRINT *,' CK = ', CK,' CL = ', CL,' DELTAT = ', DELTAT
        PRINT *,' Program aborted.'
        STOP
      ENDIF
C
C Zero the matrix AM5
C
      N2 = NH5*NR
      N1 = 3*NBN + 1
      CALL MATOP( AM5, ZERO, N1, N2, IOP )
C
      KLE  = KL
      ILNR = 2
      IRNR = NR - 1
C
C Add \nabla^2 parts onto AM5 matrix
C Lap of a toroidal field harmonic is the D_l operator
C
      FAC  = CL*DELTAT*(CFAC - 1.0d0)
C
      IHD  = 2
      DO IH5 = 1, NH5
        IS         = MP5( IH5 )
        IPARS( 1 ) = 2
        IPARS( 2 ) = ML5( IH5 )
        CALL AMLICA( N1, N2, KL, KL, KLE, IMF, IH5, IH5, INARR,
     1               IHD, NBN, ILNR, IRNR, NFDCM, NR, NDCS, IS,
     2               NDRVM, NDRVM, IPARS, AMDLT, AM5, FAC, XARR,
     3               WORK, DPARS, SVFDC )
      ENDDO
C
C Add the time derivative onto AM5 matrix
C
      FAC  = CK
      IHD  = 0
      DO IH5 = 1, NH5
        IS         = MP5( IH5 )
        IPARS( 1 ) = 1
        IPARS( 2 ) = ML5( IH5 )
        CALL AMLICA( N1, N2, KL, KL, KLE, IMF, IH5, IH5, INARR,
     1               IHD, NBN, ILNR, IRNR, NFDCM, NR, NDCS, IS,
     2               NDRVM, NDRVM, IPARS, AMDLT, AM5, FAC, XARR,
     3               WORK, DPARS, SVFDC )
      ENDDO
C
C Add the diagonal elements to AM5 matrix.
C
      FAC = 1.0d0
      CALL AMSDEA( AM5, N1, N2, KL, KL, KLE, IMF, INARR,
     1             MP5, MHIBC, 'Inner Boundary', FAC, NDCS )
      CALL AMSDEA( AM5, N1, N2, KL, KL, KLE, IMF, INARR,
     1             MP5, MHOBC, 'Outer Boundary', FAC, NDCS )
C
C Now we must attempt to do an LU decomposition of the
C AM5 matrix ... for this, we need the LAPACK routine DGBTRF
C
      CALL DGBTRF( N2, N2, KL, KL, AM5, N1, IPIV5, INFO )
C
C Check for an error from LU decomp.
C
      IF ( INFO.NE.0 ) THEN
        PRINT *,' Subroutine TFTSMF.'
        PRINT *,' The LAPACK subroutine DGBTRF has been called'
        PRINT *,' and has returned ',INFO,' as a value of INFO. '
        PRINT *,' Program aborted.'
        STOP
      ENDIF
C
C Now zero the matrix BM5 ...
C
      N1 = 2*NBN + 1
      CALL MATOP( BM5, ZERO, N1, N2, IOP )
      KLE = 0
C
C Add the \nabla^2 parts onto BM5 matrix.
C Lap of a toroidal field harmonic is the D_l operator
C
      FAC  = CL*DELTAT*CFAC
      IHD  = 2
      DO IH5 = 1, NH5
        IS         = MP5( IH5 )
        IPARS( 1 ) = 2
        IPARS( 2 ) = ML5( IH5 )
        CALL AMLICA( N1, N2, KL, KL, KLE, IMF, IH5, IH5, INARR,
     1               IHD, NBN, ILNR, IRNR, NFDCM, NR, NDCS, IS,
     2               NDRVM, NDRVM, IPARS, AMDLT, BM5, FAC, XARR,
     3               WORK, DPARS, SVFDC )
      ENDDO
C
C Add the time derivative onto BM5 matrix
C
      FAC  = CK
      IHD  = 0
      DO IH5 = 1, NH5
        IS         = MP5( IH5 )
        IPARS( 1 ) = 1
        IPARS( 2 ) = ML5( IH5 )
        CALL AMLICA( N1, N2, KL, KL, KLE, IMF, IH5, IH5, INARR,
     1               IHD, NBN, ILNR, IRNR, NFDCM, NR, NDCS, IS,
     2               NDRVM, NDRVM, IPARS, AMDLT, BM5, FAC, XARR,
     3               WORK, DPARS, SVFDC )
      ENDDO
C
      RETURN
      END
C*********************************************************************
