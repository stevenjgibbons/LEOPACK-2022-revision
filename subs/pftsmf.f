C*********************************************************************
C subroutine Poloidal Field Time-Step Matrices Form ******************
C            -        -     -    -    -        -    ******************
C Steve Gibbons Fri Nov  3 11:00:51 MET 2000                         C
C____________________________________________________________________C
C                                                                    C
C If we are time-stepping, solving for a poloidal field, p^{i+1},    C
C such that                                                          C
C                                                                    C
C       [ p^{i+1} - p^i ]         [ (1-c) \nabla^2 p^{i+1} ]         C
C  c_k  [ ------------- ] =  c_l  [    +  c \nabla^2 p^i   ]         C
C       [    delta t    ]      +    Additional_forcing_terms         C
C                                                                    C
C then PFTSMF builds the two matrices, AM4 and BM4, such that        C
C                                                                    C
C AM4 t^{i+1} = BM4 t^{i} + Forcing terms.                           C
C                                                                    C
C The solution vector must consist only of poloidal field            C
C harmonics, of which there are NH4, and the values of these         C
C radial functions, at NR grid nodes, must be given by the index     C
C                                                                    C
C   ind = ( ih - 1 )*nr + ir                                         C
C                                                                    C
C where ih is the number of the radial function and ir is the        C
C number of the radial grid node.                                    C
C                                                                    C
C The radial functions are characterised by the integer arrays       C
C ML4, MM4 (not referred to here) and MP4.                           C
C                                                                    C
C ML4( ih ) gives the spherical harmonic degree, l.                  C
C MM4( ih ) gives the spherical harmonic order, m, when ih has a     C
C (cos m phi) dependence and -m for a (sin m phi) dependence.        C
C                                                                    C
C NBN is the number of bounding nodes.                               C
C Rows are added to the matrices between grid nodes 2 and NR - 1.    C
C The matrix AM4 is LU decomposed and so additional diagonal         C
C elements must be added using the routine AMSDEA. Due to the        C
C LU decomposition, AM4 must have the dimensions                     C
C  ( 3*NBN+1, NH4*NR )                                               C
C whereas BM4 is only required for a matrix-vector multiplication    C
C and so only the dimensions ( 2*NBN+1, NH4*NR ) are needed.         C
C                                                                    C
C It zeroes both matrices on input and assumes a LAPACK format for   C
C the banded matrix ( IMF = 1 in MATIND ).                           C
C                                                                    C
C The boundary conditions imposed in the different finite diff.      C
C schemes are given by the arrays MHIBC and MHOBC. Now MP4(ih) = is  C
C where MHIBC( is ) is one of the following:-                        C
C                                                                    C
C    MHIBC( is ) = 7 --> insulating magnetic field.                  C
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
C     NH4       : Number of toroidal velocity harmonics.             C
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
C     ML4       : Dim ( NH4 ). See above.                            C
C     MP4       : Dim ( NH4 ). See above.                            C
C                                                                    C
C     IPIV4     : Dim (NH4*NR). Pivotting information for            C
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
C     AM4      : Matrix for solution of p^{i+1}. Has dimensions      C
C                 ( 3*NBN + 1, NH4*NR )                              C
C                                                                    C
C     BM4      : Matrix for multiplication of p^i. Has dimensions    C
C                 ( 2*NBN + 1, NH4*NR )                              C
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
      SUBROUTINE PFTSMF( NR, NH4, NDCS, NBN, MHIBC, MHOBC, NFDCM,
     1                   NDRVM, ML4, MP4, IPIV4, XARR, SVFDC,
     2                   AM4, BM4, CK, CL, CFAC, DELTAT )
      IMPLICIT NONE
C____________________________________________________________________C
C Variable declarations - Parameters ................................C
      INTEGER NR, NH4, NDCS, NBN, MHIBC( NDCS ), MHOBC( NDCS ),
     1        NFDCM, NDRVM, ML4( * ), MP4( * ), IPIV4( * )
      DOUBLE PRECISION XARR( * ), SVFDC( NFDCM, NR, NDRVM+1, NDCS ),
     1                AM4( 3*NBN+1, NH4*NR ), BM4( 2*NBN+1, NH4*NR ),
     2                CK, CL, CFAC, DELTAT 
C____________________________________________________________________C
C Variable declarations - Working variables .........................C
      INTEGER IOP, IMF, INARR( 3 ), ILNR, IRNR, INFO, N1, N2, KL,
     1        IH4, IHD, IPARS( 2 ), IS, KLE
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
      INARR( 3 ) = NH4
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
        PRINT *,' Subroutine PFTSMF.'
        PRINT *,' CFAC = ', CFAC
        PRINT *,' Program aborted.'
        STOP
      ENDIF
C
      IF ( DELTAT.EQ.ZERO .OR. CK.EQ.ZERO .OR. CL.EQ.ZERO ) THEN
        PRINT *,' Subroutine PFTSMF.'
        PRINT *,' CK = ', CK,' CL = ', CL,' DELTAT = ', DELTAT
        PRINT *,' Program aborted.'
        STOP
      ENDIF
C
C Zero the matrix AM4
C
      N2 = NH4*NR
      N1 = 3*NBN + 1
      CALL MATOP( AM4, ZERO, N1, N2, IOP )
C
      KLE  = KL
      ILNR = 2
      IRNR = NR - 1
C
C Add \nabla^2 parts onto AM4 matrix
C Lap of a poloidal field harmonic is the D_l operator
C
      FAC  = CL*DELTAT*(CFAC - 1.0d0)
C
      IHD  = 2
      DO IH4 = 1, NH4
        IS         = MP4( IH4 )
        IPARS( 1 ) = 2
        IPARS( 2 ) = ML4( IH4 )
        CALL AMLICA( N1, N2, KL, KL, KLE, IMF, IH4, IH4, INARR,
     1               IHD, NBN, ILNR, IRNR, NFDCM, NR, NDCS, IS,
     2               NDRVM, NDRVM, IPARS, AMDLT, AM4, FAC, XARR,
     3               WORK, DPARS, SVFDC )
      ENDDO
C
C Add the time derivative onto AM4 matrix
C
      FAC  = CK
      IHD  = 0
      DO IH4 = 1, NH4
        IS         = MP4( IH4 )
        IPARS( 1 ) = 1
        IPARS( 2 ) = ML4( IH4 )
        CALL AMLICA( N1, N2, KL, KL, KLE, IMF, IH4, IH4, INARR,
     1               IHD, NBN, ILNR, IRNR, NFDCM, NR, NDCS, IS,
     2               NDRVM, NDRVM, IPARS, AMDLT, AM4, FAC, XARR,
     3               WORK, DPARS, SVFDC )
      ENDDO
C
C Add the diagonal elements to AM4 matrix.
C
      FAC = 1.0d0
      CALL AMSDEA( AM4, N1, N2, KL, KL, KLE, IMF, INARR,
     1             MP4, MHIBC, 'Inner Boundary', FAC, NDCS )
      CALL AMSDEA( AM4, N1, N2, KL, KL, KLE, IMF, INARR,
     1             MP4, MHOBC, 'Outer Boundary', FAC, NDCS )
C
C Now we must attempt to do an LU decomposition of the
C AM4 matrix ... for this, we need the LAPACK routine DGBTRF
C
      CALL DGBTRF( N2, N2, KL, KL, AM4, N1, IPIV4, INFO )
C
C Check for an error from LU decomp.
C
      IF ( INFO.NE.0 ) THEN
        PRINT *,' Subroutine PFTSMF.'
        PRINT *,' The LAPACK subroutine DGBTRF has been called'
        PRINT *,' and has returned ',INFO,' as a value of INFO. '
        PRINT *,' Program aborted.'
        STOP
      ENDIF
C
C Now zero the matrix BM4 ...
C
      N1 = 2*NBN + 1
      CALL MATOP( BM4, ZERO, N1, N2, IOP )
      KLE = 0
C
C Add the \nabla^2 parts onto BM4 matrix.
C Lap of a poloidal field harmonic is the D_l operator
C
      FAC  = CL*DELTAT*CFAC
      IHD  = 2
      DO IH4 = 1, NH4
        IS         = MP4( IH4 )
        IPARS( 1 ) = 2
        IPARS( 2 ) = ML4( IH4 )
        CALL AMLICA( N1, N2, KL, KL, KLE, IMF, IH4, IH4, INARR,
     1               IHD, NBN, ILNR, IRNR, NFDCM, NR, NDCS, IS,
     2               NDRVM, NDRVM, IPARS, AMDLT, BM4, FAC, XARR,
     3               WORK, DPARS, SVFDC )
      ENDDO
C
C Add the time derivative onto BM4 matrix
C
      FAC  = CK
      IHD  = 0
      DO IH4 = 1, NH4
        IS         = MP4( IH4 )
        IPARS( 1 ) = 1
        IPARS( 2 ) = ML4( IH4 )
        CALL AMLICA( N1, N2, KL, KL, KLE, IMF, IH4, IH4, INARR,
     1               IHD, NBN, ILNR, IRNR, NFDCM, NR, NDCS, IS,
     2               NDRVM, NDRVM, IPARS, AMDLT, BM4, FAC, XARR,
     3               WORK, DPARS, SVFDC )
      ENDDO
C
      RETURN
      END
C*********************************************************************
