C*********************************************************************
C subroutine Magnetic Field Diffusion Matrix ForM ********************
C            -        -     -         -      -  - ********************
C Steve Gibbons Tue Nov  2 08:44:11 GMT 1999                         C
C____________________________________________________________________C
C                                                                    C
C If we are time stepping the induction equation, this routine       C
C gives us the matrix which will solve                               C
C                                                                    C
C ( 1 + dt (c-1) Lap ) B^{i+1} = ( 1 + dt c Lap ) B^{i} + Forc. Term C
C ^^^^^^^^^^^^^^^^^^^^                                               C
C It zeros the matrix, assumes that all boundary conditions are      C
C insulating, assumes all harmonics from 1 to NH are magnetic field  C
C and assumes a LAPACK format of the banded matrix (IMF = 1 in       C
C MATIND).                                                           C
C                                                                    C
C For a harmonic ih, MHT( ih ) must be 4 for poloidal and 5 for      C
C toroidal. MHL( ih ) gives the spherical harmonic degree, l.        C
C MHM( ih ) gives the spherical harmonic order, m, when ih has a     C
C (cos m phi) dependence and -m for a (sin m phi) dependence.        C
C                                                                    C
C MHP( ih ) = is where MHBCS( is ) = 2 when MHT( ih ) = 5 and        C
C MHBCS( is ) = 7 when MHT( ih ) = 4.                                C
C                                                                    C
C The magnetic field solution vector must be stored with IFORMF = 4. C
C (See indfun) i.e. INDFUN = ( IH - 1 )*NR + IR                      C
C                                                                    C
C____________________________________________________________________C
C                                                                    C
C Input variables :-                                                 C
C ===============                                                    C
C  Integer                                                           C
C  -------                                                           C
C                                                                    C
C     KL        : Number of lower diagonals in matrix.               C
C                                                                    C
C     NDCS      : Maximum distinct finite difference schemes         C
C                  stored by the array SVFDC.                        C
C                                                                    C
C     NBN       : Number of bounding nodes.                          C
C                                                                    C
C     NBAND     : Leading dimension of matrix, DMAT.                 C
C                 Must be (3*KL + 1)                                 C
C                                                                    C
C     N2        : Second dimension of matrix DMAT.                   C
C                 Must be NH*NR                                      C
C                                                                    C
C     NR        : Number of radial grid nodes.                       C
C     NH        : Number of spherical harmonic radial functions.     C
C                                                                    C
C     NDRVS     : Number of highest derivative for which             C
C                  coefficients are stored by the array SVFDC.       C
C                                                                    C
C     NDRVM     : Limit on NDRVS. Array bound for SVFDC.             C
C                                                                    C
C     NFDCM     : Leading dimension of FDCM. At least (2*NBN+1)      C
C                                                                    C
C     IPIV      : Working array for LAPACK routine DGBTRF. Dim (N2). C
C                                                                    C
C     MHBCS     : Dim ( NDCS ). See above.                           C
C     MHT       : Dim ( NH ). See above.                             C
C     MHL       : Dim ( NH ). See above.                             C
C     MHM       : Dim ( NH ). See above.                             C
C     MHP       : Dim ( NH ). See above.                             C
C                                                                    C
C  Double Precision                                                  C
C  ----------------                                                  C
C     DMAT      : Diffusion matrix. Dim ( NBAND, N2 ).               C
C                                                                    C
C     DC        : Diffusion coefficient. Must be positive.           C
C                                                                    C
C     DELTAT    : Time step.                                         C
C                                                                    C
C     CFAC      : Determines how explicit/implicit time stepping is. C
C                 (Corresponds to 'c' in above equation).            C
C                                                                    C
C     SVFDC     : Finite difference coefficient matrix.              C
C                  Dimension ( NFDCM, NR, NDRVM+1, NDCS ).           C
C                   Array is generated by the routine svfdcf         C
C                 See documentation for SVFDCF for details.          C
C                                                                    C
C     XARR      : Dim ( NR ). i^{th} element gives x_i.              C
C                                                                    C
C____________________________________________________________________C
C
C*********************************************************************
      SUBROUTINE MFDMFM( KL, NDCS, NBN, NBAND, N2, NR, NH, NDRVS,
     1                   NDRVM, NFDCM, IPIV, MHBCS, MHT, MHL, MHM,
     2                   MHP, DMAT, DC, DELTAT, CFAC, SVFDC, XARR )
      IMPLICIT NONE
C____________________________________________________________________C
C Variable declarations - Parameters ................................C
      INTEGER KL, NDCS, NBN, NBAND, N2, NR, NH, NDRVS, NDRVM, NFDCM,
     1        IPIV( N2 ), MHBCS( NDCS ),
     2        MHT( * ), MHL( * ), MHM( * ), MHP( * )
      DOUBLE PRECISION DMAT( NBAND, N2 ), DC, DELTAT, CFAC,
     1                 SVFDC( NFDCM, NR, NDRVM+1, NDCS ), XARR( NR )
C____________________________________________________________________C
C Variable declarations - Working variables .........................C
      INTEGER IMF, IOP, ILNR, IRNR, INFO, ITYPE, INARR( 3 )
      DOUBLE PRECISION ZERO, ONE, FAC, TOL
      PARAMETER ( ZERO = 0.0d0, ONE = 1.0d0, TOL = 1.0d-7 )
C____________________________________________________________________C
C START OF PROGRAM **************************************************C
C____________________________________________________________________C
C Check input parameters
C
      INARR( 1 ) = 4
      INARR( 2 ) = NR
      INARR( 3 ) = NH
C
      IF ( N2.NE.NR*NH ) THEN
        PRINT *,' Subroutine MFDMFM.'
        PRINT *,' N2 = ', N2
        PRINT *,' NR = ', NR
        PRINT *,' NH = ', NH
        PRINT *,' Program aborted.'
        STOP
      ENDIF
C
      IF ( NBAND.NE.(3*KL+1) ) THEN
        PRINT *,' Subroutine MFDMFM.'
        PRINT *,' NBAND = ', NBAND
        PRINT *,' KL = ', KL
        PRINT *,' Program aborted.'
        STOP
      ENDIF
C
C Check CFAC and DC
C
      IF ( CFAC.LT.ZERO .OR. CFAC.GT.ONE ) THEN
        PRINT *,' Subroutine MFDMFM.'
        PRINT *,' CFAC = ', CFAC
        PRINT *,' Program aborted.'
        STOP
      ENDIF
C
      IF ( DC.LT.TOL ) THEN
        PRINT *,' Subroutine MFDMFM.'
        PRINT *,' DC = ', DC
        PRINT *,' Program aborted.'
        STOP
      ENDIF
C
C Enforce matrix is built in LAPACK format
C
      IMF = 1
C
C Zero the matrix
C
      IOP = 0
      CALL MATOP( DMAT, ZERO, NBAND, N2, IOP )
C
C Now, prepare to build the diffusion matrix
C
      ILNR = 2
      IRNR = NR - 1
C
C First - poloidal part ....
C
      FAC = DC*DELTAT*(CFAC - 1.0d0)
      ITYPE = 4
      CALL AMLP( NR, INARR, MHT, MHL, MHM, MHP, ITYPE, MHT,
     1           MHL, MHM, FAC, NBN, NDRVS, NDRVM, NFDCM,
     2           ILNR, IRNR, SVFDC, XARR, NDCS, DMAT, NBAND,
     3           N2, IMF, KL, KL, KL )
C
      FAC = 1.0d0
      CALL AMTA( NR, INARR, MHT, MHL, MHM, MHP, ITYPE,
     1           MHT, MHL, MHM, ITYPE, FAC, ILNR, IRNR, DMAT,
     2           NBAND, N2, IMF, KL, KL, KL, NDRVS, NDRVM,
     3           NFDCM, NDCS, SVFDC, XARR, NBN )
C
C Now - toroidal part ....
C
      FAC = DC*DELTAT*(CFAC - 1.0d0)
      ITYPE = 5
      CALL AMLP( NR, INARR, MHT, MHL, MHM, MHP, ITYPE, MHT,
     1           MHL, MHM, FAC, NBN, NDRVS, NDRVM, NFDCM,
     2           ILNR, IRNR, SVFDC, XARR, NDCS, DMAT, NBAND,
     3           N2, IMF, KL, KL, KL )
C
      FAC = 1.0d0
      CALL AMTA( NR, INARR, MHT, MHL, MHM, MHP, ITYPE,
     1           MHT, MHL, MHM, ITYPE, FAC, ILNR, IRNR, DMAT,
     2           NBAND, N2, IMF, KL, KL, KL, NDRVS, NDRVM,
     3           NFDCM, NDCS, SVFDC, XARR, NBN )
C
C
C Now add diagonal elements to matrix
C
      FAC = 1.0d0
      CALL AMSDEA( DMAT, NBAND, N2, KL, KL, KL, IMF, INARR,
     1             MHP, MHBCS, 'Inner Boundary', FAC, NDCS )
      CALL AMSDEA( DMAT, NBAND, N2, KL, KL, KL, IMF, INARR,
     1             MHP, MHBCS, 'Outer Boundary', FAC, NDCS )
C
C
C Now we must attempt to do an LU decomposition
C of the DMAT matrix ... for this, we need the
C LAPACK routine DGBTRF
C
      CALL DGBTRF( N2, N2, KL, KL, DMAT, NBAND, IPIV, INFO )
C
C Check for an error from LU decomp.
C
      IF ( INFO.NE.0 ) THEN
        PRINT *,' Subroutine MFDMFM.'
        PRINT *,' The LAPACK subroutine DGBTRF has been '
        PRINT *,' and has returned ',INFO,' as a value of '
        PRINT *,' INFO. '
        PRINT *,' Program aborted.'
        STOP
      ENDIF
C
      RETURN
      END
C*********************************************************************
