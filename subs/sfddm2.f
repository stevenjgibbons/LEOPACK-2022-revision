C*********************************************************************
C subroutine Stress Free Dynamo Diffusion Matrix form 2 **************
C            -      -    -      -         -           - **************
C Steve Gibbons Wed Jul  5 18:47:21 BST 2000                         C
C____________________________________________________________________C
C                                                                    C
C If we are time stepping the induction, heat and vorticity          C
C equations, this routine gives us the matrix which will solve       C
C                                                                    C
C ( CA + dt.CD.(c-1) Lap ) Theta^{i+1} =                             C
C                            ( CA + dt.CD.c Lap ) Theta^{i} +        C
C ^^^^^^^^^^^^^^^^^^^^^^^^                         Forc. Term        C
C                                                                    C
C ( CE + dt.CI.(c-1) Lap ) omega^{i+1} =                             C
C                            ( CE + dt.CI.c Lap ) omega^{i} +        C
C ^^^^^^^^^^^^^^^^^^^^^^^^                         Forc. Term        C
C                                                                    C
C  (omega is the vorticity = curl u)                                 C
C                                                                    C
C ( CK + dt.CL.(c-1) Lap ) B^{i+1} =                                 C
C                            ( CK + dt.CL.c Lap ) B^{i} +            C
C ^^^^^^^^^^^^^^^^^^^^^^^^                         Forc. Term        C
C                                                                    C
C It zeros the matrix, assumes that all boundary conditions are      C
C insulating, velocity boundary conditions are rigid                 C
C and assumes a LAPACK format of the banded matrix (IMF = 1 in       C
C MATIND).                                                           C
C                                                                    C
C For a harmonic ih, MHT( ih ) must be 4 for poloidal and 5 for      C
C toroidal. MHL( ih ) gives the spherical harmonic degree, l.        C
C MHM( ih ) gives the spherical harmonic order, m, when ih has a     C
C (cos m phi) dependence and -m for a (sin m phi) dependence.        C
C                                                                    C
C MHP( ih ) = is where MHBCS( is ) = 6 when MHT( ih ) = 2 or 5,      C
C MHBCS( is ) = 7 when MHT( ih ) = 4 and MHBCS( is ) = 4 when        C
C MHT( ih ) = 1.                                                     C
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
C     N1     : Leading dimension of matrix, DMAT.                    C
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
C     NCFM     : Leading dimension of FDCM. At least (2*NBN+1)       C
C                                                                    C
C     IPIV      : Working array for LAPACK routine DGBTRF. Dim (N2). C
C                                                                    C
C     MHBCS     : Dim ( NDCS ). See above.                           C
C     MHT       : Dim ( NH ). See above.                             C
C     MHL       : Dim ( NH ). See above.                             C
C     MHM       : Dim ( NH ). See above.                             C
C     MHP       : Dim ( NH ). See above.                             C
C     MHTR      : Dim ( NH ). Type function for vorticity eqn.       C
C        a priori, you must CALL CINDSW ( NH, MHT, MHTR )            C
C                                                                    C
C     NTS       : Number of toroidal singularities.                  C
C                 Returned by the subroutine NSDSDR.                 C
C                                                                    C
C     IADC      : Dim (NTS). Contains the numbers of the rows        C
C                 which must be swapped for the angular momentum     C
C                 constraints. Not referred to if NTS = 0.           C
C                                                                    C
C  Double Precision                                                  C
C  ----------------                                                  C
C     DMAT      : Diffusion matrix. Dim ( N1, N2 ).                  C
C                                                                    C
C     CA        : Coef. of time deriv. in heat equation.             C
C     CE        : Coef. of time deriv. in vorticity equation.        C
C     CK        : Coef. of time deriv. in induction equation.        C
C                                                                    C
C     CD        : Coef. of diffusion in heat equation.               C
C     CI        : Coef. of diffusion in vorticity equation.          C
C     CL        : Coef. of diffusion in induction equation.          C
C                                                                    C
C     DELTAT    : Time step.                                         C
C                                                                    C
C     CFAC      : Determines how explicit/implicit time stepping is. C
C                 (Corresponds to 'c' in above equation).            C
C                                                                    C
C     SVFDC     : Finite difference coefficient matrix.              C
C                  Dimension ( NCFM, NR, NDRVM+1, NDCS ).            C
C                   Array is generated by the routine svfdcf         C
C                 See documentation for SVFDCF for details.          C
C                                                                    C
C     XARR      : Dim ( NR ). i^{th} element gives x_i.              C
C                                                                    C
C     V         : Dim ( N2, NTS ). Returned with vectors for         C
C                   Woodbury formula solution.                       C
C                                                                    C
C____________________________________________________________________C
C
C*********************************************************************
      SUBROUTINE SFDDM2( KL, NDCS, NBN, N1, N2, NR, NH, NDRVS,
     1                NDRVM, NCFM, IPIV, MHBCS, MHT, MHL, MHM, MHP,
     2                MHTR, DMAT, CA, CD, CE, CI, CK, CL, DELTAT,
     3                CFAC, SVFDC, XARR, NTS, IADC, V )
      IMPLICIT NONE
C____________________________________________________________________C
C Variable declarations - Parameters ................................C
      INTEGER KL, NDCS, NBN, N1, N2, NR, NH, NDRVS, NDRVM, NCFM,
     1        IPIV( N2 ), MHBCS( NDCS ), NTS, IADC( NTS ),
     2        MHT( * ), MHL( * ), MHM( * ), MHP( * ), MHTR( * )
      DOUBLE PRECISION DMAT( N1, N2 ), DELTAT, CFAC, XARR( NR ),
     1                 SVFDC( NCFM, NR, NDRVM+1, NDCS ),
     2                 CA, CD, CE, CI, CK, CL, V( N2, NTS )
C____________________________________________________________________C
C Variable declarations - Working variables .........................C
      INTEGER IMF, IOP, ILNR, IRNR, INFO, ICMP, INARR( 3 ),
     1        ILNT, IRNT, ICMPO, IRC, IND, ITS
      DOUBLE PRECISION ZERO, ONE, FAC
      PARAMETER ( ZERO = 0.0d0, ONE = 1.0d0 )
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
        PRINT *,' Subroutine SFDDM2.'
        PRINT *,' N2 = ', N2
        PRINT *,' NR = ', NR
        PRINT *,' NH = ', NH
        PRINT *,' Program aborted.'
        STOP
      ENDIF
C
      IF ( N1.NE.(3*KL+1) ) THEN
        PRINT *,' Subroutine SFDDM2.'
        PRINT *,' N1 = ', N1
        PRINT *,' KL = ', KL
        PRINT *,' Program aborted.'
        STOP
      ENDIF
C
C Check CFAC and DC
C
      IF ( CFAC.LT.ZERO .OR. CFAC.GT.ONE ) THEN
        PRINT *,' Subroutine SFDDM2.'
        PRINT *,' CFAC = ', CFAC
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
      CALL MATOP( DMAT, ZERO, N1, N2, IOP )
C
C Now, prepare to build the diffusion matrix
C
      ILNR = 2
      IRNR = NR - 1
C
      ILNT = 3
      IRNT = NR - 2
C
C Temperature
C
      FAC = CD*DELTAT*(CFAC - 1.0d0)
      ICMP = 3
      CALL AMLP( NR, INARR, MHT, MHL, MHM, MHP, ICMP, MHTR,
     1           MHL, MHM, FAC, NBN, NDRVS, NDRVM, NCFM,
     2           ILNR, IRNR, SVFDC, XARR, NDCS, DMAT, N1,
     3           N2, IMF, KL, KL, KL )
C
      FAC = CA
      CALL AMTA( NR, INARR, MHT, MHL, MHM, MHP, ICMP,
     1           MHTR, MHL, MHM, ICMP, FAC, ILNR, IRNR, DMAT,
     2           N1, N2, IMF, KL, KL, KL, NDRVS, NDRVM,
     3           NCFM, NDCS, SVFDC, XARR, NBN )
C
C Poloidal magnetic field part ....
C
      FAC = CL*DELTAT*(CFAC - 1.0d0)
      ICMP = 4
      CALL AMLP( NR, INARR, MHT, MHL, MHM, MHP, ICMP, MHTR,
     1           MHL, MHM, FAC, NBN, NDRVS, NDRVM, NCFM,
     2           ILNR, IRNR, SVFDC, XARR, NDCS, DMAT, N1,
     3           N2, IMF, KL, KL, KL )
C
      FAC = CK
      CALL AMTA( NR, INARR, MHT, MHL, MHM, MHP, ICMP,
     1           MHTR, MHL, MHM, ICMP, FAC, ILNR, IRNR, DMAT,
     2           N1, N2, IMF, KL, KL, KL, NDRVS, NDRVM,
     3           NCFM, NDCS, SVFDC, XARR, NBN )
C
C Toroidal magnetic field part ....
C
      FAC = CL*DELTAT*(CFAC - 1.0d0)
      ICMP = 5
      CALL AMLP( NR, INARR, MHT, MHL, MHM, MHP, ICMP, MHTR,
     1           MHL, MHM, FAC, NBN, NDRVS, NDRVM, NCFM,
     2           ILNR, IRNR, SVFDC, XARR, NDCS, DMAT, N1,
     3           N2, IMF, KL, KL, KL )
C
      FAC = CK
      CALL AMTA( NR, INARR, MHT, MHL, MHM, MHP, ICMP,
     1           MHTR, MHL, MHM, ICMP, FAC, ILNR, IRNR, DMAT,
     2           N1, N2, IMF, KL, KL, KL, NDRVS, NDRVM,
     3           NCFM, NDCS, SVFDC, XARR, NBN )
C
C Poloidal and toroidal velocity parts ....
C
      FAC = CI*DELTAT*(CFAC - 1.0d0)
      CALL AMLC( NR, INARR, MHT, MHL, MHM, MHP, MHTR, MHL,
     1           MHM, FAC, NBN, NDRVS, NDRVM, ILNR, IRNR,
     2           ILNT, IRNT, XARR, NCFM, DMAT, N1, N2,
     3           IMF, KL, KL, KL, NDCS, SVFDC )
C
      FAC = CE
      ICMP  = 1
      ICMPO = 2
      CALL AMCL( NR, INARR, MHT, MHL, MHM, MHP, ICMP, MHTR,
     1           MHL, MHM, ICMPO, FAC, NBN, NDRVS, NDRVM,
     2           ILNT, IRNT, XARR, NCFM, SVFDC, DMAT, N1, N2,
     3           IMF, KL, KL, KL, NDCS )
C
      ICMP  = 2
      ICMPO = 1
      CALL AMCL( NR, INARR, MHT, MHL, MHM, MHP, ICMP, MHTR,
     1           MHL, MHM, ICMPO, FAC, NBN, NDRVS, NDRVM,
     2           ILNR, IRNR, XARR, NCFM, SVFDC, DMAT, N1, N2,
     3           IMF, KL, KL, KL, NDCS )
C
C Now add diagonal elements to matrix
C
      FAC = 1.0d0
      CALL AMSDEA( DMAT, N1, N2, KL, KL, KL, IMF, INARR,
     1             MHP, MHBCS, 'Inner Boundary', FAC, NDCS )
      CALL AMSDEA( DMAT, N1, N2, KL, KL, KL, IMF, INARR,
     1             MHP, MHBCS, 'Outer Boundary', FAC, NDCS )
C
C If NTS = 0, we do not to perform any of the following
C operations and can jump directly to line 50.
C
      IF ( NTS.EQ.0 ) GOTO 50
C
C OK so NTS is non-zero.
C Must find the elements in the array V.
C
      CALL BMSCD( N1, N2, NTS, KL, KL, KL, IPIV, IADC, DMAT, V )
C
C Now must replace the matrix DMAT again
C Zero the matrix
C
      IOP = 0
      CALL MATOP( DMAT, ZERO, N1, N2, IOP )
C
C Now, prepare to build the diffusion matrix
C
      ILNR = 2
      IRNR = NR - 1
C
      ILNT = 3
      IRNT = NR - 2
C
C Temperature
C
      FAC = CD*DELTAT*(CFAC - 1.0d0)
      ICMP = 3
      CALL AMLP( NR, INARR, MHT, MHL, MHM, MHP, ICMP, MHTR,
     1           MHL, MHM, FAC, NBN, NDRVS, NDRVM, NCFM,
     2           ILNR, IRNR, SVFDC, XARR, NDCS, DMAT, N1,
     3           N2, IMF, KL, KL, KL )
C
      FAC = CA
      CALL AMTA( NR, INARR, MHT, MHL, MHM, MHP, ICMP,
     1           MHTR, MHL, MHM, ICMP, FAC, ILNR, IRNR, DMAT,
     2           N1, N2, IMF, KL, KL, KL, NDRVS, NDRVM,
     3           NCFM, NDCS, SVFDC, XARR, NBN )
C
C Poloidal magnetic field part ....
C
      FAC = CL*DELTAT*(CFAC - 1.0d0)
      ICMP = 4
      CALL AMLP( NR, INARR, MHT, MHL, MHM, MHP, ICMP, MHTR,
     1           MHL, MHM, FAC, NBN, NDRVS, NDRVM, NCFM,
     2           ILNR, IRNR, SVFDC, XARR, NDCS, DMAT, N1,
     3           N2, IMF, KL, KL, KL )
C
      FAC = CK
      CALL AMTA( NR, INARR, MHT, MHL, MHM, MHP, ICMP,
     1           MHTR, MHL, MHM, ICMP, FAC, ILNR, IRNR, DMAT,
     2           N1, N2, IMF, KL, KL, KL, NDRVS, NDRVM,
     3           NCFM, NDCS, SVFDC, XARR, NBN )
C
C Toroidal magnetic field part ....
C
      FAC = CL*DELTAT*(CFAC - 1.0d0)
      ICMP = 5
      CALL AMLP( NR, INARR, MHT, MHL, MHM, MHP, ICMP, MHTR,
     1           MHL, MHM, FAC, NBN, NDRVS, NDRVM, NCFM,
     2           ILNR, IRNR, SVFDC, XARR, NDCS, DMAT, N1,
     3           N2, IMF, KL, KL, KL )
C
      FAC = CK
      CALL AMTA( NR, INARR, MHT, MHL, MHM, MHP, ICMP,
     1           MHTR, MHL, MHM, ICMP, FAC, ILNR, IRNR, DMAT,
     2           N1, N2, IMF, KL, KL, KL, NDRVS, NDRVM,
     3           NCFM, NDCS, SVFDC, XARR, NBN )
C
C Poloidal and toroidal velocity parts ....
C
      FAC = CI*DELTAT*(CFAC - 1.0d0)
      CALL AMLC( NR, INARR, MHT, MHL, MHM, MHP, MHTR, MHL,
     1           MHM, FAC, NBN, NDRVS, NDRVM, ILNR, IRNR,
     2           ILNT, IRNT, XARR, NCFM, DMAT, N1, N2,
     3           IMF, KL, KL, KL, NDCS, SVFDC )
C
      FAC = CE
      ICMP  = 1
      ICMPO = 2
      CALL AMCL( NR, INARR, MHT, MHL, MHM, MHP, ICMP, MHTR,
     1           MHL, MHM, ICMPO, FAC, NBN, NDRVS, NDRVM,
     2           ILNT, IRNT, XARR, NCFM, SVFDC, DMAT, N1, N2,
     3           IMF, KL, KL, KL, NDCS )
C
      ICMP  = 2
      ICMPO = 1
      CALL AMCL( NR, INARR, MHT, MHL, MHM, MHP, ICMP, MHTR,
     1           MHL, MHM, ICMPO, FAC, NBN, NDRVS, NDRVM,
     2           ILNR, IRNR, XARR, NCFM, SVFDC, DMAT, N1, N2,
     3           IMF, KL, KL, KL, NDCS )
C
C Now add diagonal elements to matrix
C
      FAC = 1.0d0
      CALL AMSDEA( DMAT, N1, N2, KL, KL, KL, IMF, INARR,
     1             MHP, MHBCS, 'Inner Boundary', FAC, NDCS )
      CALL AMSDEA( DMAT, N1, N2, KL, KL, KL, IMF, INARR,
     1             MHP, MHBCS, 'Outer Boundary', FAC, NDCS )
C
C Now need to remove the rows corresponding to the
C singularities we have dealt with above.
C
      IRC = 1
      DO ITS = 1, NTS
        IND = IADC( ITS )
        CALL BMRCOP( KL, KL, KL, N2, IRC, IND, DMAT, ZERO, ZERO )
        DMAT( 2*KL + 1, IND ) = ONE
        V( IND, ITS ) = V( IND, ITS ) - ONE
      ENDDO
C
C Final LU decomposition for actual diffusion matrix.
C
 50   CONTINUE
C
C Now we must attempt to do an LU decomposition
C of the DMAT matrix ... for this, we need the
C LAPACK routine DGBTRF
C
      CALL DGBTRF( N2, N2, KL, KL, DMAT, N1, IPIV, INFO )
C
C Check for an error from LU decomp.
C
      IF ( INFO.NE.0 ) THEN
        PRINT *,' Subroutine SFDDM2.'
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
