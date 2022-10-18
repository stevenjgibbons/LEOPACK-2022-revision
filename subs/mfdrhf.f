C*********************************************************************
C subroutine Magnetic Field Diffusion Right Hand Form ****************
C            -        -     -         -     -    -    ****************
C Steve Gibbons Tue Nov  2 10:53:43 GMT 1999                         C
C____________________________________________________________________C
C                                                                    C
C If we are time stepping the induction equation, this routine       C
C adds the highlighted term to the right hand side of                C
C                                                                    C
C ( 1 + dt (c-1) Lap ) B^{i+1} = ( 1 + dt c Lap ) B^{i} + Forc. Term C
C                                ^^^^^^^^^^^^^^^^                    C
C It does NOT zero the array RHS.                                    C
C                                                                    C
C It assumes an insulating boundary at either end and therefore      C
C a degeneracy of 1 grid node at each boundary - we therefore set    C
C ILNR = 2 and IRNR = NR - 1.                                        C
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
C     NR        : Number of radial grid nodes.                       C
C     NH        : Number of spherical harmonic radial functions.     C
C                                                                    C
C     NDCS      : Maximum distinct finite difference schemes         C
C                  stored by the array SVFDC.                        C
C                                                                    C
C     MHT       : Dim ( NH ). See above.                             C
C     MHL       : Dim ( NH ). See above.                             C
C     MHM       : Dim ( NH ). See above.                             C
C     MHP       : Dim ( NH ). See above.                             C
C                                                                    C
C     NBN       : Number of bounding nodes.                          C
C                                                                    C
C     NDRVS     : Number of highest derivative for which             C
C                  coefficients are stored by the array SVFDC.       C
C                                                                    C
C     NDRVM     : Limit on NDRVS. Array bound for SVFDC.             C
C                                                                    C
C     NFDCM     : Leading dimension of FDCM. At least (2*NBN+1)      C
C                                                                    C
C  Double Precision                                                  C
C  ----------------                                                  C
C                                                                    C
C     DC        : Diffusion coefficient. Must be positive.           C
C                                                                    C
C     CFAC      : Determines how explicit/implicit time stepping is. C
C                 (Corresponds to 'c' in above equation).            C
C                                                                    C
C     DELTAT    : Time step.                                         C
C                                                                    C
C     BVEC      : Magnetic field vector at step i                    C
C     RHS       : Right hand side vector.                            C
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
      SUBROUTINE MFDRHF( NR, NH, NDCS, MHT, MHL, MHM, MHP, NBN,
     1                   NDRVS, NDRVM, NFDCM, DC, CFAC, DELTAT, BVEC,
     2                   RHS, SVFDC, XARR )
      IMPLICIT NONE
C____________________________________________________________________C
C Variable declarations - Parameters ................................C
      INTEGER NR, NH, NDCS, MHT( * ), MHL( * ), MHM( * ), MHP( * ),
     1        NBN, NDRVS, NDRVM, NFDCM
      DOUBLE PRECISION DC, CFAC, DELTAT, BVEC( * ), RHS( * ),
     1                 SVFDC( NFDCM, NR, NDRVM+1, NDCS ), XARR( NR )
C____________________________________________________________________C
C Variable declarations - Working variables .........................C
      INTEGER INARR( 3 ), ITYPE, ILNR, IRNR
      DOUBLE PRECISION FAC
C____________________________________________________________________C
C START OF PROGRAM **************************************************C
C____________________________________________________________________C
C
      INARR( 1 ) = 4
      INARR( 2 ) = NR
      INARR( 3 ) = NH
C     .
C     . Now we are ready to form the right hand
C     . side - first add Laplacian of poloidal terms
C     .
C
      ILNR = 2
      IRNR = NR - 1
C
      ITYPE = 4
      FAC = DC*CFAC*DELTAT
C
      CALL ASVLP( NR, NDCS, BVEC, INARR, MHT, MHL, MHM, ITYPE,
     1            MHP, RHS, INARR, MHT, MHL, MHM, FAC, NBN,
     2            NDRVS, NDRVM, ILNR, IRNR, SVFDC, XARR, NFDCM )
C
C  add delta_t
C
      FAC = 1.0d0
      CALL ASVTA( NR, NDCS, BVEC, INARR, MHT, MHL, MHM, ITYPE,
     1            MHP, RHS, INARR, MHT, MHL, MHM, ITYPE, FAC,
     2            ILNR, IRNR, NBN, NDRVS, NDRVM, NFDCM, SVFDC )
C
C Now the toroidal component
C
      ITYPE = 5
      FAC = DC*CFAC*DELTAT
C
      CALL ASVLP( NR, NDCS, BVEC, INARR, MHT, MHL, MHM, ITYPE,
     1            MHP, RHS, INARR, MHT, MHL, MHM, FAC, NBN,
     2            NDRVS, NDRVM, ILNR, IRNR, SVFDC, XARR, NFDCM )
C
C  add delta_t
C
      FAC = 1.0d0
      CALL ASVTA( NR, NDCS, BVEC, INARR, MHT, MHL, MHM, ITYPE,
     1            MHP, RHS, INARR, MHT, MHL, MHM, ITYPE, FAC,
     2            ILNR, IRNR, NBN, NDRVS, NDRVM, NFDCM, SVFDC )
C
      RETURN
      END
C*********************************************************************
