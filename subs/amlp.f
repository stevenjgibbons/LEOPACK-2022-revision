C*********************************************************************
C subroutine Adapted Matrix LaPlacian ********************************
C            -       -      - -       ********************************
C Steve Gibbons Wed Oct 27 07:04:48 GMT 1999                         C
C____________________________________________________________________C
C                                                                    C
C Adds FAC to all the elements of the matrix A which will add the    C
C laplacian of all harmonics IHI in solution vector VI (which        C
C satisfy MHTI( IHI ) = ICMP) to the corresponding harmonics IHO in  C
C vector VO (with MHTO( IHO ) = ICMP).                               C
C                                                                    C
C Multiplying VI by A (to form VO) should have the same result as    C
C calling NGSVLP with VI and VO.                                     C
C                                                                    C
C                                                                    C
C Now MHTI defines what each scalar function in a solution vector    C
C represents.                                                        C
C                                                                    C
C         MHTI( IH ) = 1 --> harmonic is poloidal velocity           C
C         MHTI( IH ) = 2 --> harmonic is toroidal velocity           C
C         MHTI( IH ) = 3 --> harmonic is temperature.                C
C         MHTI( IH ) = 4 --> harmonic is poloidal magnetic field     C
C         MHTI( IH ) = 5 --> harmonic is toroidal magnetic field     C
C                                                                    C
C The lapacian of a scalar with radial function f(r)                 C
C is a scalar with radial function D_l f(r) - l is degree of S.H.    C
C                                                                    C
C The lapacian of a toroidal harmonic with radial function f(r)      C
C is toroidal with radial function D_l f(r) - l is degree of S.H.    C
C                                                                    C
C The lapacian of a poloidal harmonic with radial function f(r)      C
C is poloidal with radial function D_l f(r) - l is degree of S.H.    C
C                                                                    C
C  This routine is therefore valid for all forms of harmonics        C
C (itype = 1, 2, 3, 4 and 5) and the result must go in the same      C
C type of harmonic as the original.                                  C
C                                                                    C
C____________________________________________________________________C
C                                                                    C
C Input variables :-                                                 C
C ===============                                                    C
C  Integer                                                           C
C  -------                                                           C
C                                                                    C
C     NR        : Number of radial grid nodes.                       C
C     NDCS      : Maximum distinct finite difference schemes         C
C                  stored by the array SVFDC.                        C
C     INARR     : Int. parameter array corresponding to VI and VO.   C
C                 Dim ( * ) but the following must be true.          C
C                                                                    C
C                 INARR( 1 ) = IFORMF                                C
C                 INARR( 2 ) = NRR     See INDFUN for details        C
C                 INARR( 3 ) = NH                                    C
C                                                                    C
C     MHTI      : Array length ( * ) - atleast length NH             C
C                  See above for key. (corresponds to input vec.)    C
C     MHLI      : Array length ( * ) - atleast length NH             C
C                  Sph. harm. degree, l.                             C
C     MHMI      : Array length ( * ) - atleast length NH             C
C                  Sph. harm. order, m, for cos m phi dep.           C
C                 -Sph. harm. order, m, for sin m phi dep.           C
C                                                                    C
C     MHPI      : Array length ( * ) - atleast length NHI            C
C                  Pointer array to finite difference coefficients.  C
C                  MHPI( ih ) = is, which is the 4th index of        C
C                  array SVFDC - indicates f.d. scheme used.         C
C                                                                    C
C     ICMP     : The type of harmonic we wish to take curl of.       C
C                                                 see above.         C
C                                                                    C
C     Note that MHTO, MHLO and MHMO correspond to                    C
C    exactly the above variables but for the output and not input    C
C    vectors. NRI must equal NRO and must both equal NR.             C
C                                                                    C
C     NBN       : Number of nodes on each side of point for          C
C                  central differences.                              C
C                                                                    C
C     NDRVS     : Highest derivative stored in SVFDC.                C
C                (Not needed if we are doing tor --> pol but         C
C                 must be atleast 2 if doing pol --> tor )           C
C     NDRVM     : Limit on NDRVS. Array bound for SVFDC.             C
C                                                                    C
C     NFDCM     : Leading dim of SVFDC. See SVFDCF.                  C
C                  (Must be atleast 2*NBN + 1 )                      C
C                                                                    C
C     ILNR      : First radial node to act on.                       C
C     IRNR      : Last radial node to act on.                        C
C                                                                    C
C     N1        : First dimension of matrix, A.                      C
C     N2        : Second dimension of matrix, A.                     C
C     IMF       : Matrix format flag. (See MATIND)                   C
C     KL        : Number of lower diagonals in matrix                C
C     KU        : Number of upper diagonals in matrix                C
C     KLE       : Number of additional lower diagonals in matrix     C
C                                                                    C
C  Double Precision                                                  C
C  ----------------                                                  C
C                                                                    C
C     A         : Matrix. Dim ( N1, N2 )                             C
C     FAC       : Multiplier of curl to be added.                    C
C     XARR      : Dim ( NR ). i^{th} element gives x_i.              C
C     SVFDC     : Finite difference coefficient matrix.              C
C                  Dimension ( NFDCM, NR, NDRVM+1, NDCS ).           C
C                   Array is generated by the routine svfdcf         C
C                 See documentation for SVFDCF for details.          C
C                                                                    C
C____________________________________________________________________C
C
C*********************************************************************
      SUBROUTINE AMLP( NR, INARR, MHTI, MHLI, MHMI, MHPI, ICMP, MHTO,
     1                 MHLO, MHMO, FAC, NBN, NDRVS, NDRVM, NFDCM,
     2                 ILNR, IRNR, SVFDC, XARR, NDCS, A, N1, N2, IMF,
     3                 KL, KU, KLE )
      IMPLICIT NONE
C____________________________________________________________________C
C Variable declarations - Parameters ................................C
      INTEGER NR, INARR( * ), MHTI( * ), MHLI( * ), MHMI( * ),
     1        MHPI( * ), ICMP, MHTO( * ), MHLO( * ), MHMO( * ),
     2        NBN, NDRVS, NDRVM, NFDCM, ILNR, IRNR,
     3        N1, N2, IMF, KL, KU, KLE, NDCS
      DOUBLE PRECISION FAC, XARR( NR ), A( N1, N2 ),
     1                 SVFDC( NFDCM, NR, NDRVM+1, NDCS )
C____________________________________________________________________C
C Variable declarations - Working variables .........................C
C
      INTEGER IFORMF, NRR, NH, IHI, IHO, IHD, IPARS( 2 ), IS
      DOUBLE PRECISION LOW, JUNK( 1 ), WORK( 3 )
      EXTERNAL AMDLT
      PARAMETER ( LOW = 1.0d-9 )
C____________________________________________________________________C
C START OF PROGRAM **************************************************C
C____________________________________________________________________C
C
      JUNK( 1 ) = 0.0d0
C
C Don't bother checking ILNR, IRNR, NDRVS or NBN
C as these errors will all be trapped by NGSVDR
C
C     .
C     . However do check ICMP
C     .
C
      IF ( ICMP.NE.1 .AND. ICMP.NE.2 .AND. ICMP.NE.3 .AND.
     1        ICMP.NE.4 .AND. ICMP.NE.5           )    THEN
         PRINT *,' Subroutine AMLP.'
         PRINT *,' ICMP = ', ICMP
         PRINT *,' Program aborted.'
         STOP
      ENDIF
C
      IFORMF = INARR( 1 )
      NRR    = INARR( 2 )
      NH     = INARR( 3 )
C
      IF ( IFORMF.EQ.1 .OR. IFORMF.EQ.2 ) THEN
         PRINT *,' Subroutine AMLP.'
         PRINT *,' IFORMF = ', IFORMF
         PRINT *,' Program aborted.'
         STOP
      ENDIF
C
      IF ( NRR.NE.NR ) THEN
         PRINT *,' Subroutine AMLP.'
         PRINT *,' NR   = ', NR
         PRINT *,' NRR  = ', NRR
         PRINT *,' Program aborted.'
         STOP
      ENDIF
C
C     . early exit?
C
      IF ( ABS(FAC).LT.LOW ) RETURN
      IHD = 2
C     .
C     . Loop around in harmonics
C     .
      DO IHI = 1, NH
        IF ( MHTI( IHI ).NE.ICMP ) GOTO 50
C
C ok - so our harmonic is of the correct type
C so let's look for its corresponding harmonic
C in VO
C
        IS = MHPI( IHI )
C
        DO IHO = 1, NH
          IF (   MHTO( IHO ).EQ.ICMP          .AND.
     1           MHLO( IHO ).EQ.MHLI( IHI )   .AND.
     2           MHMO( IHO ).EQ.MHMI( IHI )         ) THEN
C
C o.k. we've found the corresponding harmonic
C set INARR2 variables ready to call NGMICA
C           .
            IPARS( 1 ) = 2
            IPARS( 2 ) = MHLI( IHI )
C           .
            CALL AMLICA( N1, N2, KL, KU, KLE, IMF, IHI, IHO, INARR,
     1                   IHD, NBN, ILNR, IRNR, NFDCM, NR, NDCS, IS,
     2                   NDRVS, NDRVM, IPARS, AMDLT, A, FAC, XARR,
     3                   WORK, JUNK, SVFDC )
C           .
          ENDIF 
        ENDDO
C
 50   CONTINUE
      ENDDO
C     .
      RETURN
      END
C*********************************************************************
