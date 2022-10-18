C*********************************************************************
C subroutine Non-uniform Grid Matrix CurL ****************************
C            -           -    -      -  - ****************************
C Steve Gibbons Fri Oct 15 13:33:10 BST 1999                         C
C____________________________________________________________________C
C                                                                    C
C Adds FAC to all the elements of the matrix A which, when           C
C multiplied by the vector VI adds the curl of all harmonics IHI     C
C in solution vector VI (satisfying MHTI( IHI ) = ICMPI) to the      C
C corresponding harmonics IHO in vector VO (which satsify            C
C MHTO( IHO ) = ICMPO.                                               C
C                                                                    C
C Multiplying VI by A after running NGMCL should be equivalent to    C
C calling NGSVCL with VI.                                            C
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
C and similarly MHTO for the resulting vector.                       C
C                                                                    C
C The program aborts if ICMPI = 3 as curl has no meaning.            C
C                                                                    C
C Taking the curl of a toroidal harmonic with radial function        C
C x( r ) results in a poloidal harmonic with r.f. x(r)               C
C                                                                    C
C Hence ICMPI = 2(5) would normally correspond to ICMPO = 1(4)       C
C                                                                    C
C                                                                    C
C Taking the curl of a poloidal harmonic with radial function        C
C x( r ) results in a toroidal harmonic with r.f. -D_l x(r)          C
C                                                                    C
C Hence ICMPI = 1(4) would normally correspond to ICMPO = 2(5)       C
C                                                                    C
C Note that as this is a matrix operation (unlike in the RHS op.s),  C
C the harmonic sets must be equivalent and so only one array,        C
C INARR, must correspond to both VI and VO.                          C
C                                                                    C
C____________________________________________________________________C
C                                                                    C
C Input variables :-                                                 C
C ===============                                                    C
C  Integer                                                           C
C  -------                                                           C
C                                                                    C
C     NR        : Number of radial grid nodes.                       C
C     INARR     : Int. parameter array corresponding to vectors.     C
C                 Dim ( * ) but the following must be true.          C
C                                                                    C
C                 INARRI( 1 ) = IFORMF                               C
C                 INARRI( 2 ) = NR      See INDFUN for details       C
C                 INARRI( 3 ) = NH                                   C
C                                                                    C
C     MHTI      : Array length ( * ) - atleast length NHI            C
C                  See above for key. (corresponds to input vec.)    C
C     MHLI      : Array length ( * ) - atleast length NHI            C
C                  Sph. harm. degree, l.                             C
C     MHMI      : Array length ( * ) - atleast length NHI            C
C                  Sph. harm. order, m, for cos m phi dep.           C
C                 -Sph. harm. order, m, for sin m phi dep.           C
C                                                                    C
C     ICMPI     : The type of harmonic we wish to take curl of.      C
C                                                 see above.         C
C                                                                    C
C     Note that MHTO, MHLO, MHMO and ICMPO correspond to             C
C    exactly the above variables but for the output and not input    C
C    vectors. ICMPI cannot be 3. NRI must equal NRO and must both    C
C    equal NR.                                                       C
C                                                                    C
C     NBN       : Number of nodes on each side of point for          C
C                  central differences.                              C
C                                                                    C
C     NDRVS     : Highest derivative stored in FDCM.                 C
C                (Not needed if we are doing tor --> pol but         C
C                 must be atleast 2 if doing pol --> tor )           C
C                                                                    C
C     NFDCM     : Leading dim of FDCM. See FDCMBD.                   C
C                  (Must be atleast 2*NBN + 1 )                      C
C                                                                    C
C     ILNR      : First radial node to act on.                       C
C     IRNR      : Last radial node to act on.                        C
C                                                                    C
C     ILNC      : First radial node which may be used for f.d. form. C
C     IRNC      : Last radial node which may be used for f.d. form.  C
C                                                                    C
C  ILNC and IRNC correspond to NLMC and NRMC in the call to FDCMBD.  C
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
C     FDCM      : Dim ( NFDCM, NR, NDRVS ). Finite diff. coeff.s     C
C                (Not needed if we are doing tor --> pol).           C
C                Must be formed in advance by a call to FDCMBD.      C
C                                                                    C
C____________________________________________________________________C
C
C*********************************************************************
      SUBROUTINE NGMCL( NR, INARR, MHTI, MHLI, MHMI, ICMPI, MHTO,
     1                  MHLO, MHMO, ICMPO, FAC, NBN, NDRVS, ILNR,
     2                  IRNR, ILNC, IRNC, FDCM, XARR, NFDCM,
     3                  A, N1, N2, IMF, KL, KU, KLE )
      IMPLICIT NONE
C____________________________________________________________________C
C Variable declarations - Parameters ................................C
      INTEGER NR, INARR( * ), MHTI( * ), MHLI( * ), MHMI( * ),
     1        ICMPI, MHTO( * ), MHLO( * ), MHMO( * ),
     2        ICMPO, NBN, NDRVS, ILNR, IRNR, ILNC, IRNC, NFDCM,
     3        N1, N2, IMF, KL, KU, KLE
      DOUBLE PRECISION FAC, XARR( NR ), FDCM( NFDCM, NR, NDRVS ),
     1                 A( N1, N2 )
C____________________________________________________________________C
C Variable declarations - Working variables .........................C
C
      INTEGER IFORMF, NRR, NH, IHI, IHO,
     1        IHD, INARR2( 6 )
      DOUBLE PRECISION LOW, F2, JUNK( 1 ), WORK( 3 )
      PARAMETER ( LOW = 1.0d-9 )
      EXTERNAL MATDLT
C____________________________________________________________________C
C START OF PROGRAM **************************************************C
C____________________________________________________________________C
C
C Don't bother checking ILNR, ILNC, IRNR, IRNC, NDRVS or NBN
C as these errors will all be trapped by NGSVDR
C
C     .
C     . However do check ICMPI
C     .
C
      IF ( ICMPI.NE.1 .AND. ICMPI.NE.2 .AND. ICMPI.NE.4
     1                .AND. ICMPI.NE.5           )    THEN
         PRINT *,' Subroutine NGMCL.'
         PRINT *,' ICMPI = ', ICMPI
         PRINT *,' Program aborted.'
         STOP
      ENDIF
C
      IFORMF = INARR( 1 )
      NRR    = INARR( 2 )
      NH     = INARR( 3 )
C
      IF ( IFORMF.EQ.1 .OR. IFORMF.EQ.2 ) THEN
         PRINT *,' Subroutine NGMCL.'
         PRINT *,' IFORMF = ', IFORMF
         PRINT *,' Program aborted.'
         STOP
      ENDIF
C
      IF ( NRR.NE.NR ) THEN
         PRINT *,' Subroutine NGMCL.'
         PRINT *,' NR   = ', NR
         PRINT *,' NRR  = ', NRR
         PRINT *,' Program aborted.'
         STOP
      ENDIF
C
C     . early exit?
C
      IF ( ABS(FAC).LT.LOW ) RETURN
C     .
C     . Loop around in harmonics
C     .
      DO IHI = 1, NH
        IF ( MHTI( IHI ).NE.ICMPI ) GOTO 50
C
C ok - so our harmonic is of the correct type
C so let's look for its corresponding harmonic
C in VO
C
        DO IHO = 1, NH
          IF (   MHTO( IHO ).EQ.ICMPO          .AND.
     1           MHLO( IHO ).EQ.MHLI( IHI )    .AND.
     2           MHMO( IHO ).EQ.MHMI( IHI )         ) THEN
C
C o.k. we've found the corresponding harmonic
C Set INARR and F2 variables ready to call NGMICA
C
            INARR2( 1 ) = IFORMF
            INARR2( 2 ) = NR
            INARR2( 3 ) = NH
C           .
            IF ( MHTI( IHI ).EQ.1 .OR. MHTI( IHI ).EQ.4 ) THEN
C             .
C             . We are curling a poloidal harmonic
C             .
              IHD         = 2
              INARR2( 5 ) = 2
              INARR2( 6 ) = MHLI( IHI )
              F2          = (-1.0d0)*FAC
            ELSE
C             .
C             . We are curling a toroidal harmonic
C             .
              IHD         = 0
              INARR2( 5 ) = 1
              F2          = FAC
C             .
            ENDIF
C           .
            CALL NGMICA( N1, N2, KL, KU, KLE, IMF, IHI, IHO,
     1                   INARR2, IHD, NBN, ILNR, IRNR, ILNC,
     2                   IRNC, NFDCM, NDRVS, MATDLT, A, F2,
     3                   XARR, FDCM, JUNK, JUNK, NR, WORK )
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
