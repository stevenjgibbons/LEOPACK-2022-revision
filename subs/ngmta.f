C*********************************************************************
C subroutine Non-uniform Grid Matrix Term Addition *******************
C            -           -    -      -    -        *******************
C Steve Gibbons Fri Oct 15 08:14:18 BST 1999                         C
C____________________________________________________________________C
C                                                                    C
C Simply adds a factor of FAC to the elements of a matrix A which    C
C multiplies VI( i ) to give VO( i ). i.e. multiplying vector VI     C
C by A (after a call to NGMTA) to form VO  should be equivalent      C
C to calling NGSVTA with VI and VO. Acts on each                     C
C harmonic, IH, such that MHTI( ih ) = ICMPI and grid node such that C
C ILNR .le. ir .le. IRNR - and j is every harmonic, JH, such that    C
C MHTO( jh ) = ICMPO and grid node, jr, such that                    C
C  ILNR .le. ir .le. IRNR                                            C
C                                                                    C
C Useful maybe if adding time derivatives to the forcing term or     C
C buoyancy to toroidal component of the vorticity.                   C
C____________________________________________________________________C
C                                                                    C
C Input variables :-                                                 C
C ===============                                                    C
C  Integer                                                           C
C  -------                                                           C
C                                                                    C
C     NR        : Number of radial grid nodes.                       C
C     INARR     : Int. parameter array corresponding to VI and VO.   C
C                 Dim ( * ) but the following must be true.          C
C                                                                    C
C                 INARR( 1 ) = IFORMF                                C
C                 INARR( 2 ) = NRR     See INDFUN for details        C
C                 INARR( 3 ) = NH                                    C
C                                                                    C
C     MHTI      : Array length ( * ) - atleast length NH             C
C                  See below for key. (corresponds to input vec.)    C
C                                                                    C
C     MHTI( i ) = 1 for a poloidal velocity harmonic, i.             C
C     MHTI( i ) = 2 for a toroidal velocity harmonic, i.             C
C     MHTI( i ) = 3 for a temperature harmonic, i.                   C
C     MHTI( i ) = 4 for a poloidal magnetic field harmonic, i.       C
C     MHTI( i ) = 5 for a toroidal magnetic field harmonic, i.       C
C                                                                    C
C     MHLI      : Array length ( * ) - atleast length NH             C
C                  Sph. harm. degree, l.                             C
C     MHMI      : Array length ( * ) - atleast length NH             C
C                  Sph. harm. order, m, for cos m phi dep.           C
C                 -Sph. harm. order, m, for sin m phi dep.           C
C                                                                    C
C     Note that INARRO, MHTO, MHLO and MHMO correspond to            C
C    exactly the above variables but for the output and not input    C
C    vectors. NRI must equal NRO and must both equal NR.             C
C                                                                    C
C     ILNR      : First radial node to act on.                       C
C     IRNR      : Last radial node to act on.                        C
C                                                                    C
C     ICMPI     : The type of harmonic we wish to act upon.          C
C     ICMPO     : The type of harmonic we wish add term to.          C
C                                                                    C
C Note that ICMPI and ICMPO should generally be identical.           C
C An example of an exception would be when adding the buoyancy       C
C term to the vorticity - here, ICMPI would be 3 and ICMPO 2.        C
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
C     FAC       : Multiplication factor.                             C
C                                                                    C
C____________________________________________________________________C
C
C*********************************************************************
      SUBROUTINE NGMTA( NR, INARR, MHTI, MHLI, MHMI, ICMPI, MHTO,
     1                  MHLO, MHMO, ICMPO, FAC, ILNR, IRNR, A,
     2                  N1, N2, IMF, KL, KU, KLE )
      IMPLICIT NONE
C____________________________________________________________________C
C Variable declarations - Parameters ................................C
      INTEGER NR, INARR( * ), MHTI( * ), MHLI( * ), MHMI( * ),
     1        ICMPI, MHTO( * ), MHLO( * ), MHMO( * ),
     2        ICMPO, ILNR, IRNR, N1, N2, IMF, KL, KU, KLE
      DOUBLE PRECISION FAC, A( N1, N2 )
C____________________________________________________________________C
C Variable declarations - Working variables .........................C
C
      INTEGER INDI, INDO, IHI, IHO, NRR, NH, IR,
     1        INDFUN, IROW, ICOL
      DOUBLE PRECISION TOL
      PARAMETER ( TOL = 1.0d-9 )
C____________________________________________________________________C
C START OF PROGRAM **************************************************C
C____________________________________________________________________C
C
C Check grid nodes ...
C
      NRR    = INARR( 2 )
      NH     = INARR( 3 )
C     .
C     . Easy exit
C     .
      IF ( ABS( FAC ).LT.TOL ) RETURN
C     .
      IF ( NRR.NE.NR ) THEN
         PRINT *,' Subroutine NGMTA.'
         PRINT *,' NR   = ', NR
         PRINT *,' NRR  = ', NRR
         PRINT *,' Program aborted.'
         STOP
      ENDIF
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
C
            DO IR = ILNR, IRNR
C
C Find locations of source and destination vectors ...
C
             INDI = INDFUN( IR, IHI, INARR )
             INDO = INDFUN( IR, IHO, INARR )
             CALL MATIND( INDO, INDI, IMF, KL, KU, KLE, N1, N2,
     1                    IROW, ICOL)
             A( IROW, ICOL ) = A( IROW, ICOL ) + FAC
C
            ENDDO
C
          ENDIF
        ENDDO
C
 50   CONTINUE
      ENDDO
C     .
      RETURN
      END
C*********************************************************************

