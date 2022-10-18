C*********************************************************************
C subroutine Solution Vector 2 Scalar Harmonic Function **************
C            -        -      - -      -        -        **************
C Steve Gibbons Sat Jun  3 11:11:22 BST 2000                         C
C____________________________________________________________________C
C                                                                    C
C SV is a solution indexed as normal by INARR, MHT, MHL and MHM.     C
C SV2SHF simply takes a grid node IR and fills the spherical         C
C harmonic vector SHC with the values in the appropriate entries.    C
C                                                                    C
C SV2SHF does NOT do radial derivatives.                             C
C If you wish to evaluate a scalar function of a radial derivative   C
C then simply use CASVDR (for example) and then supply V1 into SV.   C
C                                                                    C
C____________________________________________________________________C
C                                                                    C
C Input variables :-                                                 C
C ===============                                                    C
C  Integer                                                           C
C  -------                                                           C
C                                                                    C
C     IR        : Number of radial grid node.                        C
C     INARR     : Indexing array. Dim (3).                           C
C                                                                    C
C                 INARR( 1 ) = IFORMF                                C
C                 INARR( 2 ) = NR      See INDFUN for details        C
C                 INARR( 3 ) = NH                                    C
C                                                                    C
C     LH        : Maximum spherical harmonic degree, l.              C
C     MHT       : Array length ( * ) - atleast length NH             C
C                  See below for key. (corresponds to input vec.)    C
C                                                                    C
C     MHT ( i ) = 3 for a scalar harmonic (all others are ignored).  C
C                                                                    C
C     MHL       : Array length ( * ) - atleast length NH             C
C                  Sph. harm. degree, l.                             C
C     MHM       : Array length ( * ) - atleast length NH             C
C                  Sph. harm. order, m, for cos m phi dep.           C
C                 -Sph. harm. order, m, for sin m phi dep.           C
C                                                                    C
C  Double Precision                                                  C
C  ----------------                                                  C
C                                                                    C
C     SV        : Solution vector. Dim ( * ) ( input )               C
C                 Length must be atleast NR*NH                       C
C                                                                    C
C     ZCOEF     : Coefficient of monopole term.                      C
C     SHC       : Spherical harmonic coefficients. Dim ( LH*(LH+2)). C
C                                                                    C
C____________________________________________________________________C
C
C*********************************************************************
      SUBROUTINE SV2SHC( IR, INARR, LH, MHT, MHL, MHM, SV, ZCOEF, SHC)
      IMPLICIT NONE
C____________________________________________________________________C
C Variable declarations - Parameters ................................C
C
      INTEGER IR, INARR( * ), LH, MHT( * ), MHL( * ), MHM( * )
      DOUBLE PRECISION SV( * ), ZCOEF, SHC( LH*(LH+2) )
C____________________________________________________________________C
C Variable declarations - Working variables .........................C
      INTEGER IOP, NH, IH, L, M, ICS, IHARM, INDSHC, INDFUN, IND
      DOUBLE PRECISION ZERO, D0F
      PARAMETER ( ZERO = 0.0d0, IOP = 0 )
C____________________________________________________________________C
C START OF PROGRAM **************************************************C
C____________________________________________________________________C
C
      NH    = LH*(LH+2)
      ZCOEF = ZERO
      CALL VECOP( SHC, ZERO, NH, IOP )
C
      NH    = INARR( 3 )
      DO IH = 1, NH
        IF (    MHT( IH ).NE.3     )    GOTO 50
        L  = MHL( IH )
        IF ( L.GT.LH ) THEN
          PRINT *,' Subroutine SV2SHC.'
          PRINT *,' Harmonic ',IH,' has L = ',L
          PRINT *,' LH = ', LH
          PRINT *,' Program aborted.'
          STOP
        ENDIF
        IF (          MHM( IH ).GE.0      ) THEN
          M =  MHM( IH )
          ICS = 1
        ELSE
          M = -MHM( IH )
          ICS = 2
        ENDIF
        IHARM = INDSHC( L, M, ICS )
C       .
        IND = INDFUN( IR, IH, INARR )
        D0F = SV( IND )
C       .
C       . Add to either ZCOEF or SHC
C       .
        IF ( IHARM.EQ.0 ) THEN
          ZCOEF = D0F
        ELSE
          SHC( IHARM )  = D0F
        ENDIF
C       .
 50   CONTINUE
      ENDDO
C     .
      RETURN
      END
C*********************************************************************

