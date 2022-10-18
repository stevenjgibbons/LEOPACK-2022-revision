C*********************************************************************
C subroutine INTegration Coefficient FinD ****************************
C            ---         -           -  - ****************************
C Steve Gibbons Fri Sep 24 13:52:53 BST 1999                         C
C____________________________________________________________________C
C                                                                    C
C If a function, f, of x is defined over an interval between         C
C x_{ilmn} and x_{irmn}, with x_j = XARR( j ) increasing             C
C monotonically, INTCFD returns a double precision array of          C
C dim ( NR ) such that                                               C
C                                                                    C
C \int_{x_[ilmn]}^{x_[irmn]} = \sum_{j=ilmn}^{j=irmn}                C
C                                         ( f( x_j )*DPINT( j ) )    C
C                                                                    C
C to the very roughest approximation.                                C
C____________________________________________________________________C
C Input variables :-                                                 C
C ===============                                                    C
C                                                                    C
C  Integer                                                           C
C  -------                                                           C
C     NR        : Number of radial grid nodes.                       C
C     ILMN      : Left most node to be used.                         C
C     IRMN      : Right most node to be used.                        C
C                                                                    C
C  Double Precision                                                  C
C  ----------------                                                  C
C     XARR      : Array containing the abscissa values. Dim ( NR )   C
C                                                                    C
C Output variables :-                                                C
C ================                                                   C
C                                                                    C
C     DPINTC    : Integration coeffcients. Dim ( NR )                C
C                                                                    C
C____________________________________________________________________C
C
C*********************************************************************
      SUBROUTINE INTCFD ( NR, XARR, DPINTC, ILMN, IRMN )
      IMPLICIT NONE
C____________________________________________________________________C
C Variable declarations - Parameters ................................C
      INTEGER NR
      DOUBLE PRECISION XARR( NR ), DPINTC( NR )
C____________________________________________________________________C
C Variable declarations - Working variables .........................C
      INTEGER IR, ILMN, IRMN
      DOUBLE PRECISION HI, ZERO
      PARAMETER ( ZERO = 0.0d0 )
C____________________________________________________________________C
C START OF PROGRAM **************************************************C
C____________________________________________________________________C
C Check validity of variables
C
      IF ( ILMN.LT.1 .OR. IRMN.GT.NR .OR. IRMN.LE.ILMN ) THEN
        PRINT *,' Subroutine INTCFD.'
        PRINT *,' ILMN = ', ILMN
        PRINT *,' IRMN = ', IRMN
        PRINT *,' NR   = ', NR
        PRINT *,' Program aborted.'
        STOP
      ENDIF
C
C First, simply clear all elements of DPINTC
C
      IR = 0
      CALL VECOP( DPINTC, ZERO, NR, IR )
C
      IR = ILMN
      HI = XARR( ILMN + 1 ) - XARR( ILMN )
      DPINTC( IR ) = 0.5d0*HI
C
      IR = IRMN
      HI = XARR( IRMN ) - XARR( IRMN - 1 )
      DPINTC( IR ) = 0.5d0*HI
C
C Now loop around the central points ...
C
      DO IR = ILMN + 1, IRMN - 1
        HI = XARR( IR + 1 ) - XARR( IR - 1 )
        DPINTC( IR ) = 0.5d0*HI
      ENDDO
C
      RETURN
      END
C*********************************************************************
