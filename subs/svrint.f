C*********************************************************************
C subroutine Solution Vector Radial INTerpolate **********************
C            -        -      -      ---         **********************
C Steve Gibbons Sun Oct 31 14:12:34 GMT 1999                         C
C____________________________________________________________________C
C                                                                    C
C Evaluates the zero^{th} to the (NNDS-1)^{th} derivative of a       C
C solution vector (harmonic IH) at the arbitrary radial point, R.    C
C                                                                    C
C  The double precision array element VALS( i ) is returned          C
C  with f^{i-1}( r ) where f is the function of the IH^{th}          C
C  harmonic in SV.                                                   C
C                                                                    C
C Note that ALL points may be referred to and so if solution         C
C vector contains points which are stored implicitly, these must     C
C be filled in by a call to ASVCPL.                                  C
C                                                                    C
C____________________________________________________________________C
C                                                                    C
C Input variables :-                                                 C
C ===============                                                    C
C  Integer                                                           C
C  -------                                                           C
C     INARR     : Integer array dimension ( * ).                     C
C                  Elements may be arbitrary except for              C
C                  INARR( 1 ) = IFORMF - flag for vector format.     C
C                                        See INDFUN                  C
C                  INARR( 2 ) = NR. Number of radial grid nodes.     C
C                  INARR( 3 ) = NH = total number of radial func.s   C
C     IH        : Number of radial function (harmonic).              C
C     NNDS      : Number of nodes to be used in interpolation.       C
C                 Must be atleast 2.                                 C
C     IWORK     : Dimension ( NNDS ). Work array.                    C
C                                                                    C
C  Double Precision                                                  C
C  ----------------                                                  C
C                                                                    C
C     SV        : Solution vector. Dim ( * ) but length atleast      C
C                  NR*NH.                                            C
C     XARR      : Dim ( * ) but length atleast NR. Location of       C
C                  radial grid nodes.                                C
C                                                                    C
C     VALS      : Returned containing the {i-1}^{th} derivative      C
C                  at r=R in VALS( i )                               C
C                                                                    C
C     WORK      : Dimension ( NNDS ). Work array.                    C
C     COEFM     : Dimension ( NNDS, NNDS ). Work array.              C
C                                                                    C
C____________________________________________________________________C
C
C*********************************************************************
      SUBROUTINE SVRINT( R, SV, XARR, INARR, IH, NNDS, VALS, IWORK,
     1                   WORK, COEFM )
      IMPLICIT NONE
C____________________________________________________________________C
C Variable declarations - Parameters ................................C
      INTEGER INARR( * ), NNDS, IWORK( NNDS ), IH
      DOUBLE PRECISION R, SV( * ), XARR( * ), VALS( NNDS ),
     1                 WORK( NNDS ), COEFM( NNDS, NNDS )
C____________________________________________________________________C
C Variable declarations - Working variables .........................C
      INTEGER NR, IR, IRAD, IFN, ILN, IND, I, INC, INDFUN, NH
      DOUBLE PRECISION RI, RO, ZERO, ONE, RAD, DTOL
      CHARACTER *(1) TRANS
      PARAMETER ( ZERO = 0.0d0, ONE = 1.0d0, TRANS = 'N',
     1            DTOL = 2.0d-6, INC = 1 )
C____________________________________________________________________C
C START OF PROGRAM **************************************************C
C____________________________________________________________________C
C
C Check input parameters ...
C
      NR     = INARR( 2 )
      NH     = INARR( 3 )
C
      IF ( IH.LT.1 .OR. IH.GT.NH ) THEN
        PRINT *,' Subroutine SVRINT.'
        PRINT *,' IH = ',IH
        PRINT *,' NH = ',NH
        PRINT *,' Program aborted.'
        STOP
      ENDIF
C
      RI  = XARR(  1 )
      RO  = XARR( NR )
      RAD = R
C
      IF ( RAD.LT.RI .OR. RAD.GT.RO ) THEN
C       .
C       . Ammendment. Here it is possible
C       . that RAD is only less than RI or
C       . greater than RO because of numerical roundoff.
C       . For example:
C             RAD  =   0.6666666666000000    
C             RI   =   0.6666666700000000 
C       . To rectify this, we will try and set a 
C       . new condition on RAD
C       .
        IF ( DABS( RAD - RI ).LT.DTOL ) THEN
          RAD = RI
          GOTO 51
        ENDIF
C       .
        IF ( DABS( RAD - RO ).LT.DTOL ) THEN
          RAD = RO
          GOTO 51
        ENDIF
C       .
        PRINT *,' Subroutine SVRINT.'
        PRINT *,' RAD  = ',RAD
        PRINT *,' RI   = ',RI
        PRINT *,' RO   = ',RO
        PRINT *,' Program aborted.'
        STOP
      ENDIF
 51   CONTINUE
C
      IF ( NNDS.LT.2 ) THEN
        PRINT *,' Subroutine SVRINT.'
        PRINT *,' NNDS = ', NNDS
        PRINT *,' Program aborted.'
        STOP
      ENDIF
C
      IRAD = NR
      DO IR = 1, NR - 1
        IF ( XARR( IR ).LE.RAD .AND. XARR( IR + 1 ).GT.RAD ) THEN
          IRAD = IR
          GOTO 50
        ENDIF
      ENDDO
C
 50   CONTINUE
      IFN = IRAD + 1 - (NNDS+1)/2
      ILN = IRAD + NNDS/2
C
      DO IR = 1, NNDS
        IWORK( IR ) = IFN + IR - 1
        IF ( IWORK( IR ).LT.1 ) IWORK( IR ) = ILN - IWORK( IR ) + 1
        IF ( IWORK( IR ).GT.NR ) IWORK( IR ) = NR + IFN - IWORK( IR )
        I = IWORK( IR )
        VALS( IR ) = XARR( I )
      ENDDO
C
C Calculate finite difference coefficients
C
      CALL GFDCFD( RAD, VALS, NNDS, COEFM, NNDS, IWORK, WORK )
C
      DO IR = 1, NNDS
        IWORK( IR ) = IFN + IR - 1
        IF ( IWORK( IR ).LT.1 ) IWORK( IR ) = ILN - IWORK( IR ) + 1
        IF ( IWORK( IR ).GT.NR ) IWORK( IR ) = NR + IFN - IWORK( IR )
        I = IWORK( IR )
        IND = INDFUN( I, IH, INARR )
        WORK( IR ) = SV( IND )
      ENDDO
C
C WORK now contains the function values at the nodes
C so we can now multiply COEFM by WORK to get VALS
C
      CALL DGEMV ( TRANS, NNDS, NNDS, ONE, COEFM, NNDS, WORK, INC,
     1             ZERO, VALS, INC )
C  
      RETURN
      END
C*********************************************************************
