C*********************************************************************
C subroutine RADius VaLue Find ***************************************
C            ---    - -   -    ***************************************
C Steve Gibbons Thu Sep 16 13:03:18 BST 1999                         C
C____________________________________________________________________C
C                                                                    C
C  Gives the value of radius for a given grid node, IR.              C
C                                                                    C
C____________________________________________________________________C
C Input variables :-                                                 C
C ===============                                                    C
C  Integer                                                           C
C  -------                                                           C
C     IR        : Number of radial grid node.                        C
C     INARR     : Array of integers. Dimension ( * )                 C
C                 At present the elements of INARR are as follows.   C
C                                                                    C
C                 INARR( 1 ) = IFORMF. Format flag.                  C
C                  This flag decides how the elements are arranged   C
C                  in the solution vector.                           C
C                                                                    C
C                 INARR( 2 ) = NR. Number of radial grid nodes.      C
C                 INARR( 3 ) = NH.                                   C
C                                                                    C
C INARR must be the same array which is passed to INDFUN.            C
C                                                                    C
C If IFORMF = 1 or IFORMF = 2 then                                   C
C            RAD    = RI + DBLE( IR - 1 )*H                          C
C  where     H      = ( RO - RI )/DBLE( NR - 1 )                     C
C                                                                    C
C If IFORMF = 3 or IFORMF = 4 then RAD = XARR( IR )                  C
C                                                                    C
C  Double Precision                                                  C
C  ----------------                                                  C
C                                                                    C
C     RAD       : Returned as the value of the radius.               C
C                                                                    C
C     DPARR     : Array of dimension ( * ).                          C
C                                                                    C
C      DPARR( 1 ) = RI                                               C
C      DPARR( 2 ) = RO                                               C
C                                                                    C
C     XARR      : Array of dimension ( * ).                          C
C                 May contain coefficients for abscissae in some     C
C                  non-uniform differencing scheme.                  C
C                                                                    C
C     H         : Distance between adjacent grid nodes.              C
C                 Output variable only as it is calculated in situ.  C
C                 Only valid for uniform mesh schemes.               C
C                                                                    C
C____________________________________________________________________C
C
C*********************************************************************
      SUBROUTINE RADVLF ( RAD, IR, INARR, DPARR, XARR, H )
      IMPLICIT NONE
C____________________________________________________________________C
C Variable declarations - Parameters ................................C
      INTEGER IR, INARR( * )
      DOUBLE PRECISION RAD, DPARR( * ), XARR( * ), H
C____________________________________________________________________C
C Variable declarations - Working variables .........................C
      INTEGER IFORMF, NR
      DOUBLE PRECISION RI, RO, TOL
      PARAMETER ( TOL = 1.0d-6 )
C____________________________________________________________________C
C START OF PROGRAM **************************************************C
C____________________________________________________________________C
C
      H      = 0.0d0
C
      RI     = DPARR( 1 )
      RO     = DPARR( 2 )
C
      IFORMF = INARR( 1 )
      NR     = INARR( 2 )
C
      IF ( IR.LT.1 .OR. IR.GT.NR ) THEN
        PRINT *,' Function RADVLF. NR = ', NR
        PRINT *,' IR = ', IR,' Program aborted.'
        STOP
      ENDIF
C
      IF ( NR.LT.2 ) THEN
        PRINT *,' Function RADVLF.'
        PRINT *,' NR = ', NR,' Program aborted.'
        STOP
      ENDIF
C
      IF ( IFORMF.EQ.1 .OR. IFORMF.EQ.2 ) THEN
        H = ( RO - RI )/DBLE( NR - 1 )
        IF ( H.LT.TOL ) THEN
          PRINT *,' Function RADVLF.'
          PRINT *,' H = ', H,' Program aborted.'
          STOP
        ENDIF
        RAD = RI + DBLE( IR - 1 )*H
        RETURN
      ENDIF
C
      IF ( IFORMF.EQ.3 .OR. IFORMF.EQ.4 ) THEN
        RAD = XARR( IR )
        RETURN
      ENDIF
C
      PRINT *,' Function RADVLF. IFORMF = ', IFORMF
      PRINT *,' Not current option. Program aborted.'
      STOP
      END
C*********************************************************************
