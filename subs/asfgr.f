C*********************************************************************
C subroutine Axi-Symmetric Flow Generation Routine *******************
C            -   -         -    -          -       *******************
C Steve Gibbons Tue Dec  7 12:16:56 GMT 1999                         C
C____________________________________________________________________C
C                                                                    C
C Fills a vector with an axi-symmetric velocity.                     C
C                                                                    C
C Not only are the radial functions filled, so are the               C
C index arrays, MHT, MHL and MHM. The radial functions are filled    C
C from node 1 to NR - MHP must be set after calling this routine.    C
C NHMAX is the maximum number of harmonics which may be set.         C
C                                                                    C
C The exact nature of the flow depends upon the integer flag IMODE.  C
C If IMODE = 1, the flow is simply set as a solid body rotation      C
C                                                                    C
C  i.e. v_{\phi} = s                                                 C
C                                                                    C
C If IMODE = 2, the flow is equatorially symmetric and depends       C
C upon the parameters NN1, NN2, ICS, RI and DPNF which are input     C
C in the arrays DPRARR and INTARR (see ZVDF also).                   C
C                                                                    C
C In the case IMODE = 2, the function ZVDF is called.                C
C                                                                    C
C____________________________________________________________________C
C                                                                    C
C Input variables :-                                                 C
C ===============                                                    C
C                                                                    C
C  Integer                                                           C
C  -------                                                           C
C                                                                    C
C     INARR     : See INDFUN. On input, only elements 1 and 2 need   C
C                 be filled with IFORMF and NR respectively.         C
C                                                                    C
C                 If the call is successful then on output, INARR(3) C
C                 is filled with NH.                                 C
C                                                                    C
C     MHT       : Type of harmonic formed. Only MHT( i ) = 2 is      C
C                 likely to result from this routine - toroidal vel. C
C                                                                    C
C     MHL       : Spherical harmonic degree l.                       C
C     MHM       : Spherical harmonic order, m. (Set to zero).        C
C                                                                    C
C     NHMAX     : The limit upon the number of harmonics in VEC.     C
C                                                                    C
C     LH        : Maximum spherical harmonic degree, l.              C
C                                                                    C
C     IMODE     : See above.                                         C
C                                                                    C
C     NTHPTS    : Number of points in theta. Chosen by ONTPPF        C
C                                                                    C
C     INTARR    : Integer parameter array. Dimension ( * )           C
C                 In the case of IMODE = 2,                          C
C                                                                    C
C                   INTARR( 1 ) = NN1 in above equation.             C
C                   INTARR( 2 ) = NN2 in above equation.             C
C                   INTARR( 3 ) = 1 for COS version                  C
C                                 2 for SIN version                  C
C                                                                    C
C  Double Precision                                                  C
C  ----------------                                                  C
C                                                                    C
C     DPRARR    : Double precision parameter array. Dimension ( * )  C
C                 In the case of IMODE = 2,                          C
C                                                                    C
C                   DPRARR( 1 ) = RI (inner radius value)            C
C                   DPRARR( 2 ) = DPNF (normalising factor).         C
C                                                                    C
C     PA        : Dim ( ( LH + 1 )*( LH + 2 )/2 , NTHPTS )           C
C     DPA       : Dim ( ( LH + 1 )*( LH + 2 )/2 , NTHPTS )           C
C                                                                    C
C PA and DPA are the associated Legendre functions and their         C
C derivatives and are calculated by SCHNLA.                          C
C                                                                    C
C     GAUX      : Dim ( NTHPTS ). Gauss points ( GAUWTS )            C
C     GAUW      : Dim ( NTHPTS ). Gauss weignts ( GAUWTS )           C
C                                                                    C
C     XARR      : Dim ( NR ). Values of radius.                      C
C                                                                    C
C____________________________________________________________________C
C                                                                    C
C Working Arrays   :-                                                C
C ================                                                   C
C  Double Precision                                                  C
C  ----------------                                                  C
C                                                                    C
C     QST       : Obvious ... ( dimension  (  LH*( LH+2 ) , 3 )      C
C     VF        : Vector function. dimensions ( 2, NTHPTS, 3)        C
C  (Note that as our flow is axisymmetric, we only need two points   C
C   in phi - this will save alot of time).                           C
C                                                                    C
C____________________________________________________________________C
C
C*********************************************************************
      SUBROUTINE ASFGR( VEC, INARR, MHT, MHL, MHM, NHMAX, LH, IMODE,
     1                  NTHPTS, DPRARR, INTARR, QST, PA, DPA, GAUX,
     2                  GAUW, VF, XARR )
      IMPLICIT NONE
C____________________________________________________________________C
C Variable declarations - Parameters ................................C
      INTEGER INARR( * ), MHT( * ), MHL( * ), MHM( * ), NHMAX, LH,
     1        IMODE, NTHPTS, INTARR( * )
      DOUBLE PRECISION VEC( * ), DPRARR( * ), QST( LH*LH + 2*LH, 3 ),
     1                 PA( ( LH + 1 )*( LH + 2 )/2 , NTHPTS),
     2                 DPA ( ( LH + 1 )*( LH + 2 )/2 , NTHPTS),
     3                 GAUX( NTHPTS ), GAUW( NTHPTS )
      DOUBLE PRECISION VF( 2, NTHPTS, 3), XARR( * )
C____________________________________________________________________C
C Variable declarations - Working variables .........................C
      INTEGER NPHPTS, NH, IH, NR, IR, IPHI, ITHETA, MMAX, L, M, ICS,
     1        IHARM, INDFUN, INDSHC, IND
      PARAMETER ( NPHPTS = 2, MMAX = 0 )
      DOUBLE PRECISION FTF1( 2*NPHPTS ), ZCOEF, FAC,
     1                 FTF2( 2*NPHPTS ), ZVDF, X, THE, RAD,
     2                 FTF3( 2*NPHPTS ), SQRLL1
C____________________________________________________________________C
C START OF PROGRAM **************************************************C
C____________________________________________________________________C
C Sort out invalid value of IMODE
C
      IF ( IMODE.EQ.1 ) GOTO 50
      IF ( IMODE.EQ.2 ) GOTO 50
      PRINT *,' Subroutine ASFGR.'
      PRINT *,' IMODE = ', IMODE
      PRINT *,' Program aborted.'
      STOP
C
 50   CONTINUE
      NR = INARR( 2 )
      IF ( NR.LT.10 ) THEN
        PRINT *,' Subroutine ASFGR.'
        PRINT *,' NR = ', NR
        PRINT *,' Program aborted.'
        STOP
      ENDIF
C
C Simple solid body rotation case ...
C
      IF ( IMODE.EQ.1 ) THEN
        NH = 1
        INARR( 3 ) = NH
        IF ( NH.GT.NHMAX ) THEN
          PRINT *,' Subroutine ASFGR.'
          PRINT *,' NH = ', NH,' NHMAX = ', NHMAX
          PRINT *,' Program aborted.'
          STOP
        ENDIF
        IH = 1
        MHT( IH ) = 2
        MHL( IH ) = 1
        MHM( IH ) = 0
        DO IR = 1, NR
          RAD = XARR( IR )
          IND = INDFUN( IR, IH, INARR )
          VEC( IND ) = RAD
        ENDDO
        INARR( 3 ) = NH
      ENDIF
C
C Jupiter type flow
C
      IF ( IMODE.EQ.2 ) THEN
        IF ( LH/2*2.EQ.LH ) THEN
          NH = LH/2
        ELSE
          NH = LH/2+1
        ENDIF
        INARR( 3 ) = NH
        IF ( NH.GT.NHMAX ) THEN
          PRINT *,' Subroutine ASFGR.'
          PRINT *,' NH = ', NH,' NHMAX = ', NHMAX
          PRINT *,' Program aborted.'
          STOP
        ENDIF
        DO IH = 1, NH
          MHT( IH ) = 2
          MHL( IH ) = 2*IH - 1
          MHM( IH ) = 0
        ENDDO
C
C Now let's fill in the values of VEC
C
        DO IR = 1, NR
          RAD = XARR( IR )
C         .
C         . Now fill the array VF
C         .
          DO ITHETA = 1, NTHPTS
            X   = GAUX( ITHETA )
            THE = ACOS( X )
            DO IPHI = 1, NPHPTS
              VF( IPHI, ITHETA, 1 ) = 0.0d0
              VF( IPHI, ITHETA, 2 ) = 0.0d0
              VF( IPHI, ITHETA, 3 ) = 
     1          ZVDF ( RAD, THE, DPRARR, INTARR )
            ENDDO
          ENDDO
C         .
C         . VF now contains velocity evaluated at this node.
C         . Now transform to QST coefficients.
C         .
          CALL VF2QST ( QST, VF, GAUX, GAUW, PA, DPA, FTF1, FTF2,
     1                  FTF3, ZCOEF, LH, NTHPTS, NPHPTS, MMAX )
C         .
C         . Now convert our QST coefs to toroidal
C         .
          M   = 0
          ICS = 1
          DO IH = 1, NH
            L     = MHL( IH )
            IHARM = INDSHC( L, M, ICS )
            IND   = INDFUN( IR, IH, INARR )
            FAC   = (-1.0d0)/SQRLL1( L )
            VEC( IND ) = FAC*QST( IHARM, 3 )
          ENDDO
C         .
        ENDDO
C
      ENDIF
C
      RETURN
      END
C*********************************************************************
