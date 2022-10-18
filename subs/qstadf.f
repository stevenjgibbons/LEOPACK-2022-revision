C*********************************************************************
C subroutine QST Array Derivative Find *******************************
C            --- -     -          -    *******************************
C Steve Gibbons 16.3.99                                              C
C ( adapted from svderf  4. 8.99 )                                   C
C____________________________________________________________________C
C                                                                    C
C Calculates the derivative at any grid node for a specified         C
C harmonic in the full qst vector decomposition.                     C
C____________________________________________________________________C
C                                                                    C
C Input variables :-                                                 C
C ===============                                                    C
C  Integer                                                           C
C  -------                                                           C
C                                                                    C
C     LH        : Maximum spherical harmonic degree.                 C
C     IH        : Number of harm. at which deriv. is to be evaluated C
C     NR        : Number of radial grid nodes.                       C
C     IRN       : Number of node at which deriv. is to be evaluated. C
C     IQST      : =1 for Q, =2 for S and =3 for T.                   C
C                                                                    C
C  Double Precision                                                  C
C  ----------------                                                  C
C     QSTARR    : QST decomp. Dimension ( LH*(LH+2), 3, NR )         C
C     RI        : Radius of inner boundary.                          C
C     RO        : Radius of outer boundary.                          C
C                                                                    C
C  Character                                                         C
C  ---------                                                         C
C     ORD       : *(2). Determines the order of accuracy. Options :  C
C             SS - Strictly second order                             C
C             SF - Strictly fourth order                             C
C             O5 - Optimum accuracy for bandwidth 5; this gives      C
C                  Fourth order accuracy for 1st and 2nd derivatives C
C                  and second order accuracy for 3rd and 4th der.s   C
C             O7 - Optimum accuracy for bandwidth 7; this gives      C
C                  Sixth order accuracy for 1st and 2nd derivatives  C
C                  and fourth order accuracy for 3rd and 4th der.s   C
C                                                                    C
C____________________________________________________________________C
C Output variables :-                                                C
C ================                                                   C
C  Double Precision                                                  C
C  ----------------                                                  C
C     D0F       : Value of function.                                 C
C     D1F       : First derivative.                                  C
C     D2F       : Second derivative.                                 C
C     D3F       : Third derivative.                                  C
C     D4F       : Fourth derivative.                                 C
C____________________________________________________________________C
C
C*********************************************************************
      SUBROUTINE QSTADF ( LH, IH, NR, IRN, IQST, QSTARR, ORD, 
     1                    D0F, D1F, D2F, D3F, D4F, RI, RO )
      IMPLICIT NONE
C____________________________________________________________________C
C Variable declarations - Parameters ................................C
      INTEGER LH, IH, NR, IRN, IQST
      DOUBLE PRECISION D0F, D1F, D2F, D3F, D4F,
     1                 RI, RO, QSTARR( LH*(LH+2), IQST, NR )
      CHARACTER *(2) ORD
C____________________________________________________________________C
C Variable declarations - Working variables .........................C
      DOUBLE PRECISION FM4, FM3, FM2, FM1, F00, FP1, FP2, FP3, FP4,
     1                 CM4( 4 ), CM3( 4 ), CM2( 4 ), CM1( 4 ),
     2                 C00( 4 ), CP1( 4 ), CP2( 4 ), CP3( 4 ),
     3                 CP4( 4 ), H
      INTEGER ILM, IRM, NH
      LOGICAL LB
C____________________________________________________________________C
C START OF PROGRAM **************************************************C
C____________________________________________________________________C
C
      IF ( IQST.NE.1 .AND. IQST.NE.2 .AND. IQST.NE.3 ) THEN
        PRINT *,' Subroutine QSTADF, IQST = ', IQST
        STOP
      ENDIF
C
      NH = LH*(LH + 2 )
      IF ( IH.LT.1 .OR. IH.GT.NH ) THEN
        PRINT *,' Subroutine QSTADF, IH = ', IH
        PRINT *,' LH = ', LH ,' NH = ', NH
        STOP
      ENDIF
C
      ILM = 3
      IRM = NR - 2
C
      H = ( RO - RI )/DBLE(NR - 1)
C
      IF ( ORD.EQ.'SF' .OR. ORD.EQ.'O7' ) THEN
         LB = .TRUE.
      ELSE
         LB = .FALSE.
      ENDIF
C____________________________________________________________________C
C
      IF ( IRN.EQ.1 ) THEN
         F00 = QSTARR( IH, IQST, IRN )
         FP1 = QSTARR( IH, IQST, IRN + 1)
         FP2 = QSTARR( IH, IQST, IRN + 2)
         FP3 = QSTARR( IH, IQST, IRN + 3)
         FP4 = QSTARR( IH, IQST, IRN + 4)
         CALL GSLDCF ( H, C00, CP1, CP2, CP3, CP4 )
         D0F = F00
         CALL DDERC5 ( C00, F00, CP1, FP1, CP2, FP2,
     1                 CP3, FP3, CP4, FP4, D1F, D2F, D3F, D4F )
         RETURN
      ENDIF
C
      IF ( IRN.EQ.2 ) THEN
         FM1 = QSTARR( IH, IQST, IRN - 1)
         F00 = QSTARR( IH, IQST, IRN )
         FP1 = QSTARR( IH, IQST, IRN + 1)
         FP2 = QSTARR( IH, IQST, IRN + 2)
         FP3 = QSTARR( IH, IQST, IRN + 3)
         CALL GFLDCF ( H, CM1, C00, CP1, CP2, CP3 )
         D0F = F00
         CALL DDERC5 ( CM1, FM1, C00, F00, CP1, FP1,
     1                 CP2, FP2, CP3, FP3, D1F, D2F, D3F, D4F )
         RETURN
      ENDIF
C
      IF (  IRN.EQ.ILM .OR. IRN.EQ.IRM  .OR.
     1    ( IRN.GT.ILM .AND. IRN.LT.IRM .AND. ( .NOT. LB ) )
     2                                       ) THEN
         FM2 = QSTARR( IH, IQST, IRN - 2)
         FM1 = QSTARR( IH, IQST, IRN - 1)
         F00 = QSTARR( IH, IQST, IRN )
         FP1 = QSTARR( IH, IQST, IRN + 1)
         FP2 = QSTARR( IH, IQST, IRN + 2)
         CALL CENDIF( H, 'O5', CM3, CM2, CM1, C00, CP1, CP2, CP3 )
         D0F = F00
         CALL DDERC5 ( CM2, FM2, CM1, FM1, C00, F00, CP1, FP1,
     1                 CP2, FP2, D1F, D2F, D3F, D4F )
         RETURN
      ENDIF
C
      IF ( IRN.GT.ILM .AND. IRN.LT.IRM .AND. LB ) THEN
         FM3 = QSTARR( IH, IQST, IRN - 3)
         FM2 = QSTARR( IH, IQST, IRN - 2)
         FM1 = QSTARR( IH, IQST, IRN - 1)
         F00 = QSTARR( IH, IQST, IRN )
         FP1 = QSTARR( IH, IQST, IRN + 1)
         FP2 = QSTARR( IH, IQST, IRN + 2)
         FP3 = QSTARR( IH, IQST, IRN + 3)
         CALL CENDIF( H, ORD, CM3, CM2, CM1, C00, CP1, CP2, CP3 )
         D0F = F00
         CALL DDERC7 ( CM3, FM3, CM2, FM2, CM1, FM1, C00, F00,
     1                 CP1, FP1, CP2, FP2, CP3, FP3,
     2                 D1F, D2F, D3F, D4F )
         RETURN
      ENDIF
C
      IF ( IRN.EQ.(NR-1) ) THEN
         FM3 = QSTARR( IH, IQST, IRN - 3)
         FM2 = QSTARR( IH, IQST, IRN - 2)
         FM1 = QSTARR( IH, IQST, IRN - 1)
         F00 = QSTARR( IH, IQST, IRN )
         FP1 = QSTARR( IH, IQST, IRN + 1)
         CALL GFRDCF ( H, CM3, CM2, CM1, C00, CP1 )
         D0F = F00
         CALL DDERC5 ( CM3, FM3, CM2, FM2, CM1, FM1,
     1                 C00, F00, CP1, FP1, D1F, D2F, D3F, D4F )
         RETURN
      ENDIF
C
      IF ( IRN.EQ.NR ) THEN
         FM4 = QSTARR( IH, IQST, IRN - 4 )
         FM3 = QSTARR( IH, IQST, IRN - 3 )
         FM2 = QSTARR( IH, IQST, IRN - 2 )
         FM1 = QSTARR( IH, IQST, IRN - 1 )
         F00 = QSTARR( IH, IQST, IRN )
         CALL GSRDCF ( H, CM4, CM3, CM2, CM1, C00 )
         D0F = F00
         CALL DDERC5 ( CM4, FM4, CM3, FM3, CM2, FM2, CM1, FM1, 
     1                 C00, F00, D1F, D2F, D3F, D4F )
         RETURN
      ENDIF
C
C____________________________________________________________________C
      PRINT *,' Error in QSTADF. Program aborted.'
      STOP
      END
C*********************************************************************
