C*********************************************************************
C subroutine Solution Vector DERivative Find *************************
C            -        -      ---        -    *************************
C Steve Gibbons 16.3.99                                              C
C____________________________________________________________________C
C                                                                    C
C Calculates the derivative at any grid node for a specified         C
C harmonic in the solution vector. This new routine is entirely      C
C general and makes no reference to boundary conditions etc.         C
C____________________________________________________________________C
C                                                                    C
C Input variables :-                                                 C
C ===============                                                    C
C  Integer                                                           C
C  -------                                                           C
C                                                                    C
C     NDIM      : Length of the solution vector.                     C
C     NH        : Number of harmonics.                               C
C     IH        : Number of harm. at which deriv. is to be evaluated C
C     NR        : Number of radial grid nodes.                       C
C     IRN       : Number of node at which deriv. is to be evaluated. C
C                                                                    C
C  Double Precision                                                  C
C  ----------------                                                  C
C     SV        : Solution vector - dimension ( NDIM )               C
C     RI        : Radius of inner boundary.                          C
C     RO        : Radius of outer boundary.                          C
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
      SUBROUTINE SVDERF ( NDIM, NH, IH, NR, IRN, SV, ORD, 
     1                    D0F, D1F, D2F, D3F, D4F, RI, RO )
      IMPLICIT NONE
C____________________________________________________________________C
C Variable declarations - Parameters ................................C
      INTEGER NDIM, NH, IH, NR, IRN
      DOUBLE PRECISION D0F, D1F, D2F, D3F, D4F, SV( NDIM ),
     1                 RI, RO
      CHARACTER *(2) ORD
C____________________________________________________________________C
C Variable declarations - Working variables .........................C
      DOUBLE PRECISION FM4, FM3, FM2, FM1, F00, FP1, FP2, FP3, FP4,
     1                 CM4( 4 ), CM3( 4 ), CM2( 4 ), CM1( 4 ),
     2                 C00( 4 ), CP1( 4 ), CP2( 4 ), CP3( 4 ),
     3                 CP4( 4 ), H
      INTEGER IM4, IM3, IM2, IM1, I00, IP1, IP2, IP3, IP4,
     1        ILM, IRM
      LOGICAL LB
C____________________________________________________________________C
C START OF PROGRAM **************************************************C
C____________________________________________________________________C
C
C Check the length of the vector
C
      IF ( NDIM.NE.NR*NH ) THEN
         PRINT *,' Subroutine SVDERF. Bad dimensions '
         PRINT *,' Array length = ',NDIM
         PRINT *,' NR = ', NR,' NH = ',NH
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
         I00 = ( IRN - 1 )*NH + IH
         F00 = SV( I00 )
         IP1 = ( IRN     )*NH + IH
         FP1 = SV( IP1 )
         IP2 = ( IRN + 1 )*NH + IH
         FP2 = SV( IP2 )
         IP3 = ( IRN + 2 )*NH + IH
         FP3 = SV( IP3 )
         IP4 = ( IRN + 3 )*NH + IH
         FP4 = SV( IP4 )
         CALL GSLDCF ( H, C00, CP1, CP2, CP3, CP4 )
         D0F = F00
         CALL DDERC5 ( C00, F00, CP1, FP1, CP2, FP2,
     1                 CP3, FP3, CP4, FP4, D1F, D2F, D3F, D4F )
         RETURN
      ENDIF
C
      IF ( IRN.EQ.2 ) THEN
         IM1 = ( IRN - 2 )*NH + IH
         FM1 = SV( IM1 )
         I00 = ( IRN - 1 )*NH + IH
         F00 = SV( I00 )
         IP1 = ( IRN     )*NH + IH
         FP1 = SV( IP1 )
         IP2 = ( IRN + 1 )*NH + IH
         FP2 = SV( IP2 )
         IP3 = ( IRN + 2 )*NH + IH
         FP3 = SV( IP3 )
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
         IM2 = ( IRN - 3 )*NH + IH
         FM2 = SV( IM2 )
         IM1 = ( IRN - 2 )*NH + IH
         FM1 = SV( IM1 )
         I00 = ( IRN - 1 )*NH + IH
         F00 = SV( I00 )
         IP1 = ( IRN     )*NH + IH
         FP1 = SV( IP1 )
         IP2 = ( IRN + 1 )*NH + IH
         FP2 = SV( IP2 )
         CALL CENDIF( H, 'O5', CM3, CM2, CM1, C00, CP1, CP2, CP3 )
         D0F = F00
         CALL DDERC5 ( CM2, FM2, CM1, FM1, C00, F00, CP1, FP1,
     1                 CP2, FP2, D1F, D2F, D3F, D4F )
         RETURN
      ENDIF
C
      IF ( IRN.GT.ILM .AND. IRN.LT.IRM .AND. LB ) THEN
         IM3 = ( IRN - 4 )*NH + IH
         FM3 = SV( IM3 )
         IM2 = ( IRN - 3 )*NH + IH
         FM2 = SV( IM2 )
         IM1 = ( IRN - 2 )*NH + IH
         FM1 = SV( IM1 )
         I00 = ( IRN - 1 )*NH + IH
         F00 = SV( I00 )
         IP1 = ( IRN     )*NH + IH
         FP1 = SV( IP1 )
         IP2 = ( IRN + 1 )*NH + IH
         FP2 = SV( IP2 )
         IP3 = ( IRN + 2 )*NH + IH
         FP3 = SV( IP3 )
         CALL CENDIF( H, ORD, CM3, CM2, CM1, C00, CP1, CP2, CP3 )
         D0F = F00
         CALL DDERC7 ( CM3, FM3, CM2, FM2, CM1, FM1, C00, F00,
     1                 CP1, FP1, CP2, FP2, CP3, FP3,
     2                 D1F, D2F, D3F, D4F )
         RETURN
      ENDIF
C
      IF ( IRN.EQ.(NR-1) ) THEN
         IM3 = ( IRN - 4 )*NH + IH
         FM3 = SV( IM3 )
         IM2 = ( IRN - 3 )*NH + IH
         FM2 = SV( IM2 )
         IM1 = ( IRN - 2 )*NH + IH
         FM1 = SV( IM1 )
         I00 = ( IRN - 1 )*NH + IH
         F00 = SV( I00 )
         IP1 = ( IRN     )*NH + IH
         FP1 = SV( IP1 )
         CALL GFRDCF ( H, CM3, CM2, CM1, C00, CP1 )
         D0F = F00
         CALL DDERC5 ( CM3, FM3, CM2, FM2, CM1, FM1,
     1                 C00, F00, CP1, FP1, D1F, D2F, D3F, D4F )
         RETURN
      ENDIF
C
      IF ( IRN.EQ.NR ) THEN
         IM4 = ( IRN - 5 )*NH + IH
         FM4 = SV( IM4 )
         IM3 = ( IRN - 4 )*NH + IH
         FM3 = SV( IM3 )
         IM2 = ( IRN - 3 )*NH + IH
         FM2 = SV( IM2 )
         IM1 = ( IRN - 2 )*NH + IH
         FM1 = SV( IM1 )
         I00 = ( IRN - 1 )*NH + IH
         F00 = SV( I00 )
         CALL GSRDCF ( H, CM4, CM3, CM2, CM1, C00 )
         D0F = F00
         CALL DDERC5 ( CM4, FM4, CM3, FM3, CM2, FM2, CM1, FM1, 
     1                 C00, F00, D1F, D2F, D3F, D4F )
         RETURN
      ENDIF
C
C____________________________________________________________________C
      PRINT *,' Error in SVDERF. Program aborted.'
      STOP
      END
C*********************************************************************
