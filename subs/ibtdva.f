C*********************************************************************
C subroutine Inhomogeneous Boundary Temp. Deriv.s Vector Addition. ***
C            -             -        -     -       -      -         ***
C Steve Gibbons Sat Feb  5 09:20:32 GMT 2000                         C
C____________________________________________________________________C
C                                                                    C
C If V0, V1, V2, V3 and V4 are all the derivative vectors output     C
C by CASVDR (from the homogeneous temperature functions), then       C
C IBTDVA will add on to these values the contributions from the      C
C inhomogeneous temperatures.                                        C
C                                                                    C
C____________________________________________________________________C
C                                                                    C
C The function f( r ) for ri .le. r .le. ro is defined by            C
C                                                                    C
C f( r ) =            CA sin[ pi/2 (r-ri)/(ro-ri) ]                  C
C            + CB 2 ( ri-ro )/pi cos[ pi/2 (r-ri)/(ro-ri) ]  +  CC   C
C                                                                    C
C____________________________________________________________________C
C                                                                    C
C Input variables :-                                                 C
C ===============                                                    C
C  Integer                                                           C
C  -------                                                           C
C                                                                    C
C     ILN       : Left-most node to be acted upon.                   C
C                                                                    C
C     IRN       : Right-most node to be acted upon.                  C
C                                                                    C
C     IHD       : Highest derivative requested.                      C
C                                                                    C
C     INARR     : Integer array dimension ( * ).                     C
C                  Elements may be arbitrary except for              C
C                  INARR( 1 ) = IFORMF - flag for vector format.     C
C                                        See INDFUN                  C
C                  INARR( 2 ) = NR. Number of radial grid nodes.     C
C                  INARR( 3 ) = NH = total number of radial func.s   C
C                                                                    C
C     MHT       : Type array. MHT( ih ) = 3 --> temperature harm.    C
C                                                                    C
C     MHI      : Dim ( * ) - atleast atleast length NHI.             C
C                For each temperature harmonic, IHI, MHI(IHI) gives  C
C                the index of the array CAFIT which stores the       C
C                coefficients for that radial function.              C
C                                                                    C
C f( r ) =            CA sin[ pi/2 (r-ri)/(ro-ri) ]                  C
C            + CB 2 ( ri-ro )/pi cos[ pi/2 (r-ri)/(ro-ri) ]  +  CC   C
C                                                                    C
C                If MHI( IH ) = IITH, then CA, CB and CC are         C
C                respectively stored in CAFIT( 1, IITH ),            C
C                CAFIT( 2, IITH ) and CAFIT( 3, IITH ).              C
C                                                                    C
C  Double Precision                                                  C
C  ----------------                                                  C
C                                                                    C
C     XARR      : Dim ( NR ). Radial values.                         C
C                                                                    C
C     CAFIT     : Dimension ( 3, * ). See MHI.                       C
C                                                                    C
C     V0        : Zero^th derivatives. Dim ( * ) but length atleast  C
C     V1        : First   derivatives. Dim ( * ) but length atleast  C
C     V2        : Second  derivatives. Dim ( * ) but length atleast  C
C     V3        : Third   derivatives. Dim ( * ) but length atleast  C
C     V4        : Fourth  derivatives. Dim ( * ) but length atleast  C
C                                                                    C
C IMPORTANT: ALL of V0, V1, V2, V3 and V4 are filled even if         C
C these derivatives are not requested!                               C
C                                                                    C
C                                                                    C
C____________________________________________________________________C
C
C*********************************************************************
      SUBROUTINE IBTDVA( ILN, IRN, IHD, INARR, MHT, MHI, XARR,
     1                   CAFIT, V0, V1, V2, V3, V4 )
      IMPLICIT NONE
C____________________________________________________________________C
C Variable declarations - Parameters ................................C
      INTEGER ILN, IRN, IHD, INARR( * ), MHT( * ), MHI( * )
      DOUBLE PRECISION XARR( * ), CAFIT( 3, * ), V0( * ),
     1                 V1( * ), V2( * ), V3( * ), V4( * )
C____________________________________________________________________C
C Variable declarations - Working variables .........................C
C
      INTEGER NR, NH, IH, IR, I, IITH, IND, INDFUN
      DOUBLE PRECISION RAD, RI, RO, CA, CB, CC, DERV( 5 ), DLOW
      PARAMETER ( DLOW = 1.0d-10 )
C____________________________________________________________________C
C START OF PROGRAM **************************************************C
C____________________________________________________________________C
C
      NR = INARR( 2 )
      NH = INARR( 3 )
C
      RI = XARR(  1 )
      RO = XARR( NR )
C
      IF ( IHD.LT.0 .OR. IHD.GT.4 ) THEN
        PRINT *,' Subroutine IBTDVA. IHD = ', IHD
        PRINT *,' Program aborted.'
        STOP
      ENDIF
C
      DO IH = 1, NH
        IF ( MHT( IH ).NE.3 ) GOTO 60
        IITH   = MHI( IH )
        CA     = CAFIT( 1, IITH )
        CB     = CAFIT( 2, IITH )
        CC     = CAFIT( 3, IITH )
        IF ( DABS( CA ).LT.DLOW .AND. DABS( CB ).LT.DLOW .AND.
     1       DABS( CC ).LT.DLOW ) GOTO 60
        DO IR = ILN, IRN
          RAD = XARR( IR )
C         .
          IND = INDFUN( IR, IH, INARR )
          DO I = 1, 5
            DERV( I ) = 0.0d0
          ENDDO
C
          CALL ITFA( RAD, RI, RO, CA, CB, CC, DERV, IHD )
C         .
          V0( IND ) = V0( IND ) + DERV( 1 )
          V1( IND ) = V1( IND ) + DERV( 2 )
          V2( IND ) = V2( IND ) + DERV( 3 )
          V3( IND ) = V3( IND ) + DERV( 4 )
          V4( IND ) = V4( IND ) + DERV( 5 )
C         .
        ENDDO
 60   CONTINUE
      ENDDO
C     .
      RETURN
      END
C*********************************************************************
