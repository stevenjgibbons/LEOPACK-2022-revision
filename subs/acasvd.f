C*********************************************************************
C subroutine Advanced Complete Adapted Solution Vector Derivative ****
C            -        -        -       -        -      -          ****
C Steve Gibbons Wed Mar  8 12:46:47 GMT 2000                         C
C____________________________________________________________________C
C                                                                    C
C Loops around the entire convection vector V and calculates         C
C derivatives 0 and up to 4 of every harmonic from grid node         C
C IR = ILN to grid node IRN.                                         C
C                                                                    C
C The derivatives are added to the arrays V0, V1, V2, V3 and V4.     C
C                                                                    C
C V is indexed by INIA, MHLI, MHMI, MHTI and MHPI.                   C
C V0, V1, V2, V3 and V4 are indexed by INOA, MHLO, MHMO, MHTO and    C
C MHPO. If ICLS = 1 then the arrays V0, V1, V2, V3 and V4 are all    C
C zeroed a priori. FAC multiplied by the derivatives are then added  C
C to the output arrays.                                              C
C                                                                    C
C If a harmonic in vector V is not found in the output vector,       C
C the program aborts.                                                C
C____________________________________________________________________C
C                                                                    C
C Input variables :-                                                 C
C ===============                                                    C
C  Integer                                                           C
C  -------                                                           C
C                                                                    C
C     ICLS      : Set to 1 to zero arrays V0, V1, V2, V3 and V4.     C
C     ILN       : Left-most node to be acted upon.                   C
C     IRN       : Right-most node to be acted upon.                  C
C     NBN       : Number of bounding nodes. See above.               C
C     IHD       : Highest derivative requested.                      C
C     NCFM      : Leading dimension of SVFDC. At least (2*NBN+1)     C
C     NR        : Number of radial grid nodes in each function.      C
C     NDRVS     : Number of highest derivative for which             C
C                  coefficients are stored by the array SVFDC.       C
C     NDRVM     : Limit on NDRVS. Array bound for SVFDC.             C
C     NDCS      : Maximum distinct finite difference schemes         C
C                  stored by the array SVFDC.                        C
C                                                                    C
C     INIA      : Integer array dimension ( * ). (Input vector)      C
C                  INIA( 1 ) = IFORMI - flag for vector format.      C
C                                        See INDFUN                  C
C                  INIA( 2 ) = NRI. Must be consistent with NR.      C
C                  INIA( 3 ) = NHI = total number of radial func.s   C
C                                                                    C
C     MHTI      : Harmonic type (input vector) (1,2,3,4,5). Dim(*)   C
C     MHLI      : Degree, l. (input vector). Dim(*)                  C
C     MHMI      : Order, m/-m. (input vector). Dim(*)                C
C                                                                    C
C     MHPI      : Pointer array for harmonics. If HMPI( ih ) = is    C
C                  then 'is' is the finite difference scheme used    C
C                  to take derivatives of that harm. radial func.    C
C                  If MHPI is negative, the harmonic is avoided and  C
C                  ACASVD moves on to the next harmonic.             C
C                                                                    C
C     INOA      : Integer array dimension ( * ). (Output vector)     C
C                  INOA( 1 ) = IFORMO - flag for vector format.      C
C                                        See INDFUN                  C
C                  INOA( 2 ) = NRO. Must be consistent with NR.      C
C                  INOA( 3 ) = NHO = total number of radial func.s   C
C                                                                    C
C     MHTO      : Harmonic type (output vector) (1,2,3,4,5). Dim(*)  C
C     MHLO      : Degree, l. (input vector). Dim(*)                  C
C     MHMO      : Order, m/-m. (input vector). Dim(*)                C
C                                                                    C
C  Double Precision                                                  C
C  ----------------                                                  C
C                                                                    C
C     V         : Solution vector. Dim ( * ) but length atleast      C
C                  NR*NHI.                                           C
C                                                                    C
C     SVFDC     : Finite difference coefficient matrix.              C
C                  Dimension ( NCFM, NR, NDRVM+1, NDCS ).            C
C                   Array is generated by the routine svfdcf         C
C                 See documentation for SVFDCF for details.          C
C                                                                    C
C     V0        : Zero^th derivatives. Dim ( * ) but atleast NR*NHO  C
C     V1        : First   derivatives. Dim ( * ) but atleast NR*NHO  C
C     V2        : Second  derivatives. Dim ( * ) but atleast NR*NHO  C
C     V3        : Third   derivatives. Dim ( * ) but atleast NR*NHO  C
C     V4        : Fourth  derivatives. Dim ( * ) but atleast NR*NHO  C
C                                                                    C
C IMPORTANT: ALL of V0, V1, V2, V3 and V4 are filled even if         C
C these derivatives are not requested!                               C
C                                                                    C
C     FAC       : Multiple of derivatives to be added to output      C
C                  vectors.                                          C
C____________________________________________________________________C
C
C*********************************************************************
      SUBROUTINE ACASVD ( ICLS, V, ILN, IRN, NBN, IHD, NCFM, NR,
     1                    NDRVS, NDRVM, NDCS, INIA, MHTI, MHLI, MHMI,
     2                    MHPI, INOA, MHTO, MHLO, MHMO, SVFDC, V0, V1,
     3                    V2, V3, V4, FAC )
      IMPLICIT NONE
C____________________________________________________________________C
C Variable declarations - Parameters ................................C
      INTEGER ICLS, NBN, IHD, NCFM, NR, NDRVS, NDRVM, INIA( * ),
     1        MHTI( * ), MHLI( * ), MHMI( * ), MHPI( * ), INOA( * ),
     2        MHTO( * ), MHLO( * ), MHMO( * ), NDCS, ILN, IRN
      DOUBLE PRECISION V( * ), V0( * ), V1( * ), V2( * ), V3( * ), 
     1                 V4( * ), SVFDC( NCFM, NR, NDRVM+1, NDCS ), FAC
C____________________________________________________________________C
C Variable declarations - Working variables .........................C
C
      DOUBLE PRECISION DERV( 5 ), ZERO
      INTEGER IR, IHI, NHI, IND, INDFUN, IS, N2, IOP, IL, IM, IT,
     1        IHO, NHO
      PARAMETER ( IOP = 0, ZERO = 0.0d0 )
      LOGICAL OK
C____________________________________________________________________C
C START OF PROGRAM **************************************************C
C____________________________________________________________________C
C
      N2 = 5
      CALL VECOP( DERV, ZERO, N2, IOP )
C
      IF ( IHD.LT.0 .OR. IHD.GT.4 ) THEN
        PRINT *,' Subroutine ACASVD. IHD = ', IHD
        PRINT *,' Program aborted.'
        STOP
      ENDIF
C
      IF ( INIA( 2 ).NE.NR .OR. INOA( 2 ).NE.NR ) THEN
        PRINT *,' Subroutine ACASVD.'
        PRINT *,' NR = ', NR
        PRINT *,' INIA(2) = ',INIA(2)
        PRINT *,' INOA(2) = ',INOA(2)
        PRINT *,' Program aborted.'
        STOP
      ENDIF
C
      NHO = INOA( 3 )
      NHI = INIA( 3 )
      N2 = NHO*INOA( 2 )
      IF ( ICLS.EQ.1 ) THEN
        CALL VECOP( V0, ZERO, N2, IOP )
        CALL VECOP( V1, ZERO, N2, IOP )
        CALL VECOP( V2, ZERO, N2, IOP )
        CALL VECOP( V3, ZERO, N2, IOP )
        CALL VECOP( V4, ZERO, N2, IOP )
      ENDIF
C
      DO IHI = 1, NHI
        IS  = MHPI( IHI )
        IF ( IS.LT.1 ) GOTO 60
        OK = .FALSE.
        IT = MHTI( IHI )
        IL = MHLI( IHI )
        IM = MHMI( IHI )
        DO IHO = 1, NHO
C         .
          IF (   MHTO( IHO ).NE.IT .OR. MHLO( IHO ).NE.IL .OR.
     1           MHMO( IHO ).NE.IM ) GOTO 50
          OK = .TRUE.
          DO IR = ILN, IRN
C           .
            CALL ASVDR ( V, IR, IS, IHI, NBN, IHD, NCFM, NR, NDRVS,
     1                   NDRVM, DERV, INIA, SVFDC, NDCS )
C           .
            IND = INDFUN( IR, IHO, INOA )
            V0( IND ) = V0( IND ) + FAC*DERV( 1 )
            V1( IND ) = V1( IND ) + FAC*DERV( 2 )
            V2( IND ) = V2( IND ) + FAC*DERV( 3 )
            V3( IND ) = V3( IND ) + FAC*DERV( 4 )
            V4( IND ) = V4( IND ) + FAC*DERV( 5 )
C           .
          ENDDO
 50       CONTINUE
C         .
        ENDDO
        IF ( .NOT. OK ) THEN
          PRINT *,' Subroutine ACASVD.'
          PRINT *,' No harmonic in output vector'
          PRINT *,' found with type = ', IT
          PRINT *,' l = ',IL,' m = ', IM
          PRINT *,' Program aborted.'
          STOP
        ENDIF
 60   CONTINUE
      ENDDO
C     .
      RETURN
      END
C*********************************************************************
