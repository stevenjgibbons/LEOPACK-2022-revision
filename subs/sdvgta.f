C*********************************************************************
C subroutine Stored Derivative Velocity dot Gradient of Theta Add ****
C            -      -          -            -           -     -   ****
C Steve Gibbons Mon Feb 14 18:18:21 GMT 2000                         C
C____________________________________________________________________C
C                                                                    C
C If the vector VI contains radial functions to NHI spherical harm.s C
C with type defined by MHTI, MHLI and MHMI, then SDVGTA will add FAC C
C multiplied by the V.Grad(theta) term to the appropriate locations  C
C in another solution vector VO whose NHO radial functions are       C
C defined by MHTO, MHLO and MHMO.                                    C
C                                                                    C
C Important! The velocity, v, is supplied in the array RQST and      C
C must be obtained in advance by calling ASRQST!!                    C
C SDVGTA does NOT extract the velocity from VI. This is done         C
C entirely to reduce the programming load and ease of checking.      C
C                                                                    C
C Grad theta IS calculated by SDVGTA from the harmonics in VI with   C
C MHTI( ih ) = 3.                                                    C
C                                                                    C
C                                                                    C
C____________________________________________________________________C
C                                                                    C
C Input variables :-                                                 C
C ===============                                                    C
C  Integer                                                           C
C  -------                                                           C
C                                                                    C
C     NR        : Number of radial grid nodes.                       C
C     NDCS      : Maximum distinct finite difference schemes         C
C                  stored by the array SVFDC.                        C
C     LH        : Maximum spherical harmonic degree, l.              C
C     MMAX      : Maximum spherical harmonic order, m.               C
C     MHTI      : Array length ( * ) - atleast length NHI            C
C                  See below for key. (corresponds to input vec.)    C
C                                                                    C
C     MHTI( i ) = 1 for a poloidal velocity harmonic, i.             C
C     MHTI( i ) = 2 for a toroidal velocity harmonic, i.             C
C     MHTI( i ) = 3 for a temperature harmonic, i.                   C
C     MHTI( i ) = 4 for a poloidal magnetic field harmonic, i.       C
C     MHTI( i ) = 5 for a toroidal magnetic field harmonic, i.       C
C                                                                    C
C     MHLI      : Array length ( * ) - atleast length NHI            C
C                  Sph. harm. degree, l.                             C
C     MHMI      : Array length ( * ) - atleast length NHI            C
C                  Sph. harm. order, m, for cos m phi dep.           C
C                 -Sph. harm. order, m, for sin m phi dep.           C
C                                                                    C
C     INARRI    : Int. parameter array corresponding to VI.          C
C                 Dim ( * ) but the following must be true.          C
C                                                                    C
C                 INARRI( 1 ) = IFORMF                               C
C                 INARRI( 2 ) = NRI     See INDFUN for details       C
C                 INARRI( 3 ) = NHI                                  C
C                                                                    C
C     Note that INARRO, MHTO, MHLO and MHMO correspond to            C
C    exactly the above variables but for the output and not input    C
C    vectors. NRI must equal NRO and must both equal NR.             C
C                                                                    C
C     NTHPTS    : Number of theta points.                            C
C     NPHPTS    : Number of phi points.                              C
C                                                                    C
C     ILNR      : First radial node to act on.                       C
C     IRNR      : Last radial node to act on.                        C
C                                                                    C
C  Double Precision                                                  C
C  ----------------                                                  C
C                                                                    C
C     V0        : Solution vector. Dim ( * ) ( input )               C
C                 Length must be atleast NRI*NHI                     C
C     V1        : Solution vector. Dim ( * ) ( input )               C
C                 First deriv.s of above. Calc. by CASVDR.           C
C     VO        : Solution vector. Dim ( * )  ( output )             C
C                 Length must be atleast NRO*NHO                     C
C     FAC       : Multiplication factor.                             C
C     RQST      : Input array containing scaloidal/spheroidal/tor.   C
C                  decomposition of velocity vector.                 C
C                  Has dimensions (  LH*(LH+2) ,3, NR ).             C
C              RQST (l*l+2m,1,I) = q_l^ms(r_i)                       C
C   RQST (l*l+2m-1+delta_m0,1,I) = q_l^mc(r_i)                       C
C              RQST (l*l+2m,2,I) = s_l^ms(r_i)                       C
C   RQST (l*l+2m-1+delta_m0,2,I) = s_l^mc(r_i)                       C
C              RQST (l*l+2m,3,I) = t_l^ms(r_i)                       C
C   RQST (l*l+2m-1+delta_m0,3,I) = t_l^mc(r_i)                       C
C                                                                    C
C     VF1       : Work arr. dim. ( NPHPTS, NTHPTS, 3 )               C
C     VF2       : Work arr. dim. ( NPHPTS, NTHPTS, 3 )               C
C     SF        : Work arr. dim. ( NPHPTS, NTHPTS )                  C
C     FTF1      : Work arr. dim. ( 2*NPHPTS )                        C
C     FTF2      : Work arr. dim. ( 2*NPHPTS )                        C
C     FTF3      : Work arr. dim. ( 2*NPHPTS )                        C
C                                                                    C
C     SHC       : Work array. Dim ( LH*(LH+2) ).                     C
C     DSHC      : Work array. Dim ( LH*(LH+2) ).                     C
C                                                                    C
C     XARR      : Values of radial nodes. Dim ( NR )                 C
C                                                                    C
C     ZCOEFA    : Coefficients of the monopole term in the R compon. C
C                 This is not necessarily zero in general (but ought C
C                 to be for this procedure!!! )                      C
C                 Dimension ( NR )                                   C
C                                                                    C
C     PA        : Schmidt Normalised Legendre Functions              C
C                  Dim { ( LH + 1 )*( LH + 2 )/2 , NTHPTS }          C
C                  P_l^m(X_j) = PA( L*(L+1)/2 + M + 1 , j )          C
C     DPA       : Derivatives of the above.                          C
C                                                                    C
C____________________________________________________________________C
C
C*********************************************************************
      SUBROUTINE SDVGTA ( NR, LH, MMAX, MHTI, MHLI, MHMI, V0, V1,
     1             INARRI, MHTO, MHLO, MHMO, VO, INARRO, RQST,
     2             VF1, VF2, SF, FTF1, FTF2, FTF3, SHC, DSHC, GAUX,
     3             GAUW, PA, DPA, NTHPTS, NPHPTS, FAC, ILNR, IRNR,
     4             XARR, ZCOEFA)
      IMPLICIT NONE
C____________________________________________________________________C
C Variable declarations - Parameters ................................C
C
      INTEGER NR, LH, MMAX, NTHPTS, NPHPTS, ILNR, IRNR, 
     1        MHTI( * ), MHLI( * ), MHMI( * ), INARRI( * ),
     2        MHTO( * ), MHLO( * ), MHMO( * ), INARRO( * )
      DOUBLE PRECISION V0( * ), V1( * ), VO( * ),
     1                 RQST( LH*(LH+2), 3, NR ), GAUX( NTHPTS ),
     2                 GAUW( NTHPTS ), ZCOEFA( NR ), FAC
      DOUBLE PRECISION XARR( NR ), FTF1( 2*NPHPTS ),
     1                 FTF2( 2*NPHPTS ), FTF3( 2*NPHPTS ),
     2                 PA ( ( LH + 1 )*( LH + 2 )/2 , NTHPTS),
     3                 DPA ( ( LH + 1 )*( LH + 2 )/2 , NTHPTS)
      DOUBLE PRECISION SHC( LH*(LH+2) ), DSHC( LH*(LH+2) ),
     2                 SF( NPHPTS, NTHPTS )
      DOUBLE PRECISION VF1( NPHPTS, NTHPTS, 3 ),
     1                 VF2( NPHPTS, NTHPTS, 3 )
C____________________________________________________________________C
C Variable declarations - Working variables .........................C
      INTEGER NRI, NRO, NHI, NHO, IR, IOP, NH, IHI, IHARM,
     1        IND, L, M, ICS, IHO, INDSHC, INDFUN
      DOUBLE PRECISION ZERO, D0F, D1F, DZCOEF, ZCOEF, RAD
      PARAMETER ( ZERO = 0.0d0 )
C____________________________________________________________________C
C START OF PROGRAM **************************************************C
C____________________________________________________________________C
C
      NRI    = INARRI( 2 )
      NHI    = INARRI( 3 )
C
      NRO    = INARRO( 2 )
      NHO    = INARRO( 3 )
C
      NH     = LH*(LH + 2)
C
      IF ( NRI.NE.NR .OR. NRO.NE.NR ) THEN
         PRINT *,' Subroutine SDVGTA.'
         PRINT *,' NR   = ', NR
         PRINT *,' NRI  = ', NRI
         PRINT *,' NRO  = ', NRO
         PRINT *,' Program aborted.'
         STOP
      ENDIF
C     .
C     . Loop around radial grid nodes from
C     . IR = ILNR to IR = IRNR.
C     .
      DO IR = ILNR, IRNR
C       .
        RAD = XARR( IR )
C       .
C       . Fill VF1 with the velocity evaluated at node IR
C       .
        CALL RQSTVF ( RQST, VF1, GAUX, PA, DPA, FTF1, FTF2, FTF3,
     1                LH, NTHPTS, NPHPTS, MMAX, ZCOEFA, NR, IR )
C       .
C       . First let's zero the arrays SHC and DSHC for node IR
C       .
        IOP = 0
        CALL VECOP( SHC, ZERO, NH, IOP )
        CALL VECOP( DSHC, ZERO, NH, IOP )
        DZCOEF = ZERO
C       .
C       . Let's loop around harmonics, ihi, in
C       . VI and for each one with MHTI( ihi ) = 3
C       . we put an entry into SHC and DSHC
C       .
        DO IHI = 1, NHI
          IF (    MHTI( IHI ).NE.3     )    GOTO 50
          L  = MHLI( IHI )
          IF ( L.GT.LH ) THEN
            PRINT *,' Subroutine SDVGTA.'
            PRINT *,' Harmonic ',IHI,' has L = ',L
            PRINT *,' LH = ', LH
            PRINT *,' Program aborted.'
            STOP
          ENDIF
          IF (          MHMI( IHI ).GE.0      ) THEN
            M =  MHMI( IHI )
            ICS = 1
          ELSE
            M = -MHMI( IHI )
            ICS = 2
          ENDIF
          IHARM = INDSHC( L, M, ICS )
C         .
C         . Calculate D0F and D1F
C         .
          IND = INDFUN( IR, IHI, INARRI )
          D0F = V0( IND )
          D1F = V1( IND )
C         .
C         . Add to either DZCOEF, or SHC and DSHC
C         .
          IF ( IHARM.EQ.0 ) THEN
            DZCOEF = D1F
          ELSE
            SHC( IHARM )  = D0F
            DSHC( IHARM ) = D1F
          ENDIF
C         .
 50     CONTINUE
        ENDDO
C       .
C       . SHC, DSHC and DZCOEF are now ready to be
C       . supplied to GRINVT - this will fill VF2
C       . with Grad( theta )
C       .
        CALL GRINVT( SHC, DSHC, GAUX, RAD, PA, DPA, FTF1, FTF2,
     1               FTF3, VF2, LH, NTHPTS, NPHPTS, MMAX, DZCOEF )
C       .
C       . VF1 now contains V and VF2 now contains Grad( theta )
C       . so now take scalar product.
C       .
        CALL VFDP( VF1, VF2, SF, NPHPTS, NTHPTS )
C       .
C       . SF now contains the scalar function in
C       . space - so transform back to shc coeffs
C       .
        CALL FORSST( SHC, SF, GAUW, PA, FTF1, LH, MMAX,
     1               NTHPTS, NPHPTS, ZCOEF )
C       .
C       . Now loop around 'output' harmonics
C       . and place in appropriate locations
C       .
        DO IHO = 1, NHO
          IF (    MHTO( IHO ).NE.3     )    GOTO 51
          L = MHLO( IHO )
          IF (          MHMO( IHO ).GE.0      ) THEN
            M =  MHMO( IHO )
            ICS = 1
          ELSE
            M = -MHMO( IHO )
            ICS = 2
          ENDIF
          IHARM = INDSHC( L, M, ICS )
C
C Next line means harmonic iho is beyond the scale
C resolved by the current transform
C
          IF ( IHARM.GT.NH ) GOTO 51
          IND = INDFUN( IR, IHO, INARRO )
C         .
          IF ( IHARM.EQ.0 ) THEN
            VO( IND ) = VO( IND ) + FAC*ZCOEF
          ELSE
            VO( IND ) = VO( IND ) + FAC*SHC( IHARM )
          ENDIF
C         .
 51     CONTINUE
        ENDDO
C       .
      ENDDO
C     .
      RETURN
      END
C*********************************************************************

      

