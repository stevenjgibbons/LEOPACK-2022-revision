C*********************************************************************
C subroutine Coriolis Force curl Vector Interaction Coef. Calculate **
C            -        -          -      -           -     -         **
C Steve Gibbons Mon Jan 17 17:01:14 GMT 2000                         C
C____________________________________________________________________C
C                                                                    C
C  Calculates the coefficients for the interaction                   C
C                                                                    C
C   curl ( k x v_{alpha} )                                           C
C                                                                    C
C  where k = ( cos theta, - sin theta , 0 ).                         C
C                                                                    C
C The output takes the form of 4 arrays; IHNALP, IHNGAM, CVI and     C
C TVHI. IHNALP( in ) [in = interaction number] refers to the number  C
C of the alpha harmonic - and likewise IHNGAM for the gamma harm.s   C
C                                                                    C
C If a given reaction forming the gamma harmonic is non-zero         C
C then the coefficient K^{ag}_{AG} is stored in                      C
C CVI( in ) [Coefficient of Vector Interaction].                     C
C  IHNGAM( in ) is the number defined by INDSHC for this gamma.      C
C                                                                    C
C  The corresponding element of TVHI [Type of Vector Harmonic        C
C Interaction] is given as 'KCAG' where both of A and G are          C
C replaced by Q, S or T depending upon the nature of the vector      C
C harmonics.                                                         C
C                                                                    C
C  The explicit formulae for the interactions 'KCAG' is given by     C
C the equations (B.42) through to (B.44) on page 188 of my thesis.   C
C                                                                    C
C NVI is the Number of Vector Interactions and, as an input          C
C parameter, indicates how many elements are already stored in the   C
C four arrays. (NVI=0 means that no such routine has yet been        C
C called; NVI .gt. 0 means that the coefficients found by this       C
C routine will begin at NVI + 1).                                    C
C                                                                    C
C MAXNVI is the maximum permitted number of vector interactions.     C
C                                                                    C
C____________________________________________________________________C
C                                                                    C
C Input variables :-                                                 C
C ===============                                                    C
C  Integer                                                           C
C  -------                                                           C
C                                                                    C
C     NVI      : Number of vector interactions.                      C
C     MAXNVI   : Maximum number of vector interactions.              C
C                This is an upper limit for NVI and the dimension    C
C                of the arrays IHNALP, IHNGAM, CVI and TVHI.         C
C                                                                    C
C     IHNALP   : Number of alpha harmonics. Dim ( MAXNVI )           C
C     IHNGAM   : Number of gamma harmonics. Dim ( MAXNVI )           C
C                                                                    C
C   [Key for MTI, MTO :- 1 = poloidal velocity, 2 = toroidal         C
C  velocity: 3, 4 and 5 are temp, pol mag. field and tor mag. f.     C
C  respectively but are not relevant to this routine.]               C
C                                                                    C
C     NH       : Number of alpha (input) and gamma (output)          C
C                 harmonics.                                         C
C                                                                    C
C     MTI      : Array length ( * ) - atleast length NH              C
C                  See above for key. (corresponds to alpha vec.)    C
C     MLI      : Array length ( * ) - atleast length NH              C
C                  Sph. harm. degree, l.                             C
C     MMI      : Array length ( * ) - atleast length NH              C
C                  Sph. harm. order, m, for cos m phi dep.           C
C                 -Sph. harm. order, m, for sin m phi dep.           C
C                                                                    C
C     MTO      : Array length ( * ) - atleast length NH              C
C                  See above for key. (corresponds to gamma vec.)    C
C     MLO      : Array length ( * ) - atleast length NH              C
C                  Sph. harm. degree, l.                             C
C     MMO      : Array length ( * ) - atleast length NH              C
C                  Sph. harm. order, m, for cos m phi dep.           C
C                 -Sph. harm. order, m, for sin m phi dep.           C
C                                                                    C
C     LH        : Maximum spherical harmonic degree, l.              C
C     NTHP      : The number of theta points.                        C
C     NPHP      : The number of phi points.                          C
C     MMAX      : Maximum sph. harmonic order, m.                    C
C                                                                    C
C  Character                                                         C
C  ---------                                                         C
C                                                                    C
C     TVHI      : *(4) Type of vector interaction. Dim. (MAXNVI).    C
C                 = 'KCSS', 'KCST' etc. according to the corresp.    C
C                 vector interaction.                                C
C                                                                    C
C  Double Precision                                                  C
C  ----------------                                                  C
C                                                                    C
C     CVI       : Coefficient of vector interaction. Dim ( MAXNVI )  C
C                                                                    C
C     GAUX      : Cosines of the NTHP evaluated by the routine       C
C                  gauwts. Dimension ( NTHP ).                       C
C     GAUW      : Corresponding weights. As above.                   C
C     PA        : Schmidt Normalised Legendre Functions              C
C                  Dimension ( ( LH + 1 )*( LH + 2 )/2, NTHP )       C
C                  P_l^m(X_j) = PA( L*(L+1)/2 + M + 1 , j )          C
C     DPA       : Derivatives of the above.                          C
C                                                                    C
C     FTF1      : Work array - dim. (2*NPHP)                         C
C     FTF2      : Work array - dim. (2*NPHP)                         C
C     FTF3      : Work array - dim. (2*NPHP)                         C
C     VF        : Work array - dim. ( NPHP, NTHP, 3)                 C
C     QST       : Work array - dim. ( LH*(LH+2), 3)                  C
C                                                                    C
C____________________________________________________________________C
C
C*********************************************************************
      SUBROUTINE CFVICC( NVI, MAXNVI, IHNALP, IHNGAM, NH, MTI, MLI,
     1          MMI, MTO, MLO, MMO, LH, NTHP, NPHP, MMAX, TVHI, CVI,
     2          GAUX, GAUW, PA, DPA, FTF1, FTF2, FTF3, VF, QST )
      IMPLICIT NONE
C____________________________________________________________________C
C Variable declarations - Parameters ................................C
      INTEGER NVI, MAXNVI, IHNALP( * ), IHNGAM( * ), NH, MTI( * ),
     1        MLI( * ), MMI( * ), MTO( * ), MLO( * ), MMO( * ),
     2        LH, NTHP, NPHP, MMAX
      CHARACTER *(4) TVHI( MAXNVI )
      DOUBLE PRECISION 
     1                 CVI( MAXNVI ), QST( LH*(LH+2), 3),
     2                 VF( NPHP, NTHP, 3), FTF1( 2*NPHP ),
     3                 FTF2( 2*NPHP ), FTF3( 2*NPHP )
      DOUBLE PRECISION
     1                 GAUX( NTHP ), GAUW( NTHP ),
     2                 PA ( ( LH + 1 )*( LH + 2 )/2 , NTHP),
     3                 DPA ( ( LH + 1 )*( LH + 2 )/2 , NTHP)
C____________________________________________________________________C
C Variable declarations - Working variables .........................C
C
      INTEGER IHI, IHO, LI, MI, ICSI, IQSTI, INDI,
     2        LO, MO, ICSO, INDO, INDSHC
      DOUBLE PRECISION LOW, QAB, SAB, TAB
      PARAMETER ( LOW = 1.0d-9 )
      LOGICAL OK
C____________________________________________________________________C
C START OF PROGRAM **************************************************C
C____________________________________________________________________C
C     .
      OK = .TRUE.
C     .
C     . Check that NVI is valid.
C     .
      IF ( NVI.LT.0 .OR. NVI.GT.MAXNVI ) THEN
        PRINT *,' Subroutine CFVICC.'
        PRINT *,' NVI = ', NVI,' MAXNVI = ', MAXNVI
        PRINT *,' Program aborted.'
        STOP
      ENDIF
C     .
C     . Loop around in harmonics
C     .
      DO IHI = 1, NH
        IF ( MTI( IHI ).NE.1 .AND. MTI( IHI ).NE.2 ) GOTO 50
C
C Consider separately the cases of input harmonic
C being poloidal and toroidal.
C
        IF ( MTI( IHI ).EQ.1 ) THEN
C         .
C         . Input harmonic is poloidal ...
C         .
          LI = MLI( IHI )
          IF ( MMI( IHI ).LT.0 ) THEN
            MI   = -MMI( IHI )
            ICSI = 2
          ELSE
            MI   = MMI( IHI )
            ICSI = 1
          ENDIF
          INDI   = INDSHC( LI, MI, ICSI )
C         .
C         . Call CORCOF with IQSTI = 1.
C         . This will give us K_{qs} and K_{qt} 
C         . (c.f. eq B.42 in thesis) - K_{qq} should be zero!
C         .
          IQSTI = 1
          CALL CORCOF( IQSTI, LI, MI, ICSI, QST, LH, NPHP, NTHP,
     1           FTF1, FTF2, FTF3, VF, GAUX, GAUW, PA, DPA, MMAX )
C         .
C         . Now loop around the output harmonics
C         .
          DO IHO = 1, NH
C           .
            LO = MLO( IHO )
            IF ( MMO( IHO ).LT.0 ) THEN
              MO   = -MMO( IHO )
              ICSO = 2
            ELSE
              MO   = MMO( IHO )
              ICSO = 1
            ENDIF
            INDO = INDSHC( LO, MO, ICSO )
            IF( INDO.GT.0 ) SAB   = QST( INDO, 2 )
            IF( INDO.GT.0 ) TAB   = QST( INDO, 3 )
C           .
C           . Case of poloidal harmonic
C           .
            IF ( MTO( IHO ).EQ.1 .AND. DABS( TAB ).GT.LOW ) THEN
C              .
C              . Need to do K_{qt} interaction
C              .
               CALL CNTRIC( NVI, MAXNVI, OK )
               IF ( OK ) THEN
                 IHNALP( NVI ) = INDI
                 IHNGAM( NVI ) = INDO
                 CVI( NVI )    = TAB
                 TVHI( NVI )   = 'KCQT'
               ENDIF
C              .
            ENDIF
C           .
C           . Case of toroidal harmonic
C           .
            IF ( MTO( IHO ).EQ.2 .AND. DABS( SAB ).GT.LOW ) THEN
C              .
C              . Need to do K_{qs} interaction
C              .
               CALL CNTRIC( NVI, MAXNVI, OK )
               IF ( OK ) THEN
                 IHNALP( NVI ) = INDI
                 IHNGAM( NVI ) = INDO
                 CVI( NVI )    = SAB
                 TVHI( NVI )   = 'KCQS'
               ENDIF
C              .
            ENDIF
C           .
          ENDDO
C         .
C         . Call CORCOF with IQSTI = 2.
C         . This will give us K_{sq}, K_{ss} and K_{st}
C         . (c.f. eq B.43 in thesis)
C         .
          IQSTI = 2
          CALL CORCOF( IQSTI, LI, MI, ICSI, QST, LH, NPHP, NTHP,
     1           FTF1, FTF2, FTF3, VF, GAUX, GAUW, PA, DPA, MMAX )
C         .
C         . Now loop around the output harmonics
C         .
          DO IHO = 1, NH
C           .
            LO = MLO( IHO )
            IF ( MMO( IHO ).LT.0 ) THEN
              MO   = -MMO( IHO )
              ICSO = 2
            ELSE
              MO   = MMO( IHO )
              ICSO = 1
            ENDIF
            INDO = INDSHC( LO, MO, ICSO )
            IF( INDO.GT.0 ) QAB   = QST( INDO, 1 )
            IF( INDO.GT.0 ) SAB   = QST( INDO, 2 )
            IF( INDO.GT.0 ) TAB   = QST( INDO, 3 )
C           .
C           . Case of poloidal harmonic
C           .
            IF ( MTO( IHO ).EQ.1 .AND. DABS( TAB ).GT.LOW ) THEN
C              .
C              . Need to do K_{st} interaction
C              .
               CALL CNTRIC( NVI, MAXNVI, OK )
               IF ( OK ) THEN
                 IHNALP( NVI ) = INDI
                 IHNGAM( NVI ) = INDO
                 CVI( NVI )    = TAB
                 TVHI( NVI )   = 'KCST'
               ENDIF
C              .
            ENDIF
C           .
C           . Case of toroidal harmonic
C           .
            IF ( MTO( IHO ).EQ.2 ) THEN
C              .
C              . Need to do K_{sq} interaction
C              .
               IF ( ABS( QAB ).GT.LOW ) THEN
C                .
                 CALL CNTRIC( NVI, MAXNVI, OK )
                 IF ( OK ) THEN
                   IHNALP( NVI ) = INDI
                   IHNGAM( NVI ) = INDO
                   CVI( NVI )    = QAB
                   TVHI( NVI )   = 'KCSQ'
                 ENDIF
C                .
               ENDIF
C              .
C              . Need to do K_{ss} interaction
C              .
               IF ( ABS( SAB ).GT.LOW ) THEN
C                .
                 CALL CNTRIC( NVI, MAXNVI, OK )
                 IF ( OK ) THEN
                   IHNALP( NVI ) = INDI
                   IHNGAM( NVI ) = INDO
                   CVI( NVI )    = SAB
                   TVHI( NVI )   = 'KCSS'
                 ENDIF
C                .
               ENDIF
C              .
            ENDIF
C           .
          ENDDO
C         .
        ENDIF
C       endif case of ihi poloidal
C    
        IF ( MTI( IHI ).EQ.2 ) THEN
C         .
C         . Input harmonic is toroidal ...
C         .
          LI = MLI( IHI )
          IF ( MMI( IHI ).LT.0 ) THEN
            MI   = -MMI( IHI )
            ICSI = 2
          ELSE
            MI   = MMI( IHI )
            ICSI = 1
          ENDIF
          INDI   = INDSHC( LI, MI, ICSI )
C         .
C         . Call CORCOF with IQSTI = 3.
C         . This will give us K_{tq}, K_{ts} and K_{tt}
C         . (c.f. eq B.44 in thesis)
C         .
          IQSTI = 3
          CALL CORCOF( IQSTI, LI, MI, ICSI, QST, LH, NPHP, NTHP,
     1           FTF1, FTF2, FTF3, VF, GAUX, GAUW, PA, DPA, MMAX )
C         .
C         . Now loop around the output harmonics
C         .
          DO IHO = 1, NH
C           .
            LO = MLO( IHO )
            IF ( MMO( IHO ).LT.0 ) THEN
              MO   = -MMO( IHO )
              ICSO = 2
            ELSE
              MO   = MMO( IHO )
              ICSO = 1
            ENDIF
            INDO = INDSHC( LO, MO, ICSO )
            IF( INDO.GT.0 ) QAB   = QST( INDO, 1 )
            IF( INDO.GT.0 ) SAB   = QST( INDO, 2 )
            IF( INDO.GT.0 ) TAB   = QST( INDO, 3 )
C           .
C           . Case of poloidal harmonic
C           .
            IF ( MTO( IHO ).EQ.1 .AND. DABS( TAB ).GT.LOW ) THEN
C              .
C              . Need to do K_{tt} interaction
C              .
               CALL CNTRIC( NVI, MAXNVI, OK )
               IF ( OK ) THEN
                 IHNALP( NVI ) = INDI
                 IHNGAM( NVI ) = INDO
                 CVI( NVI )    = TAB
                 TVHI( NVI )   = 'KCTT'
               ENDIF
C              .
            ENDIF
C           .
C           . Case of toroidal harmonic
C           .
            IF ( MTO( IHO ).EQ.2 ) THEN
C              .
C              . Need to do K_{tq} interaction
C              .
               IF ( ABS( QAB ).GT.LOW ) THEN
C                .
                 CALL CNTRIC( NVI, MAXNVI, OK )
                 IF ( OK ) THEN
                   IHNALP( NVI ) = INDI
                   IHNGAM( NVI ) = INDO
                   CVI( NVI )    = QAB
                   TVHI( NVI )   = 'KCTQ'
                 ENDIF
C                .
               ENDIF
C              .
C              . Need to do K_{ts} interaction
C              .
               IF ( ABS( SAB ).GT.LOW ) THEN
C                .
                 CALL CNTRIC( NVI, MAXNVI, OK )
                 IF ( OK ) THEN
                   IHNALP( NVI ) = INDI
                   IHNGAM( NVI ) = INDO
                   CVI( NVI )    = SAB
                   TVHI( NVI )   = 'KCTS'
                 ENDIF
C                .
               ENDIF
C              .
            ENDIF
C           .
          ENDDO
C         .
        ENDIF
C       endif case of ihi toroidal
C    
 50   CONTINUE
      ENDDO
C     .
      IF ( OK ) RETURN
      PRINT *,' Subroutine CFVICC.'
      PRINT *,' Number of required interactions = ', NVI
      PRINT *,' Maximum interactions allowed    = ', MAXNVI
      PRINT *,' Program aborted.'
      STOP
      END
C*********************************************************************
