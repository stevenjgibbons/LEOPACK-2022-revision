C*********************************************************************
C subroutine QST Array to Vector function array version **************
C            --- -     -- -                             **************
C Steve Gibbons 25.4.97                                              C
C____________________________________________________________________C
C                                                                    C
C Input variables :-                                                 C
C ===============                                                    C
C  Integer                                                           C
C  -------                                                           C
C     NTHPTS    : Number of theta points.                            C
C     NPHPTS    : Number of phi points.                              C
C     LH        : Level of harmonics.                                C
C     IRAD      : Current radial grid node.                          C
C     NR        : Number of radial grid nodes.                       C
C                                                                    C
C  Double Precision                                                  C
C  ----------------                                                  C
C     GAUX      : Cosines of the NTHPTS evaluated by the routine     C
C                  gauwts. Dimension ( NTHPTS ).                     C
C     QST       : Input array containing scaloidal/spheroidal        C
C                  decomposition of vector c.f. eqn (38).            C
C                  Has dimensions (  LH*(LH+2) , 3, NR ).            C
C               QST (l*l+2m,1) = q_l^ms(r_i)                         C
C    QST (l*l+2m-1+delta_m0,1) = q_l^mc(r_i)                         C
C                                (as in equation (39)                C
C               QST (l*l+2m,2) = s_l^ms(r_i)                         C
C    QST (l*l+2m-1+delta_m0,2) = s_l^mc(r_i)                         C
C                                (as in equation (40),(41)           C
C               QST (l*l+2m,3) = t_l^ms(r_i)                         C
C    QST (l*l+2m-1+delta_m0,3) = t_l^mc(r_i)                         C
C                                (as in equation (40),(41)           C
C     PA        : Schmidt Normalised Legendre Functions              C
C                  Dim { ( LH + 1 )*( LH + 2 )/2 , NTHPTS }          C
C                  P_l^m(X_j) = PA( L*(L+1)/2 + M + 1 , j )          C
C     DPA       : Derivatives of the above.                          C
C____________________________________________________________________C
C Output variables :-                                                C
C ================                                                   C
C  Double Precision                                                  C
C  ----------------                                                  C
C     VF	: Vector Function. An array of dimensions            C
C		   ( NPHPTS, NTHPTS, 3) which contain the R, THETA   C
C 		   and PHI components of a VECTOR at each point      C
C                  ... i.e. VF ( IPHI, ITHETA, 2 ) is the Theta      C
C                  compontent of the vector at (iphi, itheta).       C
C____________________________________________________________________C
C Passed in arrays :-                                                C
C ================                                                   C
C  Double Precision                                                  C
C  ----------------                                                  C
C     FTF1      : Array for fourier transforming.                    C
C     FTF2      : Array for fourier transforming.                    C
C     FTF3      : Array for fourier transforming.                    C
C                  Has dimensions ( 2*NPHPTS )                       C
C                                                                    C
C None of these arrays need any input or output values but must be   C
C in parameter list for the sake of dimensioning.                    C
C____________________________________________________________________C
C
C*********************************************************************
      SUBROUTINE QSTA2V ( QST, VF, GAUX, PA, DPA, FTF1, FTF2, FTF3,
     1                    LH, IRAD, NR, NTHPTS, NPHPTS )
      IMPLICIT NONE
C____________________________________________________________________C
C Variable declarations - Parameters ................................C
      INTEGER LH, NTHPTS, NPHPTS, IRAD, NR
      DOUBLE PRECISION QST(  LH*(LH+2) , 3, NR),
     1                 VF( NPHPTS, NTHPTS, 3),
     2                 GAUX( NTHPTS )
      DOUBLE PRECISION FTF1( 2*NPHPTS ),
     1                 FTF2( 2*NPHPTS ), FTF3( 2*NPHPTS ),
     2                 PA ( ( LH + 1 )*( LH + 2 )/2 , NTHPTS),
     3                 DPA ( ( LH + 1 )*( LH + 2 )/2 , NTHPTS)
C____________________________________________________________________C
C Variable declarations - Working variables .........................C
      LOGICAL POWT
      INTEGER IOP, ISIGN, ITHREE, ITHETA, LENV, IPHI, 
     1        IP, IND1, IND2,  L, M, INDFTC, INDFTS
      DOUBLE PRECISION ZERO,X,SINE,TERM1, TERM2
      PARAMETER ( ZERO = 0.0d0 , ITHREE=3 )
C____________________________________________________________________C
C Functions used :-
c     INTEGER INDSHC
      DOUBLE PRECISION SQRLL1
C____________________________________________________________________C
C START OF PROGRAM **************************************************C
C____________________________________________________________________C
C   
C Check the validity of arguments .....
C First that NPHPTS is a powers of 2.
c     CALL POWTWO ( NPHPTS, POWT )
c     IF ( .NOT. POWT ) THEN
c        PRINT *,' Subroutine QSTA2V.'
c        PRINT *,' NPHPTS is not a power of 2.'
c        PRINT *,' Program stopped.'
c        STOP
c     ENDIF
C ......... need to have LH < NPHPTS/2
      IF ( 2*LH.GT.NPHPTS ) THEN
         PRINT *,' Subroutine QSTA2V.'
         PRINT *,' NPHPTS must be at least 2*LH.'
         PRINT *,' Program stopped.'
         STOP
      ENDIF
C .........
C_____________________________________________________________________
      IOP = 0
      ISIGN = -1
C ............. set VF array to zero
      CALL CUBEOP ( VF, ZERO, NPHPTS, NTHPTS, ITHREE, IOP )

C ............. start to loop around theta points ...................
      DO ITHETA = 1, NTHPTS
C
         X = GAUX( ITHETA )
         SINE = DSIN( ACOS( X ) )
C
C ............. set FTF1, FTF2, FTF3 all to zero
         LENV = 2*NPHPTS
C
         CALL VECOP ( FTF1, ZERO, LENV, IOP )
         CALL VECOP ( FTF2, ZERO, LENV, IOP )
         CALL VECOP ( FTF3, ZERO, LENV, IOP )
C
C ............. start to loop around Harmonics ......................
         
         M = 0
         DO L = 1, LH
            IP = L*(L+1)/2+M+1
            FTF1( 1 ) = FTF1( 1 ) 
     1           + PA( IP , ITHETA )*QST( L*L , 1, IRAD )
            FTF2( 1 ) = FTF2( 1 ) 
     1           + DPA( IP , ITHETA )*QST( L*L, 2, IRAD)/SQRLL1(L)
            FTF3( 1 ) = FTF3( 1 ) 
     1           + DPA( IP , ITHETA )*QST( L*L, 3, IRAD)/SQRLL1(L)
         ENDDO
C ................. end of case M = 0
         DO M = 1, LH
            INDFTC = 2*M + 1
            INDFTS = 2*M + 2
            DO L = M, LH
               IP = L*(L+1)/2+M+1
c              ICS = 1
c              IND1 = INDSHC( L, M, ICS )
               IND1 = L*L + 2*M - 1
c              ICS = 2
c              IND2 = INDSHC( L, M, ICS )
               IND2 = L*L + 2*M
               TERM1 = M*PA( IP , ITHETA )/( SINE*SQRLL1( L ) )
               TERM2 = DPA( IP , ITHETA )/SQRLL1( L )
C
               FTF1( INDFTC ) = FTF1( INDFTC ) +
     1          PA( IP , ITHETA )*QST( IND1, 1, IRAD )
               FTF1( INDFTS ) = FTF1( INDFTS ) +
     1          PA( IP , ITHETA )*QST( IND2, 1, IRAD )
               FTF2( INDFTC ) = FTF2( INDFTC ) +
     1          TERM2*QST( IND1, 2, IRAD) - TERM1*QST( IND2, 3, IRAD)
               FTF2( INDFTS ) = FTF2( INDFTS ) +
     1          TERM2*QST( IND2, 2, IRAD) + TERM1*QST( IND1, 3, IRAD)
               FTF3( INDFTC ) = FTF3( INDFTC ) +
     1          TERM2*QST( IND1, 3, IRAD) + TERM1*QST( IND2, 2, IRAD)
               FTF3( INDFTS ) = FTF3( INDFTS ) +
     1          TERM2*QST( IND2, 3, IRAD) - TERM1*QST( IND1, 2, IRAD)
C
            ENDDO
         ENDDO
C ............. ended looping around Harmonics ......................
C
C ............... now perform Fourier Transforms on FTF1, FTF2, FTF3
C
         CALL FFTNEW ( FTF1, NPHPTS, ISIGN )
         CALL FFTNEW ( FTF2, NPHPTS, ISIGN )
         CALL FFTNEW ( FTF3, NPHPTS, ISIGN )
C
C ...................................................................
         DO IPHI = 1, NPHPTS
            VF ( IPHI , ITHETA, 1 ) = FTF1( 2*IPHI - 1 )
            VF ( IPHI , ITHETA, 2 ) = FTF2( 2*IPHI - 1 )
            VF ( IPHI , ITHETA, 3 ) = FTF3( 2*IPHI - 1 )
         ENDDO
C
      ENDDO
C ............. ended looping around theta points ...................
C
      RETURN
      END
C*********************************************************************

