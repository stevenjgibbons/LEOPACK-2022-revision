C*********************************************************************
C subroutine Single Harmonic VECtor transform (Optimised) ************
C            -      -        ---               -          ************
C Steve Gibbons Mon Oct  8 15:05:28 WEST 2001                        C
C____________________________________________________________________C
C Input variables :-                                                 C
C ===============                                                    C
C  Integer                                                           C
C  -------                                                           C
C                                                                    C
C     L         : Spherical harmonic degree, l.                      C
C     M         : Spherical harmonic order, m.                       C
C     ICS       : =1 for cos m phi dependence and 2 for sin m phi.   C
C     IQST	: Number of function.                                C
C                  IQST = 1 : Q - scaloidal harmonic. Evaluates      C
C                     ( Y , 0 , 0 )                                  C
C                  IQST = 2 : S - spheroidal harm. Evaluates         C
C                     ( 0 , dY/dtheta , 1/sintheta dY/dphi )*fac     C
C                  IQST = 3 : T - toroidal harmonic. Evaluates       C
C                     ( 0 , - 1/sintheta dY/dphi , dY/dtheta )*fac   C
C                                                                    C
C   where fac = 1/sqrt( l*l + l )                                    C
C                                                                    C
C     NTHPTS    : Number of theta points.                            C
C     NPHPTS    : Number of phi points.                              C
C     LH        : Maximum spherical harmonic degree, l.              C
C     M0        : Minimum non-zero wavenumber.                       C
C                                                                    C
C  Double Precision                                                  C
C  ----------------                                                  C
C                                                                    C
C     GAUX      : Cosines of the NTHPTS evaluated by the routine     C
C                  gauwts. Dimension ( NTHPTS ).                     C
C     PA         : Schmidt Normalised Legendre Functions             C
C                  Dimension ( ( LH + 1 )*( LH + 2 )/2, NTHPTS )     C
C                  P_l^m(X_j) = PA( L*(L+1)/2 + M + 1 , j )          C
C     DPA        : Derivatives of the above.                         C
C                                                                    C
C____________________________________________________________________C
C Output variables :-                                                C
C ================                                                   C
C  Double Precision                                                  C
C  ----------------                                                  C
C     VF	: Vector Function                                    C
C                 has dimensions ( NPHPTS, NTHPTS, 3 )               C
C____________________________________________________________________C
C
C*********************************************************************
      SUBROUTINE SHVECO ( L, M, ICS, IQST, VF, GAUX, PA, DPA,
     1                    NTHPTS, NPHPTS, LH, M0 )
      IMPLICIT NONE
C____________________________________________________________________C
C Variable declarations - Parameters ................................C
      INTEGER LH, NTHPTS, NPHPTS, IQST, L, M, ICS, M0
      DOUBLE PRECISION VF ( NPHPTS, NTHPTS, 3 ), GAUX( NTHPTS )
      DOUBLE PRECISION PA ( ( LH + 1 )*( LH + 2 )/2, NTHPTS ),
     1                 DPA ( ( LH + 1 )*( LH + 2 )/2, NTHPTS )
      DOUBLE PRECISION PI
      PARAMETER (PI=3.14159265358979312D0)
C____________________________________________________________________C
C Variable declarations - Working variables .........................C
      DOUBLE PRECISION X, S, THETA, ZERO, PHI, PINT, FA, SQRLL1
      INTEGER ITHETA, IOP, ITHREE, IP, IPHI
      PARAMETER ( ZERO=0.0d0, ITHREE=3 )
C____________________________________________________________________C
C Functions used :-
C____________________________________________________________________C
C START OF PROGRAM **************************************************C
C____________________________________________________________________C
C
C set VF array to zero.
      IOP = 0
      CALL CUBEOP ( VF, ZERO, NPHPTS, NTHPTS, ITHREE, IOP )
C     .
      IF ( M0.LT.1 .OR. M/M0*M0.NE.M ) THEN
        PRINT *,' Subroutine SHVECO.'
        PRINT *,' M = ', M,'  M0 = ', M0,' Program aborted.'
        STOP
      ENDIF
C     .
      IF ( ICS.NE.1 .AND. ICS.NE.2 ) THEN
        PRINT *,' Subroutine SHVECO.'
        PRINT *,' ICS = ', ICS,' Program aborted.'
        STOP
      ENDIF
C     .
      IF ( M.LT.0 .OR. M.GT.L ) THEN
        PRINT *,' Subroutine SHVECO.'
        PRINT *,' M = ',M,' L = ', L
        PRINT *,' Program aborted.'
        STOP
      ENDIF
C     .
      IF ( NTHPTS.LE.L ) THEN
         PRINT *,' Subroutine SHVECO.'
         PRINT *,' NTHPTS must be atleast ', L + 1
         PRINT *,' Program stopped.'
         STOP
      ENDIF
C
c     IF ( NPHPTS.LE.(2*M+1) ) THEN
c        PRINT *,' Subroutine SHVECO.'
c        PRINT *,' NPHPTS must be atleast ',(2*M+1)
c        PRINT *,' Program stopped.'
c        STOP
c     ENDIF
C
      IP   = L*(L+1)/2+M+1
      PINT = 2.0d0*PI/DBLE( NPHPTS*M0 )
      FA   = 1.0d0/SQRLL1( L )
C
C********************************************************************
C*********** First case IQST = 1 ... SCALOIDAL FUNCTION **************
C********************************************************************
      IF (IQST.EQ.1) THEN
C ............. start to loop around theta points ...................
        DO ITHETA = 1, NTHPTS
         X = GAUX( ITHETA )
         THETA = ACOS( X )
         S = DSIN( THETA )
C ............. First consider case M = 0
         IF (M.EQ.0) THEN
           DO IPHI = 1, NPHPTS
             VF ( IPHI, ITHETA, 1 ) = PA ( IP, ITHETA )
           ENDDO
         ELSE
           DO IPHI = 1, NPHPTS
              PHI = DBLE(IPHI - 1 )*PINT
              IF (ICS.EQ.1) THEN
                VF ( IPHI, ITHETA, 1 ) = PA ( IP, ITHETA )*
     1                                  DCOS( M*PHI )
              ELSE
                VF ( IPHI, ITHETA, 1 ) = PA ( IP, ITHETA )*
     1                                  DSIN( M*PHI )
              ENDIF
           ENDDO
         ENDIF
        ENDDO
C ............. ended looping around theta points ...................
        RETURN
      ENDIF
C
C********************************************************************
C*********** End of case IQST = 1 ( scaloidal ) **********************
C********************************************************************
C
C
C********************************************************************
C*********** First case IQST = 2 ... SPHEROIDAL FUNCTION *************
C********************************************************************
      IF (IQST.EQ.2) THEN
C ............. start to loop around theta points ...................
        DO ITHETA = 1, NTHPTS
         X = GAUX( ITHETA )
         THETA = ACOS( X )
         S = DSIN( THETA )
C ............. First consider case M = 0
        IF (M.EQ.0) THEN
           DO IPHI = 1, NPHPTS
              VF ( IPHI, ITHETA, 2 ) = DPA ( IP, ITHETA )*FA
           ENDDO
         ELSE
           DO IPHI = 1, NPHPTS
              PHI = DBLE(IPHI - 1 )*PINT
              IF (ICS.EQ.1) THEN
                VF ( IPHI, ITHETA, 2 ) = DPA ( IP, ITHETA )*FA*
     1                                  DCOS( M*PHI )
                VF ( IPHI, ITHETA, 3 ) = PA ( IP, ITHETA )*FA*
     1                             M*DSIN( M*PHI )*(-1.0d0/S)
              ELSE
                VF ( IPHI, ITHETA, 2 ) = DPA ( IP, ITHETA )*FA*
     1                                  DSIN( M*PHI )
                VF ( IPHI, ITHETA, 3 ) = PA ( IP, ITHETA )*FA*
     1                                  M*DCOS( M*PHI )/S
              ENDIF
           ENDDO
         ENDIF 
 
        ENDDO
C ............. ended looping around theta points ...................
        RETURN
      ENDIF
C
C********************************************************************
C*********** End of case IQST = 2 ( spheroidal ) *********************
C********************************************************************
C
C
C********************************************************************
C*********** First case IQST = 3 ... TOROIDAL FUNCTION ***************
C********************************************************************
      IF (IQST.EQ.3) THEN
C ............. start to loop around theta points ...................
        DO ITHETA = 1, NTHPTS
         X = GAUX( ITHETA )
         THETA = ACOS( X )
         S = DSIN( THETA )
C ............. First consider case M = 0
        IF (M.EQ.0) THEN
           DO IPHI = 1, NPHPTS
              VF ( IPHI, ITHETA, 3 ) = DPA ( IP, ITHETA )*FA
           ENDDO
         ELSE
           DO IPHI = 1, NPHPTS
              PHI = DBLE(IPHI - 1 )*PINT
              IF (ICS.EQ.1) THEN
                VF ( IPHI, ITHETA, 3 ) = DPA ( IP, ITHETA )*FA*
     1                                  DCOS( M*PHI )
                VF ( IPHI, ITHETA, 2 ) = PA ( IP, ITHETA )*FA*
     1                                  M*DSIN( M*PHI )/S
              ELSE
                VF ( IPHI, ITHETA, 3 ) = DPA ( IP, ITHETA )*FA*
     1                                  DSIN( M*PHI )
                VF ( IPHI, ITHETA, 2 ) = PA ( IP, ITHETA )*FA*
     1                             (-1.0d0)*M*DCOS( M*PHI )/S
              ENDIF
           ENDDO
         ENDIF
 
        ENDDO
C ............. ended looping around theta points ...................
        RETURN
      ENDIF
C
C********************************************************************
C*********** End of case IQST = 3 ( toroidal ) ***********************
C********************************************************************
C
      PRINT *,' IQST must be 1, 2 or 3 in SHVECO.'
      PRINT *,' Program aborted. '
      STOP
      END
C*********************************************************************
