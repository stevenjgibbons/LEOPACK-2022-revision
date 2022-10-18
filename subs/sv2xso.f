C*********************************************************************
C subroutine Solution Vector 2 Xtra Special array Optimised **********
C            -        -      - -    -             -         **********
C Steve Gibbons Fri Nov  9 11:36:11 WET 2001                         C
C____________________________________________________________________C
C                                                                    C
C Input variables :-                                                 C
C ===============                                                    C
C  Integer                                                           C
C  -------                                                           C
C                                                                    C
C     NCMX      : Maximum number of components stored in XSV         C
C     NPHP      : Number of phi points.                              C
C     NTHP      : Number of theta points.                            C
C     NR        : Number of radial grid nodes.                       C
C     ILN       : First node to begin with.                          C
C     IRN       : Last node to end with.                             C
C                                                                    C
C     LH        : Maximum spherical harmonic degree, l.              C
C                 LH is the maximum l which can 'safely' be          C
C                 computed on this theta/phi grid. If a harmonic in  C
C                 the solution vector has l exceeding LH then the    C
C                 action is determined by IRES.                      C
C                                                                    C
C     INARR     : Integer array dimension ( * ).                     C
C                  Elements may be arbitrary except for              C
C                  INARR( 1 ) = IFORMF - flag for vector format.     C
C                                        See INDFUN                  C
C                  INARR( 2 ) = NRR. Must be consistent with NR.     C
C                  INARR( 3 ) = NH = total number of radial func.s   C
C                                                                    C
C     MHT       : Array length ( * ) - atleast length NH             C
C                                                                    C
C         MHT( IH ) = 1 --> harmonic is poloidal velocity            C
C         MHT( IH ) = 2 --> harmonic is toroidal velocity            C
C         MHT( IH ) = 3 --> harmonic is temperature.                 C
C         MHT( IH ) = 4 --> harmonic is poloidal magnetic field      C
C         MHT( IH ) = 5 --> harmonic is toroidal magnetic field      C
C                                                                    C
C     MHL       : Array length ( * ) - atleast length NH             C
C                  Sph. harm. degree, l.                             C
C     MHM       : Array length ( * ) - atleast length NH             C
C                  Sph. harm. order, m, for cos m phi dep.           C
C                 -Sph. harm. order, m, for sin m phi dep.           C
C                                                                    C
C     ICM       : Component of XSA to be filled.                     C
C                                                                    C
C     ICOMP     : Vector component to be filled.                     C
C                Key to ICOMP values.                                C
C                                                                    C
C        ICOMP =  1: Poloidal velocity scalar function               C
C        ICOMP =  2: Toroidal velocity scalar function               C
C        ICOMP =  3: Temperature scalar function                     C
C        ICOMP =  4: Poloidal mag. field scalar function             C
C        ICOMP =  5: Toroidal mag. field scalar function             C
C        ICOMP =  6: Radial component of (poloidal) velocity         C
C        ICOMP =  7: Radial component of (poloidal) magnetic field   C
C        ICOMP =  8: Theta component of poloidal velocity            C
C        ICOMP =  9: Theta component of poloidal magnetic field      C
C        ICOMP = 10: Phi component of poloidal velocity              C
C        ICOMP = 11: Phi component of poloidal magnetic field        C
C        ICOMP = 12: Theta component of toroidal velocity            C
C        ICOMP = 13: Theta component of toroidal magnetic field      C
C        ICOMP = 14: Phi component of toroidal velocity              C
C        ICOMP = 15: Phi component of toroidal magnetic field        C
C                                                                    C
C     M0        : Smallest non-zero wavenumber in solution.          C
C                                                                    C
C  Double Precision                                                  C
C  ----------------                                                  C
C                                                                    C
C     GAUX      : Cosines of the NTHPTS evaluated by the routine     C
C                  gauwts. Dimension ( NTHPTS ).                     C
C                                                                    C
C     PA        : Schmidt Normalised Legendre Functions              C
C                  Dim { ( LH + 1 )*( LH + 2 )/2 , NTHPTS }          C
C                  P_l^m(X_j) = PA( L*(L+1)/2 + M + 1 , j )          C
C     DPA       : Derivatives of the above.                          C
C                                                                    C
C     FTF       : Fourier transform work array. Dim (2*NPHP)         C
C                                                                    C
C     XARR      : Dim ( NR ). i^{th} element gives x_i.              C
C                                                                    C
C     V0        : Solution vector. Dim ( * ) but length atleast      C
C                  NR*NH. Zero^{th} derivatives.                     C
C     V1        : Solution vector. Dim ( * ) but length atleast      C
C                  NR*NH. First derivatives.                         C
C                                                                    C
C     XSA       : Array dim ( NCMX, NPHP, NTHP, NR )                 C
C____________________________________________________________________C
C
C*********************************************************************
      SUBROUTINE SV2XSO( NCMX, NPHP, NTHP, NR, ILN, IRN, LH, INARR, 
     1                   MHT, MHL, MHM, ICM, ICOMP, GAUX, PA, DPA,
     2                   FTF, XARR, V0, V1, XSA, M0 )
      IMPLICIT NONE
C____________________________________________________________________C
C Variable declarations - Parameters ................................C
      INTEGER          NCMX, NPHP, NTHP, NR, ILN, IRN, LH, INARR( * ),
     1                 MHT( * ), MHL( * ), MHM( * ), ICM, ICOMP, M0
      DOUBLE PRECISION GAUX( NTHP ), FTF( 2*NPHP ),
     1                 PA( ( LH + 1 )*( LH + 2 )/2, NTHP ),
     2                 DPA( ( LH + 1 )*( LH + 2 )/2, NTHP )
      DOUBLE PRECISION XARR( NR ), V0( * ), V1( * ),
     1                 XSA( NCMX, NPHP, NTHP, NR )
C____________________________________________________________________C
C Variable declarations - Working variables .........................C
      INTEGER          NRR, NH, IR, ITHP, IPHP, IOP, LENV, ISIGN,
     1                 L, M, ICS, IP, IND, INDSAM, MM,
     2                 INDDIF, IH, IT
      DOUBLE PRECISION DZERO, RAD, X, SINE, D0F, D1F, PTERM, DPTERM,
     1                 DTTERM, QFAC, SFAC, TFAC, DLOW
      PARAMETER ( DZERO = 0.0d0, DLOW = 1.0d-6 )
C____________________________________________________________________C
C Functions used :-
      INTEGER          INDFUN
      DOUBLE PRECISION SQRLL1
      EXTERNAL         SQRLL1, INDFUN
C____________________________________________________________________C
C START OF PROGRAM **************************************************C
C____________________________________________________________________C
C
C Check the validity of arguments .....
C
      NRR = INARR( 2 )
      NH  = INARR( 3 )
C
      IF ( NRR.NE.NR ) THEN
        PRINT *,' Subroutine SV2XSO.'
        PRINT *,' NR  = ', NR
        PRINT *,' NRR = ', NRR
        PRINT *,' Program stopped.'
        STOP
      ENDIF
C
      IF ( ICM.LT.1 .OR. ICM.GT.NCMX ) THEN
        PRINT *,' Subroutine SV2XSO.'
        PRINT *,' ICM = ', ICM,' NCMX = ', NCMX
        PRINT *,' Program aborted.'
        STOP
      ENDIF
C
      IF ( ICOMP.LT.1 .OR. ICOMP.GT.15 ) THEN
        PRINT *,' Subroutine SV2XSO.'
        PRINT *,' ICOMP = ', ICOMP
        PRINT *,' Program aborted.'
        STOP
      ENDIF
C
      IF ( M0.LT.0 ) THEN
        PRINT *,' Subroutine SV2XSO.'
        PRINT *,' M0 = ', M0
        PRINT *,' Program aborted.'
        STOP
      ENDIF
C
C Zero the appropriate component of XSA
C
      DO IR = ILN, IRN
        DO ITHP = 1, NTHP
          DO IPHP = 1, NPHP
            XSA( ICM, IPHP, ITHP, IR ) = DZERO
          ENDDO
        ENDDO
      ENDDO
C
      ISIGN = -1
      LENV  = 2*NPHP
      IOP   = 0
C
C ............... begin looping around radial grid nodes .............
      DO IR = ILN, IRN
        RAD = XARR( IR )
        IF ( RAD.LT.DLOW ) RAD = DLOW
C
C ............... begin looping around theta points ..................
        DO ITHP = 1, NTHP
C
          X = GAUX( ITHP )
          SINE = DSQRT( 1.0d0 - X*X )
C
          CALL VECOP( FTF, DZERO, LENV, IOP )
C
C Fill in function
C Begin loop around the harmonics
C
          DO IH = 1, NH     
C           .
            IT = MHT( IH )
            IF ( ICOMP.EQ.1  .AND. IT.NE.1 ) GOTO 50
            IF ( ICOMP.EQ.2  .AND. IT.NE.2 ) GOTO 50
            IF ( ICOMP.EQ.3  .AND. IT.NE.3 ) GOTO 50
            IF ( ICOMP.EQ.4  .AND. IT.NE.4 ) GOTO 50
            IF ( ICOMP.EQ.5  .AND. IT.NE.5 ) GOTO 50
            IF ( ICOMP.EQ.6  .AND. IT.NE.1 ) GOTO 50
            IF ( ICOMP.EQ.7  .AND. IT.NE.4 ) GOTO 50
            IF ( ICOMP.EQ.8  .AND. IT.NE.1 ) GOTO 50
            IF ( ICOMP.EQ.9  .AND. IT.NE.4 ) GOTO 50
            IF ( ICOMP.EQ.10 .AND. IT.NE.1 ) GOTO 50
            IF ( ICOMP.EQ.11 .AND. IT.NE.4 ) GOTO 50
            IF ( ICOMP.EQ.12 .AND. IT.NE.2 ) GOTO 50
            IF ( ICOMP.EQ.13 .AND. IT.NE.5 ) GOTO 50
            IF ( ICOMP.EQ.14 .AND. IT.NE.2 ) GOTO 50
            IF ( ICOMP.EQ.15 .AND. IT.NE.5 ) GOTO 50
C           .
            L  = MHL( IH )
            IF ( MHM( IH ).LT.0 ) THEN
              ICS = 2
              M   = -MHM( IH )
            ELSE
              ICS = 1
              M   = MHM( IH )
            ENDIF
            MM = M/M0
C           .
            IP  = L*(L+1)/2+M+1
C           .
            IND = INDFUN( IR, IH, INARR )
C           .
            D0F = V0( IND )
            D1F = V1( IND )
C           .
          PTERM  = PA( IP , ITHP )
          DTTERM = DPA( IP , ITHP )/SQRLL1( L )
C         .
          IF ( ICS.EQ.1 ) THEN
            INDSAM = 2*MM + 1
            INDDIF = 2*MM + 2
            DPTERM = (-1.0d0)*M*PA( IP, ITHP )/( SINE*SQRLL1( L ))
          ENDIF
          IF ( ICS.EQ.2 ) THEN
            INDDIF = 2*MM + 1
            INDSAM = 2*MM + 2
            DPTERM = DBLE(M)*PA( IP, ITHP )/( SINE*SQRLL1( L ) )
          ENDIF
C           .
C           . icomp between 1 and 5 --> straight scalar func.
C           .
            IF ( ICOMP.GE.1 .AND. ICOMP.LE.5 ) THEN
              FTF( INDSAM ) = FTF( INDSAM ) + PTERM*D0F
            ENDIF
C           .
C           . icomp = 6 ( 7 ) --> radial component
C           .                     poloidal vel (mag. field)
C           .
            IF ( ICOMP.EQ.6 .OR. ICOMP.EQ.7 ) THEN
              QFAC = DBLE( L*L + L )*D0F/RAD
              FTF( INDSAM ) = FTF( INDSAM ) + PTERM*QFAC
            ENDIF
C           .
C           . icomp = 8 ( 9 ) --> theta component
C           .                     poloidal vel (mag. field)
C           .
            IF ( ICOMP.EQ.8 .OR. ICOMP.EQ.9 ) THEN
              SFAC = SQRLL1( L )*(D0F/RAD + D1F)
              FTF( INDSAM ) = FTF( INDSAM ) + DTTERM*SFAC
            ENDIF
C           .
C           . icomp = 10 ( 11 ) --> phi component
C           .                     poloidal vel (mag. field)
C           .
            IF ( ICOMP.EQ.10 .OR. ICOMP.EQ.11 ) THEN
              SFAC = SQRLL1( L )*(D0F/RAD + D1F)
              FTF( INDDIF ) = FTF( INDDIF ) + DPTERM*SFAC
            ENDIF
C           .
C           . icomp = 12 ( 13 ) --> theta component
C           .                     toroidal vel (mag. field)
C           .
            IF ( ICOMP.EQ.12 .OR. ICOMP.EQ.13 ) THEN
              TFAC = SQRLL1( L )*D0F*(-1.0d0)
              FTF( INDDIF ) = FTF( INDDIF ) - DPTERM*TFAC
            ENDIF
C           .
C           . icomp = 14 ( 15 ) --> phi component
C           .                     toroidal vel (mag. field)
C           .
            IF ( ICOMP.EQ.14 .OR. ICOMP.EQ.15 ) THEN
              TFAC = SQRLL1( L )*D0F*(-1.0d0)
              FTF( INDSAM ) = FTF( INDSAM ) + DTTERM*TFAC
            ENDIF
C           .
 50       CONTINUE
          ENDDO
C
C End loop around the harmonics
C Ended filling in function
C Now perform Fourier transform on FTF
C
          CALL FFTRLV( FTF, NPHP, ISIGN ) 
C
          DO IPHP = 1, NPHP
            XSA( ICM, IPHP, ITHP, IR ) = FTF( 2*IPHP - 1 )
          ENDDO
C
       ENDDO
C ............. ended looping around theta points ...................
C
      ENDDO
C ............. ended looping around radial grid nodes ..............
C
      RETURN
      END
C*********************************************************************

