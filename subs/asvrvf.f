C*********************************************************************
C subroutine Adapted Solution Vector to Radial Vector Function *******
C            -       -        -         -      -      -        *******
C Steve Gibbons Wed May 10 10:05:21 BST 2000                         C
C____________________________________________________________________C
C                                                                    C
C Input variables :-                                                 C
C ===============                                                    C
C  Integer                                                           C
C  -------                                                           C
C                                                                    C
C     NDCS      : Maximum distinct finite difference schemes         C
C                  stored by the array SVFDC.                        C
C                                                                    C
C     NTHPTS    : Number of theta points.                            C
C                                                                    C
C     NPHPTS    : Number of phi points.                              C
C                                                                    C
C     LH        : Maximum spherical harmonic degree, l.              C
C                 LH is the maximum l which can 'safely' be          C
C                 computed on this theta/phi grid. If a harmonic in  C
C                 the solution vector has l exceeding LH then the    C
C                 action is determined by IRES.                      C
C                                                                    C
C     IRES      : The resolution flag.                               C
C                 If IRES = 1, harmonics with l greater than LH      C
C                 or m greater than MMAX are simply ignored.         C
C                 If IRES = 2, the program aborts if LH or MMAX      C
C                 is exceeded.                                       C
C                                                                    C
C     NBN       : Number of bounding nodes.                          C
C                                                                    C
C     NFDCM     : Leading dimension of FDCM. At least (2*NBN+1)      C
C                                                                    C
C     ILNR      : Number of lowest radial grid node.                 C
C     IRNR      : Number of highest radial grid node.                C
C                                                                    C
C     NR        : Number of radial grid nodes in each function.      C
C                                                                    C
C     NDRVS     : Number of highest derivative for which             C
C                  coefficients are stored by the array SVFDC.       C
C                                                                    C
C     NDRVM     : Limit on NDRVS. Array bound for SVFDC.             C
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
C     MHP       : Array length ( * ) - atleast length NH             C
C                  Pointer array to finite difference coefficients.  C
C                  MHP( ih ) = is, which is the 4th index of         C
C                  array SVFDC - indicates f.d. scheme used.         C
C                                                                    C
C     MMAX      : Maximum spherical harmonic order, m.               C
C                                                                    C
C                                                                    C
C  Character                                                         C
C  ---------                                                         C
C                                                                    C
C     CHVMFF    : Velocity/magnetic field select.                    C
C                 Set to either 'VEL' or 'MAG'                       C
C                                                                    C
C  Double Precision                                                  C
C  ----------------                                                  C
C                                                                    C
C     V         : Solution vector. Dim ( * ) but length atleast      C
C                  NR*NH.                                            C
C                                                                    C
C     SVFDC     : Finite difference coefficient matrix.              C
C                  Dimension ( NFDCM, NR, NDRVM+1, NDCS ).           C
C                   Array is generated by the routine svfdcf         C
C                 See documentation for SVFDCF for details.          C
C                                                                    C
C     GAUX      : Cosines of the NTHPTS evaluated by the routine     C
C                  gauwts. Dimension ( NTHPTS ).                     C
C                                                                    C
C     PA        : Schmidt Normalised Legendre Functions              C
C                  Dim { ( LH + 1 )*( LH + 2 )/2 , NTHPTS }          C
C                  P_l^m(X_j) = PA( L*(L+1)/2 + M + 1 , j )          C
C     DPA       : Derivatives of the above.                          C
C                                                                    C
C     XARR      : Dim ( NR ). i^{th} element gives x_i.              C
C                                                                    C
C____________________________________________________________________C
C Output variables :-                                                C
C ================                                                   C
C  Double Precision                                                  C
C  ----------------                                                  C
C     RVF       : Radial Vector Function. An array of dimensions     C
C		( NR, NPHPTS, NTHPTS, 3) which contain the R, THETA  C
C 		   and PHI components of a VECTOR at each point      C
C                  ... i.e. VF ( IR, IPHI, ITHETA, 2 ) is the Theta  C
C                  compontent of the vector at (ir, iphi, itheta).   C
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
      SUBROUTINE ASVRVF( NDCS, NTHPTS, NPHPTS, LH, IRES, NBN, NFDCM,
     1                   ILNR, IRNR, NR, NDRVS, NDRVM, INARR, MHT,
     2                   MHL, MHM, MHP, MMAX, CHVMFF, V, SVFDC, GAUX,
     3                   PA, DPA, RVF, FTF1, FTF2, FTF3, XARR )
      IMPLICIT NONE
C____________________________________________________________________C
C Variable declarations - Parameters ................................C
      INTEGER NDCS, NTHPTS, NPHPTS, LH, IRES, NBN, NFDCM, NR,
     1        NDRVS, NDRVM, INARR( * ), MHT( * ), MHL( * ), MHM( * ),
     2        MHP( * ), MMAX, ILNR, IRNR
      DOUBLE PRECISION V( * ), XARR( NR ), 
     1                 SVFDC( NFDCM, NR, NDRVM+1, NDCS ),
     2                 RVF( NR, NPHPTS, NTHPTS, 3), GAUX( NTHPTS )
      DOUBLE PRECISION FTF1( 2*NPHPTS ),
     1                 FTF2( 2*NPHPTS ), FTF3( 2*NPHPTS ),
     2                 PA ( ( LH + 1 )*( LH + 2 )/2 , NTHPTS),
     3                 DPA ( ( LH + 1 )*( LH + 2 )/2 , NTHPTS)
      CHARACTER *(3) CHVMFF
C____________________________________________________________________C
C Variable declarations - Working variables .........................C
      INTEGER IOP, ISIGN, ITHREE, ITHETA, LENV, IPHI, IH,
     1        IP, ICS, L, M, NRR, NH, IPOL, ITOR,
     2        IS, IT, IHD, INDSAM, INDDIF, IR
      DOUBLE PRECISION ZERO, X, SINE, RAD, DPTERM, DTTERM,
     1        PTERM, QFAC, SFAC, TFAC, DERV( 2 ), TOL, D0F, D1F
      PARAMETER ( ZERO = 0.0d0, ITHREE = 3, TOL = 1.0d-9 )
C
C note that dtterm and dpterm are respectively
C the d/dtheta and d/dphi parts
C pterm is the term not differentiated w.r.t. theta or phi
C____________________________________________________________________C
C Functions used :-
      DOUBLE PRECISION SQRLL1
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
        PRINT *,' Subroutine ASVRVF.'
        PRINT *,' NR  = ', NR
        PRINT *,' NRR = ', NRR
        PRINT *,' Program stopped.'
        STOP
      ENDIF
C
      IF ( IRES.NE.1 .AND. IRES.NE.2 ) THEN
        PRINT *,' Subroutine ASVRVF.'
        PRINT *,' IRES = ', IRES
        PRINT *,' IRES must be either 1 or 2.'
        PRINT *,' Program stopped.'
        STOP
      ENDIF
C     .
C     . Let's check validity of CHVMFF
C     .
      IF ( CHVMFF.NE.'VEL' .AND. CHVMFF.NE.'vel' .AND.
     1     CHVMFF.NE.'Vel' .AND. CHVMFF.NE.'MAG' .AND.
     2     CHVMFF.NE.'mag' .AND. CHVMFF.NE.'Mag'  ) THEN
        PRINT *,' Subroutine ASVRVF.'
        PRINT *,' CHVMFF = ', CHVMFF
        PRINT *,' Program aborted.'
        STOP
      ENDIF
C     .
      IF ( CHVMFF.EQ.'VEL' .OR. CHVMFF.EQ.'vel' .OR.
     1     CHVMFF.EQ.'Vel' ) THEN
        IPOL = 1
        ITOR = 2
      ENDIF
C     .
      IF ( CHVMFF.EQ.'MAG' .OR. CHVMFF.EQ.'mag' .OR.
     1     CHVMFF.EQ.'Mag' ) THEN
        IPOL = 4
        ITOR = 5
      ENDIF
C 
C No need to check NPHPTS is power of 2 - this is done
C by FFTRLV ...
C ......... need to have MMAX < NPHPTS/2
C
      IF ( NTHPTS.LE.LH ) THEN
         PRINT *,' Subroutine ASVRVF.'
         PRINT *,' NTHPTS must be atleast ', LH + 1
         PRINT *,' Program stopped.'
         STOP
      ENDIF
C
      IF ( NPHPTS.LE.(2*MMAX+1) ) THEN
         PRINT *,' Subroutine ASVRVF.'
         PRINT *,' NPHPTS must be atleast ',(2*MMAX+1)
         PRINT *,' Program stopped.'
         STOP
      ENDIF
C
      IOP = 0
C ............. set RVF array to zero
      CALL QUADOP ( RVF, ZERO, NR, NPHPTS, NTHPTS, ITHREE, IOP )
C
C Begin loop around radial grid nodes
C
      DO IR = ILNR, IRNR
       RAD = XARR( IR )
       IF ( ABS( RAD ).LT.TOL ) THEN
         PRINT *,' Subroutine ASVRVF.'
         PRINT *,' RAD  = ', RAD
         PRINT *,' Division by zero will result.'
         PRINT *,' Program stopped.'
         STOP
       ENDIF
C
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
C
        DO IH = 1, NH
          IT = MHT( IH )
          IF ( IT.NE.IPOL .AND. IT.NE.ITOR ) GOTO 50
          L  = MHL( IH )
          IF ( MHM( IH ).LT.0 ) THEN
            ICS = 2
            M   = -MHM( IH )
          ELSE
            ICS = 1
            M   = MHM( IH )
          ENDIF
C         .
C         . Check the bounds on L and M
C         .
          IF ( L.GT.LH .OR. M.GT.MMAX ) THEN
            IF ( IRES.EQ.1 ) GOTO 50
            PRINT *,' Subroutine ASVRVF.'
            PRINT *,' IRES = ', IRES,' and '
            PRINT *,' L = ',L,' and LH = ',LH
            PRINT *,' M = ',M,' and MMAX = ',MMAX
            PRINT *,' Program stopped.'
            STOP
          ENDIF
C         .
          IS  = MHP( IH )
          IP  = L*(L+1)/2+M+1
C         .
C         . OK - let's take derivatives
C         .
          IHD = 1
C         .
          CALL ASVDR ( V, IR, IS, IH, NBN, IHD, NFDCM, NR, NDRVS,
     1                 NDRVM, DERV, INARR, SVFDC, NDCS )
C         .
          D0F = DERV( 1 )
          D1F = DERV( 2 )
C         .
          PTERM  = PA( IP , ITHETA )
          DTTERM = DPA( IP , ITHETA )/SQRLL1( L )
C         .
          IF ( ICS.EQ.1 ) THEN
            INDSAM = 2*M + 1
            INDDIF = 2*M + 2
            DPTERM = (-1.0d0)*M*PA( IP , ITHETA )/( SINE*SQRLL1( L ))
          ENDIF
          IF ( ICS.EQ.2 ) THEN
            INDDIF = 2*M + 1
            INDSAM = 2*M + 2
            DPTERM = DBLE(M)*PA( IP , ITHETA )/( SINE*SQRLL1( L ) )
          ENDIF
C         .
C         . First case of a poloidal harmonic
C         .
          IF ( IT.EQ.IPOL ) THEN
C           .
C           . Evaluate scaloidal and spheroidal values
C           .
            QFAC = DBLE( L*L + L )*D0F/RAD
            SFAC = SQRLL1( L )*(D0F/RAD + D1F)
C           .
            FTF1( INDSAM ) = FTF1( INDSAM ) + PTERM*QFAC
            FTF2( INDSAM ) = FTF2( INDSAM ) + DTTERM*SFAC
            FTF3( INDDIF ) = FTF3( INDDIF ) + DPTERM*SFAC
C           .
          ENDIF
C         .
C         . Now case of a toroidal harmonic
C         .
          IF ( IT.EQ.ITOR ) THEN
C           .
C           . Evaluate toroidal value
C           .
            TFAC = (-1.0d0)*SQRLL1( L )*D0F
C           .
            FTF2( INDDIF ) = FTF2( INDDIF ) - DPTERM*TFAC
            FTF3( INDSAM ) = FTF3( INDSAM ) + DTTERM*TFAC
C           .
          ENDIF
C         .
 50       CONTINUE
        ENDDO         
C ............. ended looping around Harmonics ......................
C ............... now perform Fourier Transforms on FTF1, FTF2, FTF3
C
        ISIGN = -1
        CALL FFTRLV ( FTF1, NPHPTS, ISIGN )
        CALL FFTRLV ( FTF2, NPHPTS, ISIGN )
        CALL FFTRLV ( FTF3, NPHPTS, ISIGN )
C
C ...................................................................
        DO IPHI = 1, NPHPTS
          RVF( IR, IPHI, ITHETA, 1 ) = FTF1( 2*IPHI - 1 )
          RVF( IR, IPHI, ITHETA, 2 ) = FTF2( 2*IPHI - 1 )
          RVF( IR, IPHI, ITHETA, 3 ) = FTF3( 2*IPHI - 1 )
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

