C*********************************************************************
C                                                                    C
C Steve Gibbons - Magnetic Field Component ANALysis version 1        C
C                 -        -     -         -                -        C
C                                                                    C
C Sat Apr  7 15:07:38 MET DST 2001                                   C
C                                                                    C
C*********************************************************************
      PROGRAM mfcanal1
      IMPLICIT NONE
C_____________________________________________________________________
C                                                                    C
C Variable declarations - array defining parameters.                 C
C                                                                    C
C____________________________________________________________________C
C
      INTEGER NHMAX, LHMAX, NDCS, NRMAX, ISVMAX, NBN, NCFM,
     1        NDRVM, NPHMAX, NTHMAX, NPMAX, NCMX
      PARAMETER ( NHMAX = 1700, LHMAX = 50, NDCS = LHMAX + 2,
     1            NRMAX = 100, ISVMAX = NHMAX*NRMAX, NBN = 3,
     2            NCFM  = 2*NBN + 1, NDRVM = 1, NCMX = 27,
     3            NTHMAX = 64, NPHMAX = 128 )
      PARAMETER ( NPMAX = (LHMAX+1)*(LHMAX+2)/2 )
C_____________________________________________________________________
C                                                                    C
C Variable declarations - main arrays.                               C
C                                                                    C
C____________________________________________________________________C
C
      INTEGER MHT( NHMAX ), MHL( NHMAX ), MHM( NHMAX ), MHP( NHMAX ),
     1        MHIBC( NDCS ), MHOBC( NDCS ),
     2        LARR( NDCS ), IWORK( NCFM )
C
      INTEGER INARR1( 3 ), INARR2( 3 )
C
      DOUBLE PRECISION SV( ISVMAX ), XARR( NRMAX ),
     1                 SVFDC( NCFM, NRMAX, NDRVM+1, NDCS ),
     2                 COEFM1( NCFM, NCFM ), W1( NCFM ),
     3                 COEFM2( NCFM, NCFM ), W2( NCFM )
C
      DOUBLE PRECISION 
     1                 V0( ISVMAX ), V1( ISVMAX ),
     2                 V2( ISVMAX )
C
      DOUBLE PRECISION GAUX( NTHMAX ), GAUW( NTHMAX ),
     1                 FTF( 2*NPHMAX )
C
      DOUBLE PRECISION PA( NPMAX, NTHMAX ), DPA( NPMAX, NTHMAX ),
     1                 XSF( NCMX, NPHMAX, NTHMAX, NRMAX )
C
C_____________________________________________________________________
C                                                                    C
C Variable declarations - scalars and small arrays                   C
C                                                                    C
C____________________________________________________________________C
C
      CHARACTER *(200) LINE
      CHARACTER *(80) FNAME, FNLOG
      INTEGER IWRTE, LU, I, NH, NCUDS, ILEN, NR, NR1,
     1        ILN, IRN, NDRVS, LULOG,
     2        LH, NTHP, NPHP, ICOMP, IHD, ICM
C
      DOUBLE PRECISION X1, X2, DLOW
      DOUBLE PRECISION PVKE, TVKE, PMKE, TMKE, DKE( 2 ), DINT, PMRE
      DOUBLE PRECISION BPBP, BRBP, BRBR, FSCF, PMFSF2, PMFTMF,
     1                 TMFSF2, FCCF
C_____________________________________________________________________
C                                                                    C
C Variable declarations - parameters and EXTERNAL declarations.      C
C                                                                    C
C____________________________________________________________________C
C
      PARAMETER ( IWRTE = 3, X1 = -1.0d0, X2 = 1.0d0, DLOW = 1.0d-8 )
C_____________________________________________________________________
C                  **************************************************C
C START OF PROGRAM **************************************************C
C____________________________________________________________________C
C
      LU    = 91
      LULOG = 92
 80   FORMAT(A)
C----------------------------------------------------------------
C                Enter filename of harmonics file.              |
C----------------------------------------------------------------
      PRINT *,' Enter name of harmonics file.'
 21   CONTINUE
      READ ( 5, 80 ) LINE
      IF ( LINE(1:1).EQ.'*' ) GOTO 21
C
      DO I = 1, 200
        IF ( LINE(I:I).EQ.' ' ) THEN
          ILEN = I-1
          GOTO 41
        ENDIF
      ENDDO
 41   CONTINUE
      FNAME = LINE(1:ILEN)
C
C Read in harmonics file
C First fill a diff. scheme with the no boundary condition
C property for the use of the scalar vector SV3
C
      MHIBC( 1 ) = 1
      MHOBC( 1 ) = 1
      LARR( 1 )  = 1
C
      NCUDS = 1
C     (ncuds - number of diff. schemes already in use).
      CALL HMFRD( NH, NHMAX, MHT, MHL, MHM, MHP, NCUDS, NDCS,
     1            MHIBC, MHOBC, LARR, LU, FNAME )
      INARR1( 3 ) = NH
C----------------------------------------------------------------
C                Enter filename of solution vector.             |
C----------------------------------------------------------------
      PRINT *,' Enter name of solution vector file.'
 22   CONTINUE
      READ ( 5, 80 ) LINE
      IF ( LINE(1:1).EQ.'*' ) GOTO 22
C
      DO I = 1, 200
        IF ( LINE(I:I).EQ.' ' ) THEN
          ILEN = I-1
          GOTO 42
        ENDIF
      ENDDO
 42   CONTINUE
      FNAME = LINE(1:ILEN)
C
C Read in solution vector file
C
      CALL SVFRD( INARR1, LU, NRMAX, V0, FNAME )
C
      NR1 = INARR1( 2 )
C
      IF ( INARR1( 1 ).NE.4 ) PRINT *,'INARR1(1) = ', INARR1( 1 )
C
      INARR2( 1 ) = 4
      INARR2( 2 ) = INARR1( 2 )
      INARR2( 3 ) = INARR1( 3 )
C
      CALL VECREA( V0, SV, INARR1, INARR2 )
C
C----------------------------------------------------------------
C                Enter filename of radial spacing vector.       |
C----------------------------------------------------------------
      PRINT *,' Enter name of radial spacing file.'
 23   CONTINUE
      READ ( 5, 80 ) LINE
      IF ( LINE(1:1).EQ.'*' ) GOTO 23
C
      DO I = 1, 200
        IF ( LINE(I:I).EQ.' ' ) THEN
          ILEN = I-1
          GOTO 43
        ENDIF
      ENDDO
 43   CONTINUE
      FNAME = LINE(1:ILEN)
C
C Read in radial spacing file
C
      CALL XARRRD( NR, NRMAX, XARR, LU, FNAME )
C
      IF ( NR1.NE.NR ) THEN
        PRINT *,' Solution vector and radial node '
        PRINT *,' file claim differing numbers of '
        PRINT *,' grid nodes. Program aborted.'
        STOP
      ENDIF
C----------------------------------------------------------------
C                Enter filename of output file                  |
C----------------------------------------------------------------
      PRINT *,' Enter name of output file.'
 29   CONTINUE
      READ ( 5, 80 ) LINE
      IF ( LINE(1:1).EQ.'*' ) GOTO 29
C
      DO I = 1, 200
        IF ( LINE(I:I).EQ.' ' ) THEN
          ILEN = I-1
          GOTO 44
        ENDIF
      ENDDO
 44   CONTINUE
      FNLOG = LINE(1:ILEN)
C
      CALL FOPEN ( LULOG, FNLOG, IWRTE )
C
C----------------------------------------------------------------
C          Calculate array of finite difference coefficients.   |
C----------------------------------------------------------------
C
      NDRVS = 1
      ILN   = 1
      IRN   = NR
      CALL SVFDCF( NR, NDCS, NBN, ILN, IRN, MHIBC, MHOBC,
     1             LARR, NCFM, NCFM, NDRVS, NDRVM, XARR,
     2             IWORK, SVFDC, COEFM1, COEFM2, W1, W2 )
C
C----------------------------------------------------------------
C
C First, find LH
C
      LH   = 0
      DO I = 1, NH
        IF ( MHL( I ).GT.LH ) LH   = MHL( I )
      ENDDO
      IF ( LH.GT.LHMAX ) THEN
         PRINT *,' LH = ', LH,' LHMAX = ', LHMAX
         PRINT *,' Program aborted.'
         STOP
      ENDIF
C
      CALL ONTPPF( LH, LH, NTHP, NPHP, NTHMAX, NPHMAX )
      NTHP = 3*NTHP/2
      IF ( NTHP.GT.NTHMAX ) THEN
        PRINT *,' NTHP   = ', NTHP
        PRINT *,' NTHMAX = ', NTHMAX
        PRINT *,' Program aborted.'
        STOP
      ENDIF
      CALL GAUWTS ( X1, X2, GAUX, GAUW, NTHP )
      CALL SCHNLA ( PA, DPA, GAUX, LH, NTHP )
C
      WRITE ( LULOG, * ) 'LH     =', LH
      WRITE ( LULOG, * ) 'NTHP   =', NTHP
      WRITE ( LULOG, * ) 'NPHP   =', NPHP
C
C----------------------------------------------------------------
C Calculate derivatives of radial functions
C----------------------------------------------------------------
C
      IHD = 1
      ILN =  1
      IRN = NR
      CALL CASVD2( SV, ILN, IRN, NBN, IHD, NCFM, NR, NDRVS,
     1             NDRVM, INARR2, NDCS, MHP, SVFDC, V0, V1, V2 )
C
C Calculate kinetic energy the old way
C
        PVKE   = 0.0d0
        TVKE   = 0.0d0
        PMKE   = 0.0d0
        PMRE   = 0.0d0
        TMKE   = 0.0d0
        DO I = 1, NH
          CALL SHKEER( I, NDCS, NR, INARR2, MHT, MHL, MHP, NBN,
     1               NDRVS, NDRVM, NCFM, SV, XARR, DKE, SVFDC )
          IF ( MHT( I ).EQ.1 ) PVKE = PVKE + DKE( 1 )
          IF ( MHT( I ).EQ.2 ) TVKE = TVKE + DKE( 1 )
          IF ( MHT( I ).EQ.4 ) PMKE = PMKE + DKE( 1 )
          IF ( MHT( I ).EQ.4 ) PMRE = PMRE + DKE( 2 )
          IF ( MHT( I ).EQ.5 ) TMKE = TMKE + DKE( 1 )
        ENDDO
        WRITE ( LULOG, * ) 'Old calculation: ... '
c       WRITE ( LULOG, 109 ) PVKE, TVKE, PMKE, TMKE
c       CALL FLUSH( LULOG )
c109  FORMAT(4(1PD16.7))
c     WRITE ( LULOG, * ) ' radial energy = ', PMRE
      PRINT *,' --------------------------- '
      WRITE ( LULOG, * ) 'Using SHKEER, the volume int. of magnetic '
      WRITE ( LULOG, * ) 'energy is ', PMKE + TMKE
C
      DO ICOMP = 1, 15
        ICM = ICOMP
        CALL SV2XSA( NCMX, NPHP, NTHP, NR, ILN, IRN, LH, INARR2,
     1               MHT, MHL, MHM, ICM, ICOMP, GAUX, PA, DPA,
     2               FTF, XARR, V0, V1, XSF )
      ENDDO
      CALL XSFANAL1( NCMX, NPHP, NTHP, NR, XSF )
C
C Calculate volume integral of ( magnetic energy )
C
      ICM = 20
      CALL XSVSCI( NCMX, NPHP, NTHP, NR, ICM, XARR, GAUW,
     1             XSF, DINT )
C
      WRITE ( LULOG, * ) ' --------------------------- '
      WRITE ( LULOG, * ) 'Using SV2XSA and XSVSCI, '
      WRITE ( LULOG, * ) 'the volume int. of magnetic '
      WRITE ( LULOG, * ) 'energy is ', DINT
      WRITE ( LULOG, * ) ' --------------------------- '
C
C Calculate volume integral of ( P(r,theta,phi)*2 )
C
      ICM = 22
      CALL XSVSCI( NCMX, NPHP, NTHP, NR, ICM, XARR, GAUW,
     1             XSF, DINT )
      PMFSF2 = DINT
C
C Calculate volume integral of ( T(r,theta,phi)*2 )
C
      ICM = 23
      CALL XSVSCI( NCMX, NPHP, NTHP, NR, ICM, XARR, GAUW,
     1             XSF, DINT )
      TMFSF2 = DINT
C
C Calculate volume integral of ( P(r,theta,phi)*T(r,theta,phi) )
C
      ICM = 24
      CALL XSVSCI( NCMX, NPHP, NTHP, NR, ICM, XARR, GAUW,
     1             XSF, DINT )
      PMFTMF = DINT
C
      IF ( PMFSF2.LT.DLOW .OR.
     1     TMFSF2.LT.DLOW      ) THEN
        PRINT *,' PMFSF2 = ', PMFSF2
        PRINT *,' TMFSF2 = ', TMFSF2
        PRINT *,' Program aborted.'
        STOP
      ELSE
        FSCF = PMFTMF/(PMFSF2*TMFSF2)
      ENDIF
C
C
c     WRITE ( LULOG, * ) 'integral of ( P(r,theta,phi)*2 ) = ',PMFSF2
c     WRITE ( LULOG, * ) 'integral of ( T(r,theta,phi)*2 ) = ',TMFSF2
c     WRITE ( LULOG, * ) 'integral of ( T*P              ) = ',PMFTMF
c     WRITE ( LULOG, * ) 'Correlation function             = ',FSCF
C
C Calculate volume integral of ( |B_rad| )
C
      ICM = 25
      CALL XSVSCI( NCMX, NPHP, NTHP, NR, ICM, XARR, GAUW,
     1             XSF, DINT )
      BRBR   = DINT
C
C Calculate volume integral of ( |B_phi| )
C
      ICM = 26
      CALL XSVSCI( NCMX, NPHP, NTHP, NR, ICM, XARR, GAUW,
     1             XSF, DINT )
      BPBP   = DINT
C
C Calculate volume integral of ( B_rad*B_phi )
C
      ICM = 27
      CALL XSVSCI( NCMX, NPHP, NTHP, NR, ICM, XARR, GAUW,
     1             XSF, DINT )
      BRBP   = DINT
C
      IF ( BRBR.LT.DLOW .OR.
     1     BPBP.LT.DLOW      ) THEN
        PRINT *,' BPBP = ', BPBP
        PRINT *,' BRBR = ', BRBR
        PRINT *,' Program aborted.'
        STOP
      ELSE
        FCCF = BRBP/( BRBR*BPBP )
      ENDIF
C
      WRITE ( LULOG, * ) 'integral of ( |B_rad| ) = ',BRBR
      WRITE ( LULOG, * ) 'integral of ( |B_phi| ) = ',BPBP
      WRITE ( LULOG, * ) 'integral of ( B_rad*B_phi ) = ',BRBP
      WRITE ( LULOG, * ) 'Correlation of ( B_rad*B_phi ) = ',FCCF
C
      CALL FCLOSE ( LULOG, FNLOG, 'Error' )
C
      STOP
      END
C*********************************************************************

C*********************************************************************
      SUBROUTINE VECREA( VEC1, VEC2, INARR1, INARR2 )
      IMPLICIT NONE
C
      INTEGER          INARR1( 3 ), INARR2( 3 )
      DOUBLE PRECISION VEC1( * ), VEC2( * )
C
      INTEGER          IH, IR, IND1, IND2, INDFUN,
     1                 NR1, NR2, NH1, NH2
      EXTERNAL         INDFUN
C
      NR1    = INARR1( 2 )
      NH1    = INARR1( 3 )
C
      NR2    = INARR2( 2 )
      NH2    = INARR2( 3 )
C
      IF ( NR1.NE.NR2 .OR. NH1.NE.NH2 ) THEN
        PRINT *,' Subroutine VECREA '
        PRINT *,' NR1 = ', NR1,' NH1 = ', NH1
        PRINT *,' NR2 = ', NR2,' NH2 = ', NH2
        PRINT *,' Program aborted.'
        STOP
      ENDIF
C
      DO IR = 1, NR1
        DO IH = 1, NH1
          IND1 = INDFUN( IR, IH, INARR1 )
          IND2 = INDFUN( IR, IH, INARR2 )
          VEC2( IND2 ) = VEC1( IND1 )
        ENDDO
      ENDDO
C
      RETURN
      END
C*********************************************************************

C*********************************************************************
      SUBROUTINE XSFANAL1( NCMX, NPHP, NTHP, NR, XSF )
      IMPLICIT NONE
C____________________________________________________________________C
C Variable declarations - Parameters ................................C
      INTEGER          NCMX, NPHP, NTHP, NR
      DOUBLE PRECISION XSF( NCMX, NPHP, NTHP, NR )
C____________________________________________________________________C
C Variable declarations - Working variables .........................C
      INTEGER          IR, ITHP, IPHP
      DOUBLE PRECISION DSCAL1, DSCAL2, DSCAL3, DSCAL4, DSCAL5,
     1                 DVRPOL, DMRPOL, DVTPOL, DMTPOL, DVPPOL, DMPPOL,
     2                                 DVTTOR, DMTTOR, DVPTOR, DMPTOR
      DOUBLE PRECISION POLME, TORME, POLKE, TORKE, POLTOR,
     1                 BRAD, BTHE, BPHI, VRAD, VTHE, VPHI,
     2                 TOTKE, TOTME
C____________________________________________________________________C
C START OF PROGRAM **************************************************C
C____________________________________________________________________C
C
      IF ( NCMX.LT.27 ) THEN
        PRINT *,' Subroutine XSFANAL1.'
        PRINT *,' NCMX = ', NCMX
        PRINT *,' Program aborted.'
        STOP
      ENDIF
C
      DO IR = 1, NR
        DO ITHP = 1, NTHP
          DO IPHP = 1, NPHP
            DSCAL1 = XSF(  1, IPHP, ITHP, IR )
            DSCAL2 = XSF(  2, IPHP, ITHP, IR )
            DSCAL3 = XSF(  3, IPHP, ITHP, IR )
            DSCAL4 = XSF(  4, IPHP, ITHP, IR )
            DSCAL5 = XSF(  5, IPHP, ITHP, IR )
            DVRPOL = XSF(  6, IPHP, ITHP, IR )
            DMRPOL = XSF(  7, IPHP, ITHP, IR )
            DVTPOL = XSF(  8, IPHP, ITHP, IR )
            DMTPOL = XSF(  9, IPHP, ITHP, IR )
            DVPPOL = XSF( 10, IPHP, ITHP, IR )
            DMPPOL = XSF( 11, IPHP, ITHP, IR )
            DVTTOR = XSF( 12, IPHP, ITHP, IR )
            DMTTOR = XSF( 13, IPHP, ITHP, IR )
            DVPTOR = XSF( 14, IPHP, ITHP, IR )
            DMPTOR = XSF( 15, IPHP, ITHP, IR )
C
            VRAD   = DVRPOL
            VTHE   = DVTPOL + DVTTOR
            VPHI   = DVPPOL + DVPTOR
C
            BRAD   = DMRPOL
            BTHE   = DMTPOL + DMTTOR
            BPHI   = DMPPOL + DMPTOR
C
            POLME  = DMRPOL*DMRPOL +
     1               DMTPOL*DMTPOL +
     2               DMPPOL*DMPPOL
C
            TORME  =
     1               DMTTOR*DMTTOR +
     2               DMPTOR*DMPTOR
C
            POLKE  = DVRPOL*DVRPOL +
     1               DVTPOL*DVTPOL +
     2               DVPPOL*DVPPOL
C
            TORKE  =
     1               DVTTOR*DVTTOR +
     2               DVPTOR*DVPTOR
C
            POLTOR = DMTTOR*DMTPOL +
     2               DMPTOR*DMPPOL
C
            TOTKE  = VRAD*VRAD + VTHE*VTHE + VPHI*VPHI
C
            TOTME  = BRAD*BRAD + BTHE*BTHE + BPHI*BPHI
C
            XSF( 16, IPHP, ITHP, IR ) = POLME*0.5d0
            XSF( 17, IPHP, ITHP, IR ) = TORME*0.5d0
            XSF( 18, IPHP, ITHP, IR ) = POLKE*0.5d0
            XSF( 19, IPHP, ITHP, IR ) = TORKE*0.5d0
            XSF( 20, IPHP, ITHP, IR ) = TOTME*0.5d0
            XSF( 21, IPHP, ITHP, IR ) = TOTKE*0.5d0
            XSF( 22, IPHP, ITHP, IR ) = DSQRT( POLME )
            XSF( 23, IPHP, ITHP, IR ) = DSQRT( TORME )
            XSF( 24, IPHP, ITHP, IR ) = POLME*TORME
            XSF( 25, IPHP, ITHP, IR ) = DSQRT( BRAD*BRAD )
            XSF( 26, IPHP, ITHP, IR ) = DSQRT( BPHI*BPHI )
            XSF( 27, IPHP, ITHP, IR ) = BRAD*BPHI
C
          ENDDO
        ENDDO
      ENDDO
C
      RETURN
      END
C*********************************************************************
C*********************************************************************
C subroutine HarMonic File ReaD **************************************
C            -  -     -    -  - **************************************
C Steve Gibbons Sat Nov 13 12:50:33 GMT 1999                         C
C____________________________________________________________________C
C                                                                    C
C Reads in the integer indices of spherical harmonic sets incl.      C
C the appropriate boundary conditions from a file.                   C
C It carefully assesses the boundary conditions and creates or       C
C appends the MHIBC, MHOBC and LARR arrays which must be sent to     C
C SVFDCF to calculate the finite difference coefficients.            C
C                                                                    C
C____________________________________________________________________C
C                                                                    C
C Input variables :-                                                 C
C ===============                                                    C
C  Integer                                                           C
C  -------                                                           C
C                                                                    C
C     NH        : Number of vector spherical harmonics. (Output)     C
C     NHMAX     : Maximum number of vector spherical harmonics.      C
C                                                                    C
C     MHT       : Array length ( * ) - atleast length NHMAX          C
C                                                                    C
C         MHT( IH ) = 1 --> harmonic is poloidal velocity            C
C         MHT( IH ) = 2 --> harmonic is toroidal velocity            C
C         MHT( IH ) = 3 --> harmonic is temperature.                 C
C         MHT( IH ) = 4 --> harmonic is poloidal magnetic field      C
C         MHT( IH ) = 5 --> harmonic is toroidal magnetic field      C
C                                                                    C
C     MHL       : Array length ( * ) - atleast length NHMAX          C
C                  Sph. harm. degree, l.                             C
C     MHM       : Array length ( * ) - atleast length NHMAX          C
C                  Sph. harm. order, m, for cos m phi dep.           C
C                 -Sph. harm. order, m, for sin m phi dep.           C
C                                                                    C
C     MHP       : Array length ( * ) - atleast length NHMAX          C
C                  Pointer array to finite difference coefficients.  C
C                  MHP( ih ) = is, which is the 4th index of         C
C                  array SVFDC - indicates f.d. scheme used.         C
C                                                                    C
C     NCUDS     : Number of currently used finite diff. schemes.     C
C     NDCS      : Number of distinct finite difference schemes.      C
C                                                                    C
C     MHIBC     : Dimension ( NDCS ). Governs behaviour at inner     C
C                  boundary for finite diff. scheme ( is )           C
C                                                                    C
C  MHIBC( is ) = 1 --> No condition is to be assumed at the bndry.   C
C  MHIBC( is ) = 2 --> Function must vanish at the bndry.            C
C  MHIBC( is ) = 3 --> First derivative must vanish at the bndry.    C
C  MHIBC( is ) = 4 --> Function must vanish at the bndry AND         C
C                       first derivative must vanish at the bndry.   C
C  MHIBC( is ) = 5 --> Function must vanish at the bndry AND         C
C                       second derivative must vanish at the bndry.  C
C  MHIBC( is ) = 6 --> r df/dr - f(r) = 0   at the bndry.            C
C  MHIBC( is ) = 7 --> r df/dr - l f(r) = 0 at the bndry.            C
C                        where L = MHL( ih )                         C
C                                                                    C
C     MHOBC     : Dimension ( NDCS ). Governs behaviour at outer     C
C                  boundary for finite diff. scheme ( is )           C
C                                                                    C
C  MHOBC( is ) = 1 --> No condition is to be assumed at the bndry.   C
C  MHOBC( is ) = 2 --> Function must vanish at the bndry.            C
C  MHOBC( is ) = 3 --> First derivative must vanish at the bndry.    C
C  MHOBC( is ) = 4 --> Function must vanish at the bndry AND         C
C                       first derivative must vanish at the bndry.   C
C  MHOBC( is ) = 5 --> Function must vanish at the bndry AND         C
C                       second derivative must vanish at the bndry.  C
C  MHOBC( is ) = 6 --> r df/dr - f(r) = 0   at the bndry.            C
C  MHOBC( is ) = 7 --> r df/dr + (l+1) f(r) = 0 at the bndry.        C
C                        where L = MHL( ih )                         C
C                                                                    C
C     LARR      : Spherical harmonic degree, L. Dim. ( NDCS ).       C
C                 Array to be passed to SVFDCF.                      C
C                                                                    C
C     LU        : Logical file unit number.                          C
C                                                                    C
C  Character                                                         C
C  ---------                                                         C
C                                                                    C
C     FNAME     : *(*) File name.                                    C
C                                                                    C
C                                                                    C
C____________________________________________________________________C
C
C*********************************************************************
      SUBROUTINE HMFRD( NH, NHMAX, MHT, MHL, MHM, MHP, NCUDS, NDCS,
     1                  MHIBC, MHOBC, LARR, LU, FNAME )
      IMPLICIT NONE
C____________________________________________________________________C
C Variable declarations - Parameters ................................C
      INTEGER NH, NHMAX, MHT( NHMAX ), MHL( NHMAX ), MHM( NHMAX ),
     1        MHP( NHMAX ), NCUDS, NDCS,
     1        MHIBC( NDCS ), MHOBC( NDCS ), LARR( NDCS ), LU
      CHARACTER *(*) FNAME
C____________________________________________________________________C
C Variable declarations - Working variables .........................C
      INTEGER IH, IWR, IS, IIBF, IOBF, ND1
C____________________________________________________________________C
C START OF PROGRAM **************************************************C
C____________________________________________________________________C
C
C Open file for reading
C
      IWR = 1
      CALL FOPEN ( LU, FNAME, IWR )
C
C Read in number of spherical harmonics
C
      READ ( LU, * ) NH
C
C Check value of NH
C
      IF ( NH.GT.NHMAX ) THEN
        PRINT *, ' Subroutine HMFRD.'
        PRINT *, ' From file, NH = ', NH
        PRINT *, ' NHMAX = ', NHMAX
        PRINT *, ' Program aborted.'
        STOP
      ENDIF
C
C     Now loop around the harmonics and read in their
C     properties
C
      DO IH = 1, NH
        READ ( LU, * ) MHT( IH ), MHL( IH ), MHM( IH ), IIBF, IOBF
C       .
C       . We now need to see whether we already have a
C       . suitable finite difference scheme for this outcome.
C       .
        ND1 = MIN( NCUDS, NDCS )
        DO IS = 1, ND1
         IF ( MHIBC( IS ).EQ.IIBF .AND. MHOBC( IS ).EQ.IOBF ) THEN
C          .
C          . Simple case where LARR needn't be applied
C          .
           IF ( IIBF.NE.7 .AND. IOBF.NE.7 ) THEN
             MHP( IH ) = IS
             GOTO 60
           ENDIF
C          .
C          . OK - check to see if LARR value is correct
C          .
           IF ( LARR( IS ).EQ.MHL( IH ) ) THEN
             MHP( IH ) = IS
             GOTO 60
           ENDIF
C          .
         ENDIF
        ENDDO
C       .
C       . OK - we don't already have this diff. scheme
C       .
        NCUDS = NCUDS + 1
        IF ( NCUDS.GT.NDCS ) THEN
          PRINT *, ' Subroutine HMFRD.'
          PRINT *, NDCS,' schemes are available.'
          PRINT *, ' and these have now been exhausted.'
          PRINT *, ' Program aborted.'
          STOP
        ENDIF
C       .
        MHIBC( NCUDS ) = IIBF
        MHOBC( NCUDS ) = IOBF
        IF ( IIBF.NE.7 .AND. IOBF.NE.7 ) THEN
          LARR( NCUDS ) = 0
        ELSE
          LARR( NCUDS ) = MHL( IH )
        ENDIF
        MHP( IH ) = NCUDS
C       .
 60     CONTINUE
C       .
      ENDDO
C
C Close file
C
      CALL FCLOSE ( LU, FNAME, 'Error closing file.' )
C
C Abort program if more finite difference schemes were
C requested by the harmonic sets than were allowed
C
      IF ( NCUDS.LT.NDCS ) THEN
        DO IS = NCUDS + 1, NDCS
          LARR( IS ) = -1
        ENDDO
      ENDIF
C
      RETURN
      END
C*********************************************************************

C*********************************************************************
C subroutine Solution Vector File ReaD *******************************
C            -        -      -    -  - *******************************
C Steve Gibbons Sat Nov 13 14:30:24 GMT 1999                         C
C____________________________________________________________________C
C                                                                    C
C Reads in a solution vector from a file.                            C
C IMPORTANT. The set of spherical harmonics must ALREADY BE KNOWN    C
C before calling SVFRD. If SVFRD reads in a file with NH different   C
C to the NH input in INARR( 3 ), an error will be reported.          C
C                                                                    C
C The only check on NR is that it is not greater than NRMAX.         C
C____________________________________________________________________C
C                                                                    C
C Input variables :-                                                 C
C ===============                                                    C
C  Integer                                                           C
C  -------                                                           C
C                                                                    C
C     INARR     : Array of integers. Dimension ( * )                 C
C                 At present the elements of INARR are as follows.   C
C                                                                    C
C                 INARR( 1 ) = IFORMF. Format flag.                  C
C                  This flag decides how the elements are arranged   C
C                  in the solution vector.                           C
C                                                                    C
C                  The current options are:-                         C
C                                                                    C
C                   IFORMF = 1, 3. INDFUN = ( IR - 1 )*NH + IH       C
C                   IFORMF = 2, 4. INDFUN = ( IH - 1 )*NR + IR       C
C                                                                    C
C  where IR and IH are the current grid node and harmonic resp.      C
C  and NR and NH are the total numbers of nodes and harmonics        C
C  in the solution vector.                                           C
C                                                                    C
C                 INARR( 2 ) = NR. Number of radial grid nodes.      C
C                 INARR( 3 ) = NH. Number of harmonics in sol. vect. C
C                                                                    C
C     LU        : Logical file unit number.                          C
C                                                                    C
C     NRMAX     : Maximum permitted radial grid nodes.               C
C                                                                    C
C  Double Precision                                                  C
C  ----------------                                                  C
C                                                                    C
C     SV        : Dim ( * ) but length atleast NR*NH.                C
C                  Solution vector defined by INARR.                 C
C                                                                    C
C  Character                                                         C
C  ---------                                                         C
C                                                                    C
C     FNAME     : *(*) File name.                                    C
C____________________________________________________________________C
C
C*********************************************************************
      SUBROUTINE SVFRD( INARR, LU, NRMAX, SV, FNAME )
      IMPLICIT NONE
C____________________________________________________________________C
C Variable declarations - Parameters ................................C
      INTEGER INARR( * ), LU, NRMAX
      CHARACTER *(*) FNAME
      DOUBLE PRECISION SV( * )
C____________________________________________________________________C
C Variable declarations - Working variables .........................C
      INTEGER I, IWR, ILEN, IFORMF, NR, NH, NH1, IFORM
C____________________________________________________________________C
C START OF PROGRAM **************************************************C
C____________________________________________________________________C
C
      NH     = INARR( 3 )
C
C Open file for reading
C
      IWR = 1
      CALL FOPEN ( LU, FNAME, IWR )
C
C Read iformf, nr, nh, iform
C  
       READ ( LU, * ) IFORMF, NR, NH1, IFORM
C
C Check that value of IFORM is legal
C
      IF ( IFORM.NE.1 ) THEN
        PRINT *,' Subroutine SVFRD.'
        PRINT *,' IFORM = ', IFORM
        PRINT *,' Currently, 1 is the only permissible value.'
        PRINT *,' Program aborted.'
        STOP
      ENDIF
C
C Check that value of NH is legal
C
      IF ( NH.NE.NH1 ) THEN
        PRINT *,' Subroutine SVFRD.'
        PRINT *,' INARR( 3 ) = ', NH
        PRINT *,' File contains ',NH1,' harmonics.'
        PRINT *,' Load in correct indices before calling SVFRD.'
        PRINT *,' Program aborted.'
        STOP
      ENDIF
C
C Check that value of NR is legal
C
      IF ( NR.GT.NRMAX ) THEN
        PRINT *,' Subroutine SVFRD.'
        PRINT *,' In file, NR = ', NR
        PRINT *,' NRMAX = ', NRMAX
        PRINT *,' Program aborted.'
        STOP
      ENDIF
C
      INARR( 1 ) = IFORMF
      INARR( 2 ) = NR
C
      ILEN   = NR*NH
C
C OK, so read SV values ...
C
      IF ( IFORM.EQ.1 ) READ ( LU, 41 ) ( SV( I ), I = 1, ILEN )
C
C Close file
C
      CALL FCLOSE ( LU, FNAME, 'Error closing file.' )
C
 41   FORMAT(5(1PD16.7))
      RETURN
      END
C*********************************************************************

C*********************************************************************
C subroutine X value ARRay ReaD **************************************
C            -       ---   -  - **************************************
C Steve Gibbons Fri Nov 12 07:58:51 GMT 1999                         C
C____________________________________________________________________C
C                                                                    C
C Reads in the XARR array of abscissae from a file.                  C
C____________________________________________________________________C
C                                                                    C
C Input variables :-                                                 C
C ===============                                                    C
C  Integer                                                           C
C  -------                                                           C
C                                                                    C
C     NR        : Number of radial grid nodes.                       C
C                 Note that NR is not checked for correspondence to  C
C                 any other value - merely for being not greater     C
C                 than NRMAX.                                        C
C                                                                    C
C     NRMAX     : Maximum number of radial grid nodes.               C
C                                                                    C
C     LU        : Logical file unit number.                          C
C                                                                    C
C  Double Precision                                                  C
C  ----------------                                                  C
C                                                                    C
C     XARR      : Dim ( * ) but length atleast NR. Location of       C
C                  radial grid nodes.                                C
C                                                                    C
C  Character                                                         C
C  ---------                                                         C
C                                                                    C
C     FNAME     : *(*) File name.                                    C
C____________________________________________________________________C
C
C*********************************************************************
      SUBROUTINE XARRRD( NR, NRMAX, XARR, LU, FNAME )
      IMPLICIT NONE
C____________________________________________________________________C
C Variable declarations - Parameters ................................C
      INTEGER NR, NRMAX, LU
      CHARACTER *(*) FNAME
      DOUBLE PRECISION XARR( * )
C____________________________________________________________________C
C Variable declarations - Working variables .........................C
      INTEGER I, IWR, IFORM
C____________________________________________________________________C
C START OF PROGRAM **************************************************C
C____________________________________________________________________C
C
C Open file for reading
C
      IWR = 1
      CALL FOPEN ( LU, FNAME, IWR )
C
C Read number of radial grid nodes
C
      READ ( LU, * ) NR, IFORM
C
C Check NR
C
      IF ( NR.GT.NRMAX ) THEN
        PRINT *,' Subroutine XARRRD.'
        PRINT *,' From file, NR = ', NR
        PRINT *,' NRMAX = ', NRMAX
        PRINT *,' Program aborted.'
        STOP
      ENDIF
C
C Check that value of IFORM is legal
C
      IF ( IFORM.NE.1 ) THEN
        PRINT *,' Subroutine XARRRD.'
        PRINT *,' From file, IFORM = ', IFORM
        PRINT *,' Currently, 1 is the only permissible value.'
        PRINT *,' Program aborted.'
        STOP
      ENDIF
C
C OK, so read X values ...
C
      IF ( IFORM.EQ.1 ) READ ( LU, 41 ) ( XARR( I ), I = 1, NR )
C
C Close file
C
      CALL FCLOSE ( LU, FNAME, 'Error closing file.' )
C
 41   FORMAT(5(1PD16.7))
      RETURN
      END
C*********************************************************************
C********************************************************************
C SUBROUTINE File OPEN **********************************************
C            -    ---- **********************************************
C Steve Gibbons 14.4.97 (Adapted from Dan Gordon's Code)            C
C Routine modified 15.3.99
C___________________________________________________________________C
C Opens a file with number LU, name FNAME, access OACCES  and       C
C a flag IRW to indicate whether the file is to be read or written  C
C to. ( IRW=1 ==> read only, IRW=2 ==> write but only if the file   C
C doesn't already exist, IRW=3 ==> write regardless of whether file C
C exists or not.)						    C
C___________________________________________________________________C
C Input Variables :-						    C
C ===============   						    C
C  Integer							    C
C  -------							    C
C     LU	: File number					    C
C     IRW	: Read / Write Flag 				    C
C                  = 1 for read only		                    C
C                  = 2 for write (provided that file doesn't exist. C
C                  = 3 for write (regardless of existence of file.  C
C                  = 4 for append status.                           C
C  Character							    C
C  ---------							    C
C     FNAME	: File name					    C
C___________________________________________________________________C
C Working Variables :-						    C
C =================   						    C
C  Character							    C
C  ---------							    C
C     OACCES 	: Access flag - should be set to 'OLD' for read     C
C                          and 'UNKNOWN' for write                  C
C     CONTYN    : For a yes/no to IWR = 2 option.		    C
C     LABEL     : Null string to pass into FNAMER option.	    C
C  Logical							    C
C  -------							    C
C     LEXIST	: Existence of file. File present <==> LEXIST=.TRUE.C
C___________________________________________________________________C
C Subroutines Used :-                                               C
C ================                                                  C
C     FNAMER	: For the case of IRW = 2; trying to write to an    C
C		   existing file. Used if alternative filename is   C
C                   asked for.					    C
C___________________________________________________________________C
C
C********************************************************************
      SUBROUTINE FOPEN ( LU, FNAME, IRW)
      IMPLICIT NONE
C___________________________________________________________________C
C Variable declarations - Parameters ...............................C
      INTEGER LU, IRW
      CHARACTER *(*) FNAME
C___________________________________________________________________C
C Variable declarations - Working Variables ........................C
      LOGICAL LEXIST
      CHARACTER *(7) OACCES
      CHARACTER *(1) CONTYN
      CHARACTER *(1) LABEL
C___________________________________________________________________C
      IF ( LU.EQ.0 ) THEN
         PRINT *,' Subroutine FOPEN'
         PRINT *,' I bet you ve forgotten to set LU ...??'
         PRINT *,' Think again and come back when you have'
         PRINT *,' remembered that LU must be a non zero integer!'
         PRINT *,' See you later. Bye for now!!'
         STOP
      ENDIF
C************************
C temporary code : SJG Thu Jun  1 08:07:06 BST 2000
C The Linux compiler will not allow file opening with
C lu.ge.100, so the following line should prevent it.
C 
      IF ( LU.GT.99 ) THEN
         PRINT *,' Subroutine FOPEN'
         PRINT *,' LU = ', LU,' too large.'
         PRINT *,' Program aborted.'
         STOP
      ENDIF
C************************
C--------------
 600  CONTINUE
      INQUIRE (FILE=FNAME, EXIST=LEXIST)
C Case of read only - LEXIST must be = .TRUE.
      IF ( IRW.EQ.1 ) THEN
         OACCES = 'OLD'
         IF ( .NOT. LEXIST ) THEN
            PRINT *,' Subroutine FOPEN. You are trying to open an'
            PRINT *,' old file which does not exist.'
            PRINT *,' Filename = ',FNAME,' LU= ',LU
            PRINT *,' Program aborted.'
            STOP
         ELSE
            GOTO 500
         ENDIF
      ENDIF
C Case of write to file provided that it doesn't exist
      IF ( IRW.EQ.2 ) THEN
         OACCES = 'UNKNOWN'
         IF ( LEXIST ) THEN
            PRINT *, ' Subroutine FOPEN. You are trying to write'
            PRINT *, ' to an existing file with IRW set to 2.'
            PRINT *,' Filename = ',FNAME,' LU= ',LU
            PRINT *,' Do you wish to give an alternative FNAME?'
            PRINT *,' Type y or n.'
            READ ( 5, 267) CONTYN
 267         FORMAT (A)
            IF (CONTYN.NE.'y'.AND.CONTYN.NE.'Y') THEN
               PRINT *, ' Program Aborted.'
               STOP
            ELSE
               LABEL=' '
               CALL FNAMER ( FNAME, LABEL )
               GOTO 600
            ENDIF
         ELSE
            GOTO 500
         ENDIF
      ENDIF
C Case of write to file regardless of the existence of file.
      IF ( IRW.EQ.3 ) THEN
         OACCES = 'UNKNOWN'
         GOTO 500
      ENDIF
C Treat appendment case
      IF ( IRW.EQ.4 ) THEN
         OACCES = 'UNKNOWN'
         OPEN ( UNIT=LU , FILE=FNAME , STATUS=OACCES,
     1          ACCESS='APPEND', ERR=999 )
         RETURN
      ENDIF
C___________________________________________________________________C
C All the IRW cases as of 14.4.97 have now been covered.
      PRINT *,' Subroutine FOPEN. IRW must be set to 1, 2, 3 or 4.'
      PRINT *,' Program aborted.'
      STOP

 500  CONTINUE
      OPEN ( UNIT=LU , FILE=FNAME , STATUS=OACCES, ERR=999 )
      RETURN

 999  PRINT *,' Subroutine FOPEN. Error in opening file ',FNAME
      STOP

      END
C********************************************************************
C___________________________________________________________________C
C*********************************************************************
C subroutine Solution Vector Finite Difference Coefficients Form *****
C            -        -      -      -          -            -    *****
C Steve Gibbons Fri Oct 22 09:33:36 BST 1999                         C
C____________________________________________________________________C
C                                                                    C
C  If XARR is an array of length ( NR ) such that the j^{th}         C
C  element is the value of x_j, then SVFDCF builds an array          C
C  SVFDC of dimension ( NFDCM, NR, NDRVM+1, NDCS ) such that for a   C
C  node number, j, the ND^{th} derivative of radial function         C
C  given f_{IH} ( x ) will be given by                               C
C                                                                    C
C  f_{IH}^{ND}( x_j ) =                                              C
C         \sum_{i=LN}^{RN} SVFDC ( IRAD, j, ND+1, K ) f_{IH} ( x_i ) C
C                                                                    C
C  where LN (the left node)  = MAX( NLMR, j - NBN ) and              C
C        RN (the right node) = MIN( NRMC, j + NBN ),                 C
C                                                                    C
C  IRAD = i - j + NBN + 1 and K = MHP( ih ).                         C
C                                                                    C
C  NDCS is the number of distinct sets of coefficients required for  C
C  different types of harmonics (in generally will be considerably   C
C  smaller than the number of harmonics, NH).                        C
C                                                                    C
C  The arrays MHIBC and MHOBC instruct SVFDCF how to manipulate      C
C  the finite difference coefficients at the boundaries.             C
C                                                                    C
C  MHIBC( ih ) = 1 --> No condition is to be assumed at the bndry.   C
C  MHIBC( ih ) = 2 --> Function must vanish at the bndry.            C
C  MHIBC( ih ) = 3 --> First derivative must vanish at the bndry.    C
C  MHIBC( ih ) = 4 --> Function must vanish at the bndry AND         C
C                       first derivative must vanish at the bndry.   C
C  MHIBC( ih ) = 5 --> Function must vanish at the bndry AND         C
C                       second derivative must vanish at the bndry.  C
C  MHIBC( ih ) = 6 --> r df/dr - f(r) = 0   at the bndry.            C
C  MHIBC( ih ) = 7 --> r df/dr - l f(r) = 0 at the bndry.            C
C                        where L = LARR( ih )                        C
C                                                                    C
C  Similarly, at the outer boundary:-                                C
C                                                                    C
C  MHOBC( ih ) = 1 --> No condition is to be assumed at the bndry.   C
C  MHOBC( ih ) = 2 --> Function must vanish at the bndry.            C
C  MHOBC( ih ) = 3 --> First derivative must vanish at the bndry.    C
C  MHOBC( ih ) = 4 --> Function must vanish at the bndry AND         C
C                       first derivative must vanish at the bndry.   C
C  MHOBC( ih ) = 5 --> Function must vanish at the bndry AND         C
C                       second derivative must vanish at the bndry.  C
C  MHOBC( ih ) = 6 --> r df/dr - f(r) = 0   at the bndry.            C
C  MHOBC( ih ) = 7 --> r df/dr + (l+1) f(r) = 0 at the bndry.        C
C                        where L = LARR( ih )                        C
C                                                                    C
C  The elements of this array are filled in from j = NLMR            C
C  to j = NRMR ( number of the left most node and number of the      C
C  right most node ) - other rows are left unreferred to.            C
C  This is incase a higher order derivative is required for          C
C  central nodes than boundary nodes; in which case SVFDCF must      C
C  be called for the remaining nodes with modified parameters.       C
C                                                                    C
C____________________________________________________________________C
C                                                                    C
C Input variables :-                                                 C
C ===============                                                    C
C                                                                    C
C  Integer                                                           C
C  -------                                                           C
C                                                                    C
C     NR        : Number of grid nodes.                              C
C     NDCS       : Number of distinct differencing coeff.s           C
C                  represented in SVFDC.                             C
C                                                                    C
C     NBN       : Maximum number of nodes on either side for         C
C                  central differencing.                             C
C                                                                    C
C     NLMR      : This is the lowest j for which the terms are       C
C                  calculated for SVFDC( i, j, ND+1, K )             C
C     NRMR      : This is the highest j for which the terms are      C
C                  calculated for SVFDC( i, j, ND+1, K )             C
C                                                                    C
C     MHIBC     : Dimension ( NDCS ). See above.                     C
C     MHOBC     : Dimension ( NDCS ). See above.                     C
C     LARR      : Spherical harmonic degree, L. Dim. ( NDCS ).       C
C                 This value is only useful for harmonics            C
C                 when calculating derivatives for magnetic fields.  C
C                 If LARR( ih ) = -1, the harmonic is ignored        C
C                 completely.                                        C
C                                                                    C
C     NCFM      : Leading order of working coefficient matrix.       C
C                 Must be atleast (2*NBN + 1) where NBN is the       C
C                 maximum number of nodes on either side of the      C
C                 central node.                                      C
C     NFDCM     : Leading order of the array SVFDC.                  C
C                 This must be atleast (2*NBN + 1)                   C
C     NDRVS     : Number of derivatives required.                    C
C                  This will be limited by the available bandwidth.  C
C                                                                    C
C                  Let NLCS = NLMR - 1    and let                    C
C                      NRCS = NR - NRMR                              C
C                                                                    C
C                  Now, let I = MIN( NLCS, NRCS) + NBN               C
C                                                                    C
C                  then NDRVS must be no greater than I.             C
C                  This is checked for.                              C
C     NDRVM     : Maximum number of derivatives allowed.             C
C                                                                    C
C     IWORK     : Integer work array. Dimension ( NCFM )             C
C                                                                    C
C  Double Precision                                                  C
C  ----------------                                                  C
C                                                                    C
C     XARR      : Array of dimension (  NR  ).                       C
C                  XARR( i ) contains the value of x or r at the     C
C                   i^{th} radial grid node.                         C
C                                                                    C
C     SVFDC     : Finite difference coefficient matrix.              C
C                  Dimension ( NFDCM, NR, NDRVM+1, NDCS ).           C
C                                                                    C
C     COEFM1    : Coefficient work array. Dimension ( NCFM, NCFM )   C
C     COEFM2    : Coefficient work array. Dimension ( NCFM, NCFM )   C
C     WORK1     : Coefficient work array. Dimension ( NCFM )         C
C     WORK2     : Coefficient work array. Dimension ( NCFM )         C
C____________________________________________________________________C
C
C*********************************************************************
      SUBROUTINE SVFDCF( NR, NDCS, NBN, NLMR, NRMR, MHIBC, MHOBC,
     1                   LARR, NCFM, NFDCM, NDRVS, NDRVM, XARR,
     2                   IWORK, SVFDC, COEFM1, COEFM2, WORK1, WORK2 )
      IMPLICIT NONE
C____________________________________________________________________C
C Variable declarations - Parameters ................................C
      INTEGER NR, NDCS, NBN, NLMR, NRMR, MHIBC( NDCS ), MHOBC( NDCS ),
     1        LARR( NDCS ), NCFM, NFDCM, NDRVS, NDRVM, 
     1        IWORK( NCFM )
      DOUBLE PRECISION XARR( NR ), SVFDC( NFDCM, NR, NDRVM+1, NDCS ),
     1                 COEFM1( NCFM, NCFM ), COEFM2( NCFM, NCFM ),
     2                 WORK1( NCFM ), WORK2( NCFM )
C____________________________________________________________________C
C Variable declarations - Working variables .........................C
C
      INTEGER NDER, IRAD, NLCS, NRCS, NLN, NRN, IDCS, L,
     1        NNDS, INDS, I, INODE, NSNIB, NSNOB, NALF, NARF,
     2        IIBC, IOBC, ND1
C
C nsnib is the number of special nodes on the inner boundary
C nsnob is the number of special nodes on the outer boundary
C
      DOUBLE PRECISION DZERO, X0, EMMULT, FAC
      PARAMETER ( DZERO = 0.0d0 )
C
      LOGICAL OCHNGE
C
C ochnge is .TRUE. when the boundary conditions
C play a part in the finite difference coefficients
C and .FALSE. otherwise
C____________________________________________________________________C
C START OF PROGRAM **************************************************C
C____________________________________________________________________C
C
C     .
C     . Check the values of integers ...
C     .
      IF ( NDRVS.GT.NDRVM ) THEN
         PRINT *,' Subroutine SVFDCF.'
         PRINT *,' NDRVS  = ',NDRVS
         PRINT *,' NDRVM  = ',NDRVM
         STOP
      ENDIF
C     .
      NLCS = NLMR - 1
      NRCS = NR - NRMR
C     . 
C     . Check that sufficient points are allowed
C     . for the derivatives ...
C     . 
      IF ( (NR-1).LT.(NBN+1) ) THEN
         PRINT *,' Subroutine SVFDCF.'
         PRINT *,' NBN  = ', NBN
         PRINT *,' NLMR = ', NLMR
         PRINT *,' NRMR = ', NRMR
         PRINT *,' Insufficient nodes for differencing.'
         STOP
      ENDIF
C     . 
      IF ( NLCS.LT.0 ) THEN
         PRINT *,' Subroutine SVFDCF.'
         PRINT *,' NLMR = ', NLMR
         STOP
      ENDIF
C     . 
      IF ( NRCS.LT.0 ) THEN
         PRINT *,' Subroutine SVFDCF.'
         PRINT *,' NRMR = ', NRMR
         STOP
      ENDIF
C     .
      I = MIN( NLCS, NRCS) + NBN
C     .
      IF ( NDRVS.GT.I ) THEN
         PRINT *,' Subroutine SVFDCF.'
         PRINT *,' You have requested deriv.s to order ',NDRVS
         PRINT *,' At one node, you have only', I
         PRINT *,' side nodes to use for differencing.'
         STOP
      ENDIF
C     .
      I = 2*NBN + 1
C     .
      IF ( NCFM.LT.I ) THEN
         PRINT *,' Subroutine SVFDCF.'
         PRINT *,' NCFM = ', NCFM
         PRINT *,' NBN  = ', NBN
         STOP
      ENDIF
C     .
      IF ( NFDCM.LT.I ) THEN
         PRINT *,' Subroutine SVFDCF.'
         PRINT *,' NFDCM = ', NFDCM
         PRINT *,' NBN  = ', NBN 
         STOP
      ENDIF
C     .
      IF ( NLMR.GT.NRMR ) THEN
         PRINT *,' Subroutine SVFDCF.'
         PRINT *,' NLMR = ', NLMR
         PRINT *,' NRMR = ', NRMR
         STOP
      ENDIF
C     .
C     . Whew ... all input parameters seem to be o.k.
C     . Loop around the requested nodes
C     .
      DO IRAD = NLMR, NRMR
C      .
C      . Loop around the different 'harmonic forms'
C      .
       DO IDCS = 1, NDCS
C
C Check for flag to ignore this harmonic
C
        IF ( LARR( IDCS ).EQ.-1 ) GOTO 50
C
C We now need to check which boundary condition is
C required. Make sure that it is valid.
C
        IF ( MHIBC( IDCS ).LT.1 .AND. MHIBC( IDCS ).GT.7 ) THEN
          PRINT *,' Subroutine SVFDCF.'
          PRINT *,' MHIBC(',IDCS,') = ', MHIBC( IDCS )
          STOP
        ENDIF
C
C O.k. inner b.c. is fine. Now need to 
C see how many points this effects.
C
        IF ( MHIBC( IDCS ).EQ.1 ) THEN
          NSNIB = 0
        ENDIF
C
        IF ( MHIBC( IDCS ).EQ.2 .OR. MHIBC( IDCS ).EQ.3 .OR.
     1       MHIBC( IDCS ).EQ.6 .OR. MHIBC( IDCS ).EQ.7     ) THEN
          NSNIB = 1
        ENDIF
C
        IF ( MHIBC( IDCS ).EQ.4 .OR. MHIBC( IDCS ).EQ.5 ) THEN
          NSNIB = 2
        ENDIF
C
        IF ( MHOBC( IDCS ).LT.1 .AND. MHOBC( IDCS ).GT.7 ) THEN
          PRINT *,' Subroutine SVFDCF.'
          PRINT *,' MHOBC(',IDCS,') = ', MHOBC( IDCS )
          STOP
        ENDIF
C
C O.k. outer b.c. is fine. Now need to
C see how many points this effects.
C
        IF ( MHOBC( IDCS ).EQ.1 ) THEN
          NSNOB = 0
        ENDIF
C
        IF ( MHOBC( IDCS ).EQ.2 .OR. MHOBC( IDCS ).EQ.3 .OR.
     1       MHOBC( IDCS ).EQ.6 .OR. MHOBC( IDCS ).EQ.7     ) THEN
          NSNOB = 1 
        ENDIF
C
        IF ( MHOBC( IDCS ).EQ.4 .OR. MHOBC( IDCS ).EQ.5 ) THEN
          NSNOB = 2 
        ENDIF
C
        DO NDER = 0, NDRVS
          ND1 = NDER + 1
          DO I = 1, NFDCM
            SVFDC( I, IRAD, ND1, IDCS ) = DZERO
          ENDDO
        ENDDO
C
C irad is the node for which we want to calculate
C our coefficients
C
        NLCS = IRAD - 1
        NRCS = NR - IRAD
C
C we wish to calculate NLN (number of left nodes)
C and NRN ( number of right nodes )
C NNDS ( total number of nodes) is then NLN + NRN + 1 ...
C
        NLN = MIN( NBN, NLCS )
        NRN = MIN( NBN, NRCS )
C
        NNDS = NLN + NRN + 1
C
C We must work out how many nodes are affected
C to the left and the right. NALF and NARF
C are respectively the number of affected nodes
C to the left and right.
C
        IF ( (IRAD-NLN).GT.NSNIB ) NALF = 0
        IF ( (IRAD-NLN).EQ.NSNIB ) NALF = 1
        IF ( (IRAD-NLN).LT.NSNIB ) NALF = 2
C
        IF ( (IRAD+NRN).LT.(NR+1-NSNOB) ) NARF = 0
        IF ( (IRAD+NRN).EQ.(NR+1-NSNOB) ) NARF = 1
        IF ( (IRAD+NRN).GT.(NR+1-NSNOB) ) NARF = 2
C
        IF ( NALF.EQ.0 .AND. NARF.EQ.0 ) THEN
          OCHNGE = .FALSE.
        ELSE
          OCHNGE = .TRUE.
        ENDIF
        IF ( .NOT. OCHNGE ) GOTO 51
C       .
C       . OK - we need to form a matrix COEFM2 such that
C       . the correct coeffcients are given when
C       . COEFM1 is multiplied by COEFM2
C       .
        L    = LARR( IDCS )
        IIBC = MHIBC( IDCS )
        IOBC = MHOBC( IDCS )
C       .
        CALL LDGNMF( NR, NNDS, NALF, NARF, L, IIBC, IOBC, NCFM,
     1             XARR, COEFM2, COEFM1, WORK1, WORK2, IWORK )
C       .
 51     CONTINUE
        X0 = XARR( IRAD )
        DO INDS = 1, NNDS
          INODE = IRAD - NLN - 1 + INDS
          WORK1( INDS ) = XARR( INODE )
        ENDDO
C
C Now ready to calculate the coefficients
C
        CALL GFDCFD( X0, WORK1, NNDS, COEFM1, NCFM, 
     1               IWORK, WORK2 )
C
C coefm matrix should now contain the coeff.s
C
        IF ( OCHNGE ) THEN
C        .
C        . Our coefficients are modified
C        . by the boundary conditions
C        .
         DO NDER = 0, NDRVS
          ND1 = NDER + 1
          DO INDS = 1, NNDS
            INODE = INDS - NLN + NBN
            FAC = EMMULT( ND1, INDS, NCFM, NCFM, NNDS,
     1                      COEFM1, COEFM2 )
            SVFDC( INODE, IRAD, ND1, IDCS ) = FAC
          ENDDO
         ENDDO
        ELSE
C        .
C        . Our coefficients are not modified
C        . by the boundary conditions
C        .
         DO NDER = 0, NDRVS
          ND1 = NDER + 1
          DO INDS = 1, NNDS
            INODE = INDS - NLN + NBN
            SVFDC( INODE, IRAD, ND1, IDCS ) = 
     1                      COEFM1( ND1, INDS )
          ENDDO
         ENDDO
C        .
        ENDIF
C
 50    CONTINUE
       ENDDO
C      .
      ENDDO
C     .
      RETURN
      END
C*********************************************************************
C*********************************************************************
C subroutine Optimum Number of Theta and Phi Points Find *************
C            -       -         -         -   -      -    *************
C Steve Gibbons 9.9.97                                               C
C____________________________________________________________________C
C For a given level of harmonics LH; this routine will return a      C
C number of theta points ( NTHPTS ) which is greater then LH         C
C and less than or equal to NTHMAX. Also a number of PHI points      C
C ( NPHPTS ) which is greater than 2*MMAX and is also a              C
C power of 2 and is also smaller than NPHMAX. Failiure to do either  C
C will be reported.                                                  C
C____________________________________________________________________C
C Input variables :-                                                 C
C ===============                                                    C
C  Integer                                                           C
C  -------                                                           C
C     LH        : Highest degree, l, of spherical harmonic.          C
C     MMAX      : Highest order, m, of spherical harmonic.           C
C                                                                    C
C     For a fully 3-D problem, MMAX will equal LH but may be         C
C     less if certain symmetries are applied in the azimuthal        C
C     direction.                                                     C
C                                                                    C
C     NTHPTS    : Number of theta points.                            C
C     NPHPTS    : Number of phi points.                              C
C     NTHMAX    : Maximum number of theta points.                    C
C     NPHMAX    : Maximum number of phi points.                      C
C                                                                    C
C____________________________________________________________________C
C
C*********************************************************************
      SUBROUTINE ONTPPF ( LH, MMAX, NTHPTS, NPHPTS, NTHMAX, NPHMAX )
      IMPLICIT NONE
C____________________________________________________________________C
C Variable declarations - Parameters ................................C
      INTEGER LH, MMAX, NTHPTS, NPHPTS, NTHMAX, NPHMAX
C____________________________________________________________________C
C START OF PROGRAM **************************************************C
C____________________________________________________________________C
C
      IF ( LH.EQ.0 ) THEN
        PRINT *,' Subroutine ONTPPF. LH = 0.'
        STOP
      ENDIF
C
      IF ( MMAX.LT.0 .OR. MMAX.GT.LH ) THEN
        PRINT *,' Subroutine ONTPPF. MMAX = ', MMAX
        PRINT *,' Must be between 0 and LH (= ',LH,')'
        STOP
      ENDIF
C
      NPHPTS = 2
 500  CONTINUE
      IF ( 2*MMAX.GE.NPHPTS ) THEN
         NPHPTS = NPHPTS*2
         GOTO 500
      ENDIF
      IF ( NPHPTS.GT.NPHMAX ) THEN
         PRINT *,' Subroutine ONTPPF.'
         PRINT *,' NPHPTS must be atleast ', NPHPTS
         PRINT *,' NPHMAX = ', NPHMAX
         PRINT *,' Program aborted.'
         STOP
      ENDIF
C
      NTHPTS = LH+1
      IF ( NTHPTS/2*2.NE.NTHPTS ) NTHPTS = NTHPTS + 1
      IF ( NTHPTS.GT.NTHMAX ) THEN
         PRINT *,' Subroutine ONTPPF.'
         PRINT *,' NTHPTS must be atleast ', LH+1
         PRINT *,' NTHMAX = ', NTHMAX
         PRINT *,' Program aborted.'
         STOP
      ENDIF
C
      RETURN
      END
C*********************************************************************
C*********************************************************************
C subroutine GAUWTS **************************************************
C Adapted 22.4.97 from Numerical Recipes routine GAULEG              C
C____________________________________________________________________C
C Input variables :-                                                 C
C ===============                                                    C
C  Integer                                                           C
C  -------                                                           C
C     NTHPTS	: Number of theta points.                            C
C                                                                    C
C  Double Precision                                                  C
C  ----------------                                                  C
C     X1	: Starting value for integration.                    C
C     X2	: Ending value for integration.                      C
C____________________________________________________________________C
C Output variables :-                                                C
C ================                                                   C
C  Double Precision                                                  C
C  ----------------                                                  C
C     GAUX	: Array containing abscissae of the Gauss-Legendre   C
C                  NTHPTS-points quadrature formula.                 C
C     GAUW      : Array containing the weights for the above points. C
C                                                                    C
C ( Both GAUX and GAUW have dimension NTHPTS ).                      C
C____________________________________________________________________C
C
C*********************************************************************
      SUBROUTINE GAUWTS ( X1, X2, GAUX, GAUW, NTHPTS )
      IMPLICIT NONE
C____________________________________________________________________C
C Variable declarations - Parameters ................................C
      INTEGER NTHPTS
      DOUBLE PRECISION X1, X2, GAUX( NTHPTS ), GAUW( NTHPTS )
C____________________________________________________________________C
C Variable declarations - Working variables .........................C
      INTEGER M,I,J
      DOUBLE PRECISION XM,XL,P1,P2,P3,EPS,PP,Z,Z1
      PARAMETER (EPS=1.0d-13)
      DOUBLE PRECISION PI
      PARAMETER (PI=3.14159265358979312D0)
C____________________________________________________________________C
C START OF PROGRAM **************************************************C
C____________________________________________________________________C
C 
      IF ( ABS( X1 - X2 ).LT.EPS ) THEN
        PRINT *,' Subroutine GAUWTS,'
        PRINT *,' X1 = ', X1
        PRINT *,' X2 = ', X2
        PRINT *,' Program aborted.'
        STOP
      ENDIF
C ....... the roots are symmetric in the interval so need only find
C half of them.
      M = ( NTHPTS + 1 )/2
      XM = 0.5d0 * ( X2 + X1 )
      XL = 0.5d0 * ( X2 - X1 )
C ........start looping over the desired roots .......................
      DO I = 1, M
         Z = DCOS( PI*( DBLE(I) - 0.25d0)/( DBLE(NTHPTS) + 0.5d0 ))
C           ..... starting with this approximation to the Ith root, we
C                enter the main loop of refinement by Newton's method.
 100     CONTINUE
            P1 = 1.0D0
            P2 = 0.0D0
C           ........... Loop up the recurrence relation to get the
C                      legendre Polynomial evaluated at Z.
            DO J = 1, NTHPTS
               P3 = P2
               P2 = P1
               P1 = ((2.0d0*J-1.0d0)*Z*P2 - (J-1.0d0)*P3)/DBLE( J )
            ENDDO
C           ..................... finish recurrence relation loop ...
C ... P1 is now the desired Legendre Polynomial. We now compute PP,
C    its derivative by a standard relation involving also P2, the 
C    polynomial of one order lower.
            PP = NTHPTS*(Z*P1-P2)/(Z*Z-1.0d0)
            Z1 = Z
            Z = Z1 - P1/PP
         IF ( ABS(Z-Z1).GT.EPS ) GOTO 100
C ...........scale the root to the desired interval .................
         GAUX( I ) = XM - XL*Z
C ...........and add its symmetric counterpart ......................
         GAUX( NTHPTS+1-I ) = XM + XL*Z
C ...........calculate the weight ...................................
         GAUW( I ) = 2.0d0*XL/((1.0d0-Z*Z)*PP*PP)
C ...........and add its symmetric counterpart ......................
         GAUW( NTHPTS + 1 - I ) = GAUW( I )
      ENDDO
C ......... end looping over the desired roots .......................
      RETURN
      END
C*********************************************************************

C*********************************************************************
C subroutine SCHmidt Normalised Legendre function Array **************
C            ---     -          -                 -     **************
C Steve Gibbons 22.4.97                                              C
C____________________________________________________________________C
C Does the same as SCHNLF except that instead of a single valued X   C
C for one theta point, it fills arrays PA and DPA with the           C
C legendre Functions etc. for each of the NTHPTS values of cos(theta)C
C in the array GAUX.                                                 C
C____________________________________________________________________C
C Input variables :-                                                 C
C ===============                                                    C
C  Integer                                                           C
C  -------                                                           C
C     LH	: Highest degree, l, of spherical harmonic.          C
C     NTHPTS	: Number of theta points.                            C
C  Double Precision                                                  C
C  ----------------                                                  C
C     PA	: Schmidt Normalised Legendre Functions Dimension.   C
C		   {  ( LH + 1 )*( LH + 2 )/2 , NTHPTS }             C
C		   P_l^m(X_j) = PA( L*(L+1)/2 + M + 1 , j )          C
C     DPA	: Derivatives of the above.                          C
C     GAUX	: Array of cosines to the NTHPTS angles.             C
C                  Dimension ( NTHPTS ).                             C
C____________________________________________________________________C
C Functions ... Calling Proceedures :-                               C
C  Double Precision                                                  C
C  ----------------                                                  C
C PMM ( M, S )					                     C
C DPMM ( M, C, S )				                     C
C PMM1 ( M, X, PMM0 )                                                C
C PLM ( L, M, X, PLMIN1, PLMIN2 )				     C
C DPMM1 ( M , X , S, PMM , DPMM)                                     C
C DPLM ( L, M , X , S, PMM1 , DPMM1, DPMM2 )                         C
C____________________________________________________________________C
C
C*********************************************************************
      SUBROUTINE SCHNLA ( PA, DPA, GAUX, LH, NTHPTS )
      IMPLICIT NONE
C____________________________________________________________________C
C Variable declarations - Parameters ................................C
      INTEGER LH, NTHPTS
      DOUBLE PRECISION GAUX( NTHPTS ),
     1                 PA ( ( LH + 1 )*( LH + 2 )/2 , NTHPTS),
     2                 DPA ( ( LH + 1 )*( LH + 2 )/2 , NTHPTS)
C____________________________________________________________________C
C Variable declarations - Working variables .........................C
      INTEGER INDEX,L,M,IOLD1,IOLD2,NTHETA
      DOUBLE PRECISION SINE,PMIN1,PMIN2,TOL,DPMIN1,
     1                 DPMIN2,X
      PARAMETER (TOL=1.0d-6)
C____________________________________________________________________C
C Variable declarations - Functions called ..........................C
      DOUBLE PRECISION PMM,PMM1,PLM,DPMM,DPMM1,DPLM
C____________________________________________________________________C
C START OF PROGRAM **************************************************C
C____________________________________________________________________C
C Check validity of arguments ....
C
C ........................... now loop around theta points
      DO NTHETA = 1, NTHPTS
C
        X = GAUX ( NTHETA )
        SINE = X*X
        IF ( SINE.GT.1.0d0 ) THEN
           PRINT *,' Subroutine SCHNLA.'
           PRINT *,' Illegal Cos(theta) has been entered.'
           PRINT *,' ( For NTHETA = ',NTHETA,' )'
           PRINT *,' Program stopped.'
           STOP
        ENDIF
C
C Set SINE (theta) in terms of X
        SINE = DSQRT ( (1.0d0 + X)*(1.0d0 - X) )
        IF ( SINE.LT.TOL ) THEN
           PRINT *,' Subroutine SCHNLA.'
           PRINT *,' SINE is too small. Division by zero imminent.'
           PRINT *,' ( For NTHETA = ',NTHETA,' )'
           PRINT *,' Program stopped.'
           STOP
        ENDIF
C..................... first calculate the P_l^m ..........
        DO M = 0, LH - 2
C                        ............. Calculate P_M^M .....
           L = M
           INDEX = L*(L+1)/2+M+1
           PA ( INDEX , NTHETA) = PMM ( M , SINE )
           DPA ( INDEX , NTHETA) = DPMM ( M , X, SINE )
C                         ............. Calculate P_(M+1)^M .
           PMIN1 = PA ( INDEX , NTHETA)
           DPMIN1 = DPA ( INDEX , NTHETA)
           IOLD2 = INDEX
           L = L + 1
           INDEX = L*(L+1)/2+M+1
           PA (INDEX , NTHETA) = PMM1 ( M , X , PMIN1 )
           DPA (INDEX , NTHETA) = DPMM1 (M,X,SINE,PMIN1, DPMIN1)
           IOLD1 = INDEX
C                         ......... Calculate P_L^M general .
           DO L = M + 2, LH
              PMIN2 = PA ( IOLD2 , NTHETA)
              PMIN1 = PA ( IOLD1 , NTHETA)
              DPMIN2 = DPA ( IOLD2 , NTHETA)
              DPMIN1 = DPA ( IOLD1 , NTHETA)
              INDEX = L*(L+1)/2+M+1
              PA ( INDEX , NTHETA) = PLM ( L,M,X,PMIN1,PMIN2 )
              DPA ( INDEX , NTHETA) = DPLM (L,M,X,SINE , PMIN1, 
     1                              DPMIN1, DPMIN2 )
              IOLD2 = IOLD1
              IOLD1 = INDEX
           ENDDO
        ENDDO
        M = LH - 1
        L = M
        INDEX = L*(L+1)/2+M+1
        PA( INDEX , NTHETA) = PMM ( M , SINE )
        DPA ( INDEX , NTHETA) = DPMM ( M , X, SINE )
        PMIN1 = PA( INDEX , NTHETA)
        DPMIN1 = DPA( INDEX , NTHETA)
        L = LH
        INDEX = L*(L+1)/2+M+1
        PA( INDEX , NTHETA) = PMM1 ( M , X , PMIN1 )
        DPA(INDEX ,NTHETA) = DPMM1 (M ,X ,SINE,PMIN1,DPMIN1 )
        M = LH
        INDEX = L*(L+1)/2+M+1
        PA( INDEX , NTHETA) = PMM ( M , SINE )
        DPA ( INDEX , NTHETA) = DPMM ( M , X, SINE )
C......................finished calculating P_l^m .........
      ENDDO
C......................finished looping around theta points

      RETURN
      END
C*********************************************************************
C*********************************************************************
C subroutine Complete Adapted Solution Vector Derivative 2 ***********
C            -        -       -        -      -          - ***********
C Steve Gibbons Sat Feb  5 08:34:24 GMT 2000                         C
C____________________________________________________________________C
C                                                                    C
C Loops around the entire convection vector and calculates           C
C derivatives 0 and up to 2 of every harmonic from grid node         C
C IR = ILN to grid node IRN.                                         C
C                                                                    C
C____________________________________________________________________C
C                                                                    C
C Input variables :-                                                 C
C ===============                                                    C
C  Integer                                                           C
C  -------                                                           C
C                                                                    C
C     ILN       : Left-most node to be acted upon.                   C
C     IRN       : Right-most node to be acted upon.                  C
C     NBN       : Number of bounding nodes. See above.               C
C     IHD       : Highest derivative requested.                      C
C     NCFM      : Leading dimension of SVFDC. At least (2*NBN+1)     C
C     NR        : Number of radial grid nodes in each function.      C
C     NDRVS     : Number of highest derivative for which             C
C                  coefficients are stored by the array SVFDC.       C
C     NDRVM     : Limit on NDRVS. Array bound for SVFDC.             C
C     INARR     : Integer array dimension ( * ).                     C
C                  Elements may be arbitrary except for              C
C                  INARR( 1 ) = IFORMF - flag for vector format.     C
C                                        See INDFUN                  C
C                  INARR( 2 ) = NRR. Must be consistent with NR.     C
C                  INARR( 3 ) = NH = total number of radial func.s   C
C     NDCS      : Maximum distinct finite difference schemes         C
C                  stored by the array SVFDC.                        C
C     MHP       : Pointer array for harmonics. If HMP( ih ) = is     C
C                  then 'is' is the finite difference scheme used    C
C                   to take derivatives of that harm. radial func.   C
C                   If MHP is negative, the harmonic is avoided and  C
C                   CASVD2 moves on to the next harmonic.            C
C                                                                    C
C  Double Precision                                                  C
C  ----------------                                                  C
C                                                                    C
C     V         : Solution vector. Dim ( * ) but length atleast      C
C                  NR*NH.                                            C
C                                                                    C
C     SVFDC     : Finite difference coefficient matrix.              C
C                  Dimension ( NCFM, NR, NDRVM+1, NDCS ).            C
C                   Array is generated by the routine svfdcf         C
C                 See documentation for SVFDCF for details.          C
C                                                                    C
C     V0        : Zero^th derivatives. Dim ( * ) but length atleast  C
C     V1        : First   derivatives. Dim ( * ) but length atleast  C
C     V2        : Second  derivatives. Dim ( * ) but length atleast  C
C                                                                    C
C IMPORTANT: ALL of V0, V1 and V2 are filled even if                 C
C these derivatives are not requested!                               C
C                                                                    C
C____________________________________________________________________C
C
C*********************************************************************
      SUBROUTINE CASVD2 ( V, ILN, IRN, NBN, IHD, NCFM, NR, NDRVS,
     1                   NDRVM, INARR, NDCS, MHP, SVFDC, V0, V1, V2 )
      IMPLICIT NONE
C____________________________________________________________________C
C Variable declarations - Parameters ................................C
      INTEGER NBN, IHD, NCFM, NR, NDRVS, NDRVM, INARR( * ), NDCS,
     1        MHP( * ), ILN, IRN
      DOUBLE PRECISION V( * ), V0( * ), V1( * ), V2( * ),
     1                 SVFDC( NCFM, NR, NDRVM+1, NDCS )
C____________________________________________________________________C
C Variable declarations - Working variables .........................C
C
      DOUBLE PRECISION DERV( 3 )
      INTEGER IR, IH, NH, IND, INDFUN, IS, I
C____________________________________________________________________C
C START OF PROGRAM **************************************************C
C____________________________________________________________________C
C
      DO I = 1, 3
        DERV( I ) = 0.0d0
      ENDDO
C
      IF ( IHD.LT.0 .OR. IHD.GT.2 ) THEN
        PRINT *,' Subroutine CASVD2. IHD = ', IHD
        PRINT *,' Program aborted.'
        STOP
      ENDIF
C
      NH = INARR( 3 )
      DO IH = 1, NH
        IS  = MHP( IH )
        IF ( IS.LT.1 ) GOTO 60
        DO IR = ILN, IRN
C         .
          IND = INDFUN( IR, IH, INARR )
          CALL ASVDR ( V, IR, IS, IH, NBN, IHD, NCFM, NR, NDRVS,
     1                   NDRVM, DERV, INARR, SVFDC, NDCS )
C         .
          V0( IND ) = DERV( 1 )
          V1( IND ) = DERV( 2 )
          V2( IND ) = DERV( 3 )
C         .
        ENDDO
 60   CONTINUE
      ENDDO
C     .
      RETURN
      END
C*********************************************************************
C*********************************************************************
C subroutine Single Harmonic Kinetic Energy Evaluation Routine *******
C            -      -        -       -      -          -       *******
C Steve Gibbons Thu Oct 28 08:51:22 GMT 1999                         C
C____________________________________________________________________C
C                                                                    C
C Returns 0.5d0 * \int_{volume} A.A dV where A is a single vector    C
C harmonic - the type being indicated by MHT( ih ).                  C
C                                                                    C
C This can be either magnetic energy or kinetic energy.              C
C                                                                    C
C The value for A = ( A_r, A_theta, A_phi ) is returned in DKE( 1 ). C
C                                                                    C
C The value for A = ( A_r, 0 , 0 ) is returned in DKE( 2 ).          C
C                                                                    C
C If MHT( ih ) does not correspond to a poloidal or toroidal         C
C vector harmonic, DKE( 1 ) and DKE( 2 ) are both returned zero.     C
C                                                                    C
C____________________________________________________________________C
C                                                                    C
C Input variables :-                                                 C
C ===============                                                    C
C  Integer                                                           C
C  -------                                                           C
C                                                                    C
C     IH        : Number of harmonic to be evaluated.                C
C                                                                    C
C     NDCS       : Number of distinct differencing coeff.s           C
C                  represented in SVFDC.                             C
C                                                                    C
C     NR        : Number of radial grid nodes.                       C
C                                                                    C
C     INARR     : Array of integers. Dimension ( * )                 C
C                 At present the elements of INARR are as follows.   C
C                                                                    C
C                 INARR( 1 ) = IFORMF. Format flag.                  C
C                  This flag decides how the elements are arranged   C
C                  in the solution vector.                           C
C                                                                    C
C                  The current options are:-                         C
C                                                                    C
C                   IFORMF = 1, 3. INDFUN = ( IR - 1 )*NH + IH       C
C                   IFORMF = 2, 4. INDFUN = ( IH - 1 )*NR + IR       C
C                                                                    C
C  where IR and IH are the current grid node and harmonic resp.      C
C  and NR and NH are the total numbers of nodes and harmonics        C
C  in the solution vector.                                           C
C                                                                    C
C                 INARR( 2 ) = NR. Number of radial grid nodes.      C
C                 INARR( 3 ) = NH. Number of harmonics in sol. vect. C
C                                                                    C
C     MHT       : Array length ( * ) - atleast length NH             C
C                                                                    C
C MHT defines what each scalar function in a solution vector         C
C represents.                                                        C
C                                                                    C
C         MHT( IH ) = 1 --> harmonic is poloidal velocity            C
C         MHT( IH ) = 2 --> harmonic is toroidal velocity            C
C         MHT( IH ) = 3 --> harmonic is temperature.                 C
C         MHT( IH ) = 4 --> harmonic is poloidal magnetic field      C
C         MHT( IH ) = 5 --> harmonic is toroidal magnetic field      C
C                                                                    C
C     MHL       : Array length ( * ) - atleast length NH             C
C                  Sph. harm. degree, l.                             C
C                                                                    C
C     MHP       : Array length ( * ) - atleast length NH             C
C                  Pointer array to finite difference coefficients.  C
C                  MHP( ih ) = is, which is the 4th index of         C
C                  array SVFDC - indicates f.d. scheme used.         C
C                                                                    C
C     NBN       : Number of nodes on each side of point for          C
C                  central differences.                              C
C                                                                    C
C     NDRVS     : Highest derivative stored in SVFDC.                C
C                (Not needed if we are doing tor --> pol but         C
C                 must be atleast 2 if doing pol --> tor )           C
C     NDRVM     : Limit on NDRVS. Array bound for SVFDC.             C
C                                                                    C
C     NFDCM     : Leading dim of SVFDC. See SVFDCF.                  C
C                  (Must be atleast 2*NBN + 1 )                      C
C                                                                    C
C  Double Precision                                                  C
C  ----------------                                                  C
C                                                                    C
C     SV        : Solution vector. Dim ( * ) ( input )               C
C                 Length must be atleast NR*NH                       C
C                                                                    C
C     XARR      : Array of dimension ( NR )                          C
C                 XARR( j ) = element x_j                            C
C                                                                    C
C     SVFDC     : Finite difference coefficient matrix.              C
C                  Dimension ( NFDCM, NR, NDRVM+1, NDCS ).           C
C                   Array is generated by the routine svfdcf         C
C                 See documentation for SVFDCF for details.          C
C                                                                    C
C Output variables :-                                                C
C ================                                                   C
C                                                                    C
C  Double Precision                                                  C
C  ----------------                                                  C
C                                                                    C
C     DKE       : Dimension ( 2 ).                                   C
C                                                                    C
C                 DKE( 1 ) is returned with the whole contribution   C
C                 to the kinetic energy from harmonic IH             C
C                                                                    C
C                 DKE( 2 ) is returned with kinetic energy of the    C
C                 radial component of harmonic IH                    C
C                                                                    C
C                                                                    C
C____________________________________________________________________C
C
C*********************************************************************
      SUBROUTINE SHKEER( IH, NDCS, NR, INARR, MHT, MHL, MHP, NBN,
     1                   NDRVS, NDRVM, NFDCM, SV, XARR, DKE, SVFDC )
      IMPLICIT NONE
C____________________________________________________________________C
C Variable declarations - Parameters ................................C
      INTEGER IH, NDCS, NR, INARR( * ), MHT( * ), MHL( * ), MHP( * ),
     1        NBN, NDRVS, NDRVM, NFDCM
      DOUBLE PRECISION XARR( NR ), SV( * ), DKE( 2 ),
     1                 SVFDC( NFDCM, NR, NDRVM+1, NDCS )
C____________________________________________________________________C
C Variable declarations - Working variables .........................C
C
      INTEGER IRAD, IHD, ITYPE, L, NR2, IS
      DOUBLE PRECISION COEF, DK, SQRLL1, PI, DERV( 2 ), D0F, D1F,
     1                 RAD, ER, ETOT, FAC
      PARAMETER (PI=3.14159265358979312D0)
C____________________________________________________________________C
C START OF PROGRAM **************************************************C
C____________________________________________________________________C
C
      NR2 = INARR( 2 )
      IF ( NR2.NE.NR ) THEN
        PRINT *,' Subroutine SHKEER.'
        PRINT *,' NR         = ', NR
        PRINT *,' INARR( 2 ) = ', NR2
        PRINT *,' Program aborted.'
        STOP
      ENDIF
C
      DKE( 1 ) = 0.0d0
      DKE( 2 ) = 0.0d0
C     .
      ITYPE = MHT( IH )
      L     = MHL( IH )
      IS    = MHP( IH )
      DK    = SQRLL1( L )
C     .
C     . First do poloidal harmonics
C     .
      IF ( ITYPE.EQ.1 .OR. ITYPE.EQ.4 ) THEN
        IHD = 1
        DO IRAD = 1, NR
          RAD = XARR( IRAD )
C         .
          IF ( IRAD.EQ.1 ) THEN
            COEF = 0.5d0*(  XARR( 2 ) - XARR( 1 )  )
          ENDIF
C         .
          IF ( IRAD.GT.1 .AND. IRAD.LT.NR ) THEN
            COEF = 0.5d0*(  XARR( IRAD+1 ) - XARR( IRAD-1 )  )
          ENDIF
C         .
          IF ( IRAD.EQ.NR ) THEN
            COEF = 0.5d0*(  XARR( NR ) - XARR( NR-1 )  )
          ENDIF
C         .
          CALL ASVDR( SV, IRAD, IS, IH, NBN, IHD, NFDCM, NR, NDRVS,
     1                NDRVM, DERV, INARR, SVFDC, NDCS )
          D0F = DERV( 1 )
          D1F = DERV( 2 )
C         .
C         . Add contribution to ER and ETOT
C         .
          ER   = D0F*D0F
          ETOT = DK*DK*D0F*D0F + (D0F + RAD*D1F)*(D0F + RAD*D1F)
C         .
C         . Add to cumulative integral
C         .
          DKE( 1 ) = DKE( 1 ) + ETOT*COEF
          DKE( 2 ) = DKE( 2 ) + ER*COEF
C         .
        ENDDO
C       .
C       . Finally, multiply by leading factors
C       .
        FAC = 2.0d0*PI*DK*DK/(2.0d0*DBLE( L ) + 1.0d0)
        DKE( 1 ) = DKE( 1 )*FAC
        DKE( 2 ) = DKE( 2 )*FAC*DK*DK
C       .
      ENDIF
C     .
C     . Now do toroidal harmonics
C     .
      IF ( ITYPE.EQ.2 .OR. ITYPE.EQ.5 ) THEN
        IHD = 0
        DO IRAD = 1, NR
          RAD = XARR( IRAD )
C         .
          IF ( IRAD.EQ.1 ) THEN
            COEF = 0.5d0*(  XARR( 2 ) - XARR( 1 )  )
          ENDIF
C         .
          IF ( IRAD.GT.1 .AND. IRAD.LT.NR ) THEN
            COEF = 0.5d0*(  XARR( IRAD+1 ) - XARR( IRAD-1 )  )
          ENDIF
C         .
          IF ( IRAD.EQ.NR ) THEN
            COEF = 0.5d0*(  XARR( NR ) - XARR( NR-1 )  )
          ENDIF
C         .
          CALL ASVDR( SV, IRAD, IS, IH, NBN, IHD, NFDCM, NR, NDRVS,
     1                NDRVM, DERV, INARR, SVFDC, NDCS )
          D0F = DERV( 1 )
C         .
C         . Add contribution to ER and ETOT
C         .
          ETOT = RAD*RAD*D0F*D0F
C         .
C         . Add to cumulative integral
C         .
          DKE( 1 ) = DKE( 1 ) + ETOT*COEF
C         .
        ENDDO
C       .
C       . Finally, multiply by leading factors
C       .
        FAC = 2.0d0*PI*DK*DK/(2.0d0*DBLE( L ) + 1.0d0)
        DKE( 1 ) = DKE( 1 )*FAC
      ENDIF
C     .
      RETURN
      END
C*********************************************************************
C*********************************************************************
C subroutine Solution Vector 2 Xtra Special Array ********************
C            -        -      - -    -       -     ********************
C Steve Gibbons Sun Apr  8 11:43:40 MET DST 2001                     C
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
      SUBROUTINE SV2XSA( NCMX, NPHP, NTHP, NR, ILN, IRN, LH, INARR, 
     1                   MHT, MHL, MHM, ICM, ICOMP, GAUX, PA, DPA,
     2                   FTF, XARR, V0, V1, XSA )
      IMPLICIT NONE
C____________________________________________________________________C
C Variable declarations - Parameters ................................C
      INTEGER          NCMX, NPHP, NTHP, NR, ILN, IRN, LH, INARR( * ),
     1                 MHT( * ), MHL( * ), MHM( * ), ICM, ICOMP
      DOUBLE PRECISION GAUX( NTHP ), FTF( 2*NPHP ),
     1                 PA( ( LH + 1 )*( LH + 2 )/2, NTHP ),
     2                 DPA( ( LH + 1 )*( LH + 2 )/2, NTHP )
      DOUBLE PRECISION XARR( NR ), V0( * ), V1( * ),
     1                 XSA( NCMX, NPHP, NTHP, NR )
C____________________________________________________________________C
C Variable declarations - Working variables .........................C
      INTEGER          NRR, NH, IR, ITHP, IPHP, IOP, LENV, ISIGN,
     1                 L, M, ICS, IP, IND, INDSAM,
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
        PRINT *,' Subroutine SV2XSA.'
        PRINT *,' NR  = ', NR
        PRINT *,' NRR = ', NRR
        PRINT *,' Program stopped.'
        STOP
      ENDIF
C
      IF ( ICM.LT.1 .OR. ICM.GT.NCMX ) THEN
        PRINT *,' Subroutine SV2XSA.'
        PRINT *,' ICM = ', ICM,' NCMX = ', NCMX
        PRINT *,' Program aborted.'
        STOP
      ENDIF
C
      IF ( ICOMP.LT.1 .OR. ICOMP.GT.15 ) THEN
        PRINT *,' Subroutine SV2XSA.'
        PRINT *,' ICOMP = ', ICOMP
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
            INDSAM = 2*M + 1
            INDDIF = 2*M + 2
            DPTERM = (-1.0d0)*M*PA( IP, ITHP )/( SINE*SQRLL1( L ))
          ENDIF
          IF ( ICS.EQ.2 ) THEN
            INDDIF = 2*M + 1
            INDSAM = 2*M + 2
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

C*********************************************************************
C subroutine Xtra Special Vector Single Component Integrate **********
C            -    -       -      -      -         -         **********
C Steve Gibbons Mon Apr  9 12:18:19 MET DST 2001                     C
C____________________________________________________________________C
C                                                                    C
C The function XSV( NCMX, NPHP, NTHP, NR ) defines, in real space,   C
C the values of a function. SF( ICMP, iphi, ithe, ir ) is the value  C
C at  RAD = xarr( ir )                                               C
C                                                                    C
C  THETA = ACOS[ GAUX( ithe ) ]                                      C
C                                 and                                C
C  PHI   = (iphi-1)*DELTAP                                           C
C                                 with                               C
C  deltap = 2*pi/(NPHP*M0).                                          C
C                                                                    C
C XSVSCI will return the value DINT, the volume integral (over the   C
C full spherical shell or sphere) of the function stored in the      C
C component ICMP.                                                    C
C                                                                    C
C____________________________________________________________________C
C                                                                    C
C Input variables :-                                                 C
C ===============                                                    C
C                                                                    C
C  Integer                                                           C
C  -------                                                           C
C                                                                    C
C     NCMX      : Maximum number of components stored in XSV         C
C     NPHP      : Number of phi points.                              C
C     NTHP      : Number of theta points.                            C
C     NR        : Total number of radial grid nodes.                 C
C     ICM       : Component of XSV which stores scalar function.     C
C                                                                    C
C  Double Precision                                                  C
C  ----------------                                                  C
C                                                                    C
C     XARR      : Array of radial spacing values. Dim ( NR )         C
C                                                                    C
C     GAUW      : Gauss weights as calculated by GAUWTS. ( NTHP )    C
C                                                                    C
C     XSV       : Xtra Special Function: is an array containing a    C
C                  function over a set of theta, phi and r points.   C
C                    Dimensions are                                  C
C                      ( NCMX, NPHP, NTHP, NR )                      C
C                                                                    C
C     DINT      : Spherical integral of function.                    C
C                                                                    C
C____________________________________________________________________C
C
C*********************************************************************
      SUBROUTINE XSVSCI( NCMX, NPHP, NTHP, NR, ICM, XARR, GAUW,
     1                   XSV, DINT )
      IMPLICIT NONE
C____________________________________________________________________C
C Variable declarations - Parameters ................................C
      INTEGER          NCMX, NPHP, NTHP, NR, ICM
      DOUBLE PRECISION GAUW( NTHP ), DINT, XARR( NR ),
     1                 XSV( NCMX, NPHP, NTHP, NR )
C____________________________________________________________________C
C Variable declarations - Working variables .........................C
      INTEGER          ITHE, IPHI, IR
      DOUBLE PRECISION PI, ZERO, ZCOEF, WEIGHT, RAD, FAC1, FAC2,
     1                 VFAC, DRAD
      PARAMETER        ( ZERO = 0.0d0, PI = 3.14159265358979312D0 )
C____________________________________________________________________C
C START OF PROGRAM **************************************************C
C____________________________________________________________________C
C
C Check value of ICM
C
      IF ( ICM.LT.1 .OR. ICM.GT.NCMX ) THEN
        PRINT *,' Subroutine XSVSCI.'
        PRINT *,' ICM = ', ICM,' NCMX = ', NCMX
        PRINT *,' Program aborted.'
        STOP
      ENDIF
C
      DINT = ZERO
C
      FAC1 = 0.5d0/DBLE( NPHP )
C
C ............... begin looping around radial grid nodes .............
      DO IR = 1, NR
        RAD  = XARR( IR )
        VFAC = 4.0d0*PI*RAD**2
        FAC2 = ZERO
C
C ............... begin looping around theta points ..................
        DO ITHE = 1, NTHP
          WEIGHT = GAUW( ITHE )*FAC1
          ZCOEF  = ZERO
          DO IPHI = 1, NPHP
            ZCOEF = ZCOEF + XSV( ICM, IPHI, ITHE, IR )
          ENDDO
          FAC2 = FAC2 + WEIGHT*ZCOEF
        ENDDO
C       .
C       . End looping around theta points
C       .
        IF ( IR.EQ.1 ) THEN
          DRAD = 0.5d0*(  XARR( 2 ) - XARR( 1 )  )
        ENDIF
C       .
        IF ( IR.GT.1 .AND. IR.LT.NR ) THEN
          DRAD = 0.5d0*(  XARR( IR+1 ) - XARR( IR-1 )  )
        ENDIF
C       .
        IF ( IR.EQ.NR ) THEN
          DRAD = 0.5d0*(  XARR( NR ) - XARR( NR-1 )  )
        ENDIF
C       .
        DINT = DINT + DRAD*FAC2*VFAC
C       .
      ENDDO
C     .
C     . End looping around radial grid nodes
C     .
      RETURN
      END
C*********************************************************************
C********************************************************************
C subroutine File CLOSE *********************************************
C            -    ----- *********************************************
C Steve Gibbons 14.4.97                                             C
C  ( note that this is essentially the routine of Dan Gordon        C
C___________________________________________________________________C
C Closes file with integer logical unit LU, filename FNAME.         C
C LABEL contains any other information regarding the nature of the  C
C the file.							    C
C___________________________________________________________________C
C Input Variables :-						    C
C ===============						    C
C  Integer							    C
C  -------							    C
C     LU	: Number of file.				    C
C  Character							    C
C  ---------							    C
C     FNAME	: Name of file. Undefined length		    C
C     LABEL	: Any further information. Undefined length         C
C___________________________________________________________________C
C
C********************************************************************
      SUBROUTINE FCLOSE ( LU, FNAME, LABEL )
      IMPLICIT NONE
C___________________________________________________________________C
C Variable Declarations - Parameters ...............................C
      INTEGER LU
      CHARACTER *(*) FNAME
      CHARACTER *(*) LABEL
C___________________________________________________________________C
C START OF PROGRAM *************************************************C
C___________________________________________________________________C
      IF ( LU.EQ.0 ) THEN
         PRINT *,' Subroutine FCLOSE '
         PRINT *,' I bet you ve forgotten to set LU ...??'
         PRINT *,' Think again and come back when you have'
         PRINT *,' remembered that LU must be a non zero integer!'
         PRINT *,' See you later. Bye for now!!'
         STOP
      ENDIF
C----------------------------
      CLOSE (UNIT=LU, STATUS='KEEP', ERR=989 )
      RETURN
C
 989  PRINT *,' Error.  Failed to close ', LABEL, ' file ', FNAME
      STOP
      END
C********************************************************************
C*********************************************************************
C subroutine File NAME giveR *****************************************
C            -    ----     - *****************************************
C Steve Gibbons 14.4.97 (Adapted from Dan Gordon's course.)          C
C____________________________________________________________________C
C Asks the user for a file name FNAME. Pretty simple really ...      C
C____________________________________________________________________C
C Input Variable :-						     C
C ==============   						     C
C  Character							     C
C  ---------							     C
C     LABEL	: Message arbitrary length  			     C
C____________________________________________________________________C
C Output Variable :-						     C
C ===============   						     C
C  Character							     C
C  ---------							     C
C     FNAME	: Filename arbitrary length			     C
C____________________________________________________________________C
C
C*********************************************************************
      SUBROUTINE FNAMER ( FNAME, LABEL )
      IMPLICIT NONE
      CHARACTER *(*) FNAME
      CHARACTER *(*) LABEL

      PRINT *,' Please enter a filename for FNAME.'
      PRINT *, LABEL
      READ (5, 200) FNAME
 200  FORMAT (A)
      RETURN
      END
C*********************************************************************
C*********************************************************************
C subroutine Linear Dependence of Grid Node Matrix Form **************
C            -      -             -    -    -      -    **************
C Steve Gibbons Sat Oct 23 15:01:52 BST 1999                         C
C____________________________________________________________________C
C                                                                    C
C Let f be a function of x and let f_j denote the value of f at      C
C the grid node j (x value is x_j given by XARR( j ) ... )           C
C                                                                    C
C If f has to satisfy a particular boundary condition then           C
C all the f_j may not be linearly independent.                       C
C                                                                    C
C For a group of n ( = NNDS ) grid nodes, from s = k+1 to k+n        C
C for some integer, k,                                               C
C                                                                    C
C  f_{k+i} = \sum_{s=1}^n DMAT( i, s ) f_{k+s},                      C
C                                                                    C
C and the routine LDGNMF returns the matrix DMAT.                    C
C                                                                    C
C The number of linearly dependent nodes to the left and right are   C
C given by NALF and NARF respectively.                               C
C                                                                    C
C The boundary conditions at the inner and outer boundaries are      C
C specified by the integers IIBC and IOBC which may take the         C
C following values                                                   C
C                                                                    C
C    IIBC              Inner Boundary Condition                      C
C    ====              ========================                      C
C                                                                    C
C      1       No condition to be applied                            C
C      2       Function must vanish                                  C
C      3       First derivative must vanish                          C
C      4       Function AND first derivative must vanish             C
C      5       Function AND second derivative must vanish            C
C      6       rdf/dr - f(r) = 0                                     C
C      7       r df/dr - l f(r) = 0                                  C
C                                                                    C
C    IOBC              Outer Boundary Condition                      C
C    ====              ========================                      C
C                                                                    C
C      1       No condition to be applied                            C
C      2       Function must vanish                                  C
C      3       First derivative must vanish                          C
C      4       Function AND first derivative must vanish             C
C      5       Function AND second derivative must vanish            C
C      6       rdf/dr - f(r) = 0                                     C
C      7       r df/dr + (l+1) f(r) = 0                              C
C                                                                    C
C____________________________________________________________________C
C                                                                    C
C Input variables :-                                                 C
C ===============                                                    C
C                                                                    C
C  Integer                                                           C
C  -------                                                           C
C                                                                    C
C     NR        : Total number of radial grid nodes.                 C
C     NNDS      : Number of nodes needed to take derivative.         C
C     NALF      : Number of nodes to left which are a linear         C
C                 combination of the other nodes.                    C
C     NARF      : Number of nodes to right which are a linear        C
C                 combination of the other nodes.                    C
C     L         : Spherical harmonic degree, l.                      C
C     IIBC      : Inner boundary flag - see above.                   C
C     IOBC      : Outer boundary flag - see above.                   C
C     NCFM      : Leading dimension of array DMAT etc.               C
C     IPCM      : Dimension ( NCFM ). Working array.                 C
C                                                                    C
C  Double Precision                                                  C
C  ----------------                                                  C
C                                                                    C
C     XARR      : Array of dimension ( NR ).                         C
C                  XARR( i ) contains the value of r at the          C
C                   i^{th} grid node.                                C
C                                                                    C
C     DMAT      : Dimension ( NCFM, NCFM ). See above.               C
C     WMAT      : Dimension ( NCFM, NCFM ). Working array.           C
C     WORK1     : Dimension ( NCFM ). Working array.                 C
C     WORK2     : Dimension ( NCFM ). Working array.                 C
C                                                                    C
C____________________________________________________________________C
C
C*********************************************************************
      SUBROUTINE LDGNMF( NR, NNDS, NALF, NARF, L, IIBC, IOBC, NCFM,
     1                   XARR, DMAT, WMAT, WORK1, WORK2, IPCM )
      IMPLICIT NONE
C____________________________________________________________________C
C Variable declarations - Parameters ................................C
      INTEGER NR, NNDS, NALF, NARF, L, IIBC, IOBC, NCFM,
     1        IPCM( NCFM )
      DOUBLE PRECISION XARR( NR ), DMAT( NCFM, NCFM),
     1                 WMAT( NCFM, NCFM), WORK1( NCFM ),
     2                 WORK2( NCFM )
C____________________________________________________________________C
C Variable declarations - Working variables .........................C
C
      INTEGER I, IFN, ILN, NDS2, ID1, IHND
      DOUBLE PRECISION ZERO, TOL, X0, FAC, QUOT
      PARAMETER ( ZERO = 0.0d0, TOL = 1.0d-9 )
C ifn is first linearly independent node
C iln is last linearly independent node
C ihnd is the highest number derivative which
C will be needed to be calculated during this
C routine (by GFDCFD)
C____________________________________________________________________C
C START OF PROGRAM **************************************************C
C____________________________________________________________________C
C
C First check on values of NALF and NARF
C
      I = NALF*NARF
      IF ( I.NE.0 .OR. (NALF.LT.0) .OR. (NALF.GT.2) .OR.
     1                 (NARF.LT.0) .OR. (NARF.GT.2)       ) THEN
        PRINT *,' Subroutine LDGNMF.'
        PRINT *,' NALF = ', NALF
        PRINT *,' NARF = ', NARF
        PRINT *,' Program aborted.'
        STOP
      ENDIF
C     .
C     . Ok - the number of affected nodes is valid
C     . so let's proceed. First, zero the matrix DMAT
C     .
      I = 0
      CALL MATOP( DMAT, ZERO, NCFM, NCFM, I )
C     .
C     . Now add the diagonal elements for the
C     . linearly independent variables
C     .
      IF ( NALF.EQ.0 ) IFN = 1
      IF ( NALF.EQ.1 ) IFN = 2
      IF ( NALF.EQ.2 ) IFN = 3
C     .
      IF ( NARF.EQ.0 ) ILN = NNDS
      IF ( NARF.EQ.1 ) ILN = NNDS - 1
      IF ( NARF.EQ.2 ) ILN = NNDS - 2
C     .
      DO I = IFN, ILN
        DMAT( I, I ) = 1.0d0
      ENDDO
C     .
C     . We can now return if the identity matrix is required
C     .
      IF ( NALF.EQ.0 .AND. NARF.EQ.0 ) RETURN
C     .
C     .
C     . OK - now decide whether we are doing inner or
C     . outer boundary
C     .
      IF ( NARF.EQ.0 ) THEN
C       .
C       . We are considering the inner boundary
C       . First, return if IIBC = 1 or 2
C       .
        IF ( IIBC.EQ.1 .OR. IIBC.EQ.2 ) RETURN
C       .
C       . So we need to calculate deriv.s at
C       . the inner boundary, XARR( 1 ).
C       . NDS2 is the number of nodes we can use at
C       . the inner boundary.
C       .
        IF ( IIBC.EQ.5 ) THEN
          IHND = 2
        ELSE
          IHND = 1
        ENDIF
C       .
        NDS2 = 1 + NNDS - NALF
C       .
        IF ( (NDS2-1).LT.IHND ) THEN
           PRINT *,' Subroutine LDGNMF '
           PRINT *,' NDS2                = ', NDS2
           PRINT *,' Required derivative = ', IHND
           PRINT *,' Program aborted.'
           STOP
        ENDIF
C       .
        DO I = 1, NDS2
          WORK1( I ) = XARR( I )
        ENDDO
        X0 = XARR( 1 )
        CALL GFDCFD ( X0, WORK1, NDS2, WMAT, NCFM, IPCM, WORK2 )
C       .
C       . WMAT( m + 1, i ) now contains the coefficient
C       . with which you multiply f_i to get the m^{th}
C       . derivative of f at x_0.
C       .
C       . Consider the case IIBC = 3
C       . We require the first derivative = 0
C       .
        IF ( IIBC.EQ.3 ) THEN
C          .
C          . Just quickly check that NALF = 1
C          . It should not be anything else at this stage.
C          .
           IF ( NALF.NE.1 ) THEN
             PRINT *,' Subroutine LDGNMF '
             PRINT *,' IIBC = ',IIBC,' and NALF = ',NALF
             PRINT *,' Program aborted.'
             STOP
           ENDIF
C          .
C          . Check for imminent division by zero
C          .
           IF ( ABS( WMAT( 2, 1 ) ).LT.TOL ) THEN
             PRINT *,' Subroutine LDGNMF '
             PRINT *,' WMAT( 2, 1 ) = ',WMAT( 2, 1 )
             PRINT *,' Program aborted.'
             STOP
           ENDIF
C          .
           DO I = 2, NDS2
             DMAT( 1, I ) = (-1.0d0)*WMAT( 2, I )/WMAT( 2, 1 )
           ENDDO
           RETURN
C          .
        ENDIF
C       .
C       . Consider the case IIBC = 4(5)
C       . We need both the function and the first
C       . (second) derivative to vanish. We do not need to
C       . to change anything to make the function
C       . vanish so just need do the first (second) deriv. cond.
C       .
        IF ( IIBC.EQ.4 .OR. IIBC.EQ.5 ) THEN
C          .
C          . Check for imminent division by zero
C          .
           IF ( IIBC.EQ.4 ) ID1 = 2
           IF ( IIBC.EQ.5 ) ID1 = 3
C          .
           IF ( ABS( WMAT( ID1, 2 ) ).LT.TOL ) THEN
             PRINT *,' Subroutine LDGNMF '
             PRINT *,' WMAT(',ID1,', 2 ) = ',WMAT( ID1, 2 )
             PRINT *,' Program aborted.'
             STOP
           ENDIF
C          .
           DO I = 3, NDS2
             DMAT( NALF, NALF - 2 + I ) =
     1                    (-1.0d0)*WMAT( ID1, I )/WMAT( ID1, 2 )
           ENDDO
           RETURN
C          .
        ENDIF
C       .
C       . We may treat the cases IIBC.EQ.6 and IIBC.EQ.7
C       . together as they are both conditions of the form
C       . rdf/dr + FAC f(r) = 0
C       .
        IF ( IIBC.EQ.6 .OR. IIBC.EQ.7 ) THEN
C          .
C          . First we make an early escape if r_{inner bnd} = 0
C          . this is equivlent to setting the function to zero
C          . which is what the current status of the matrix is
C          .
           IF ( ABS( X0 ).LT.TOL ) RETURN
C          .
C          . OK - so it's non-trivial!
C          . First, let's allocate the correct value of FAC
C          .
           IF ( IIBC.EQ.6 ) FAC = -1.0d0
           IF ( IIBC.EQ.7 ) FAC = -1.0d0*DBLE( L )
C          .
C          . Just quickly check that NALF = 1
C          . It should not be anything else at this stage.
C          .
           IF ( NALF.NE.1 ) THEN
             PRINT *,' Subroutine LDGNMF '
             PRINT *,' IIBC = ',IIBC,' and NALF = ',NALF
             PRINT *,' Program aborted.'
             STOP
           ENDIF
C          .
C          . Check for imminent division by zero
C          . We can safely divide by X0 now since the
C          . X0 = 0.0 case has been discounted above.
C          .
           QUOT = FAC/X0 + WMAT( 2, 1 )
           IF ( ABS( QUOT ).LT.TOL ) THEN
             PRINT *,' Subroutine LDGNMF '
             PRINT *,' QUOT = ', QUOT
             PRINT *,' Program aborted.'
             STOP
           ENDIF
C          .
           DO I = 2, NDS2
             DMAT( 1, I ) = (-1.0d0)*WMAT( 2, I )/QUOT
           ENDDO
           RETURN
C          .
        ENDIF
C       .
      ELSE
C       .
C       . We are considering the outer boundary
C       . First, return if IOBC = 1 or 2
C       .
        IF ( IOBC.EQ.1 .OR. IOBC.EQ.2 ) RETURN
C       .
C       . So we need to calculate deriv.s at
C       . the outer boundary, XARR( NR ).
C       . NDS2 is the number of nodes we can use at
C       . the outer boundary.
C       .
        IF ( IOBC.EQ.5 ) THEN
          IHND = 2
        ELSE
          IHND = 1
        ENDIF
C       .
        NDS2 = 1 + NNDS - NARF
C       .
        IF ( (NDS2-1).LT.IHND ) THEN
           PRINT *,' Subroutine LDGNMF '
           PRINT *,' NDS2                = ', NDS2
           PRINT *,' Required derivative = ', IHND
           PRINT *,' Program aborted.'
           STOP
        ENDIF
C       .
        DO I = 1, NDS2
          WORK1( I ) = XARR( NR - NDS2 + I )
        ENDDO
        X0 = XARR( NR )
        CALL GFDCFD ( X0, WORK1, NDS2, WMAT, NCFM, IPCM, WORK2 )
C       .
C       . WMAT( m + 1, i ) now contains the coefficient
C       . with which you multiply f_{nr-nds2+i} to get the m^{th}
C       . derivative of f at x_0.
C       .
C       . Consider the case IOBC = 3
C       . We require the first derivative = 0
C       .
        IF ( IOBC.EQ.3 ) THEN
C          .
C          . Just quickly check that NARF = 1
C          . It should not be anything else at this stage.
C          .
           IF ( NARF.NE.1 ) THEN
             PRINT *,' Subroutine LDGNMF '
             PRINT *,' IOBC = ',IOBC,' and NARF = ',NARF
             PRINT *,' Program aborted.'
             STOP
           ENDIF
C          .
C          . Check for imminent division by zero
C          .
           IF ( ABS( WMAT( 2, NDS2 ) ).LT.TOL ) THEN
             PRINT *,' Subroutine LDGNMF '
             PRINT *,' WMAT( 2,',NDS2,') = ',WMAT( 2, NDS2 )
             PRINT *,' Program aborted.'
             STOP
           ENDIF
C          .
C          .
           DO I = 1, NDS2-1
             DMAT( NNDS, NNDS - NDS2 + I ) =
     1             (-1.0d0)*WMAT( 2, I )/WMAT( 2, NDS2 )
           ENDDO
           RETURN
C          .
        ENDIF
C       .
C       . Consider the case IOBC = 4(5)
C       . We need both the function and the first
C       . (second) derivative to vanish. We do not need to
C       . to change anything to make the function
C       . vanish so just need do the first (second) deriv. cond.
C       .
        IF ( IOBC.EQ.4 .OR. IOBC.EQ.5 ) THEN
C          .
C          . Check for imminent division by zero
C          .
           IF ( IOBC.EQ.4 ) ID1 = 2
           IF ( IOBC.EQ.5 ) ID1 = 3
C          .
           IF ( ABS( WMAT( 2, NDS2-1 ) ).LT.TOL ) THEN
             PRINT *,' Subroutine LDGNMF '
             PRINT *,' WMAT(',ID1,',',NDS2-1,') = ',WMAT( ID1,NDS2-1)
             PRINT *,' Program aborted.'
             STOP
           ENDIF
C          .
           DO I = 1, NDS2-2
             DMAT( NNDS + 1 - NARF, NNDS - NDS2 - NARF + 2 + I ) =
     1                    (-1.0d0)*WMAT( ID1, I )/WMAT( ID1, NDS2-1)
           ENDDO
           RETURN
C          .
        ENDIF
C       .
C       . We may treat the cases IOBC.EQ.6 and IOBC.EQ.7
C       . together as they are both conditions of the form
C       . rdf/dr + FAC f(r) = 0
C       .
        IF ( IOBC.EQ.6 .OR. IOBC.EQ.7 ) THEN
C          .
C          . Just quickly check that NARF = 1
C          . It should not be anything else at this stage.
C          .
           IF ( NARF.NE.1 ) THEN
             PRINT *,' Subroutine LDGNMF '
             PRINT *,' IOBC = ',IOBC,' and NARF = ',NARF
             PRINT *,' Program aborted.'
             STOP
           ENDIF
C          .
C          . First, let's allocate the correct value of FAC
C          .
           IF ( IOBC.EQ.6 ) FAC = -1.0d0
           IF ( IOBC.EQ.7 ) FAC = DBLE( L + 1 )
C          .
C          . Check for imminent division by zero
C          . (Once again, we safely divide by X0)
C          .
           QUOT = FAC/X0 + WMAT( 2, NDS2 )
           IF ( ABS( QUOT ).LT.TOL ) THEN
             PRINT *,' Subroutine LDGNMF '
             PRINT *,' QUOT = ', QUOT
             PRINT *,' Program aborted.'
             STOP
           ENDIF
C          .
C          .
           DO I = 1, NDS2-1
             DMAT( NNDS, NNDS - NDS2 + I ) =
     1             (-1.0d0)*WMAT( 2, I )/QUOT
           ENDDO
           RETURN
C          .
        ENDIF
C       .
      ENDIF
C     .
      PRINT *,' Subroutine LDGNMF '
      PRINT *,' Problem with inputs.'
      PRINT *,' IIBC = ', IIBC
      PRINT *,' IOBC = ', IOBC
      PRINT *,' Program aborted.'
      STOP
      END
C*********************************************************************
C*********************************************************************
C subroutine General Finite Difference Coefficient Find **************
C            -       -      -          -           -    **************
C Steve Gibbons Mon Sep 20 16:57:54 BST 1999                         C
C____________________________________________________________________C
C                                                                    C
C Given a value of X and the values of x_i at NNDS distinct points,  C
C (X need not necessarily be one of the x_i) then the array COEFM    C
C is returned with the finite difference coefficients such that      C
C the ND^{th} derivative of a function f, evaluated at x = X,        C
C is given by                                                        C
C                                                                    C
C f^{ ND }( X ) = \sum_{j = 1, NNDS} COEFM( ND + 1, j )*f( x_j )     C
C                                                                    C
C This is a general version of the routine FDCINV which is valid     C
C only for equally spaced grid nodes.                                C
C                                                                    C
C Coefficients for up to the (NNDS-1)^{th} derivative are given      C
C although care must be taken to ensure the highest derivatives      C
C are sufficiently accurate.                                         C
C                                                                    C
C____________________________________________________________________C
C                                                                    C
C Input variables :-                                                 C
C ===============                                                    C
C  Double Precision                                                  C
C  ----------------                                                  C
C     X         : Abscissa at which derivatives are to be            C
C                  evaluated.                                        C
C     XARR      : Array of dimension ( NNDS ).                       C
C                  XARR( i ) contains the value of x/r at the        C
C                   i^{th} grid node.                                C
C                                                                    C
C     COEFM     : Dimension ( NCFM, NCFM).                           C
C     WORK      : Workspace array for LAPACK inversion routine.      C
C                 Dimension ( NCFM )                                 C
C                                                                    C
C  Integer                                                           C
C  -------                                                           C
C                                                                    C
C     NNDS      : Number of grid nodes.                              C
C     NCFM      : Leading order of coefficient matrix.               C
C                                                                    C
C     IPCM      : Work array for LAPACK routines to perform          C
C                 pivotting in the matrix inversion.                 C
C                 Dimension ( NCFM )                                 C
C                                                                    C
C____________________________________________________________________C
C
C*********************************************************************
      SUBROUTINE GFDCFD ( X, XARR, NNDS, COEFM, NCFM, IPCM, WORK)
      IMPLICIT NONE
C____________________________________________________________________C
C Variable declarations - Parameters ................................C
      INTEGER NNDS, NCFM, IPCM( NCFM )
      DOUBLE PRECISION X, COEFM( NCFM, NCFM ), WORK( NCFM ),
     1                 XARR( NNDS )
C____________________________________________________________________C
C Variable declarations - Working variables .........................C
C
      INTEGER I, J, NDER, INODE, INFO, ICOL, IROW
      DOUBLE PRECISION DZERO, LOW, FAC
      PARAMETER ( DZERO = 0.0d0, LOW = 1.0d-8 )
C____________________________________________________________________C
C START OF PROGRAM **************************************************C
C____________________________________________________________________C
C
C Check bounds of integer parameters
C
      IF ( NNDS.GT.NCFM ) THEN
         PRINT *,' Subroutine GFDCFD: '
         PRINT *,' NNDS = ', NNDS,'. NCFM = ', NCFM
         PRINT *,' Program aborted.'
         STOP
      ENDIF
C
C Calculate the distances h_i = ( x_i - X )
C These h_i must be distinct and this will be
C checked for - otherwise matrix is singular.
C We can store the h_i in WORK as this will not
C be needed until the inversion, by which time
C it will not be needed by us.
C
      DO I = 1, NNDS
        WORK( I ) = XARR( I ) - X
      ENDDO
C
C Now check for the uniqueness of the points ...
C
      DO I = 1, NNDS - 1
        DO J = I + 1, NNDS
          IF ( ABS( WORK( I ) - WORK( J ) ).LT.LOW ) THEN
            PRINT *,' Subroutine GFDCFD.'
            PRINT *,' X values ',I,' and ',J,' are'
            PRINT *,' identical.'
            PRINT *,' Program aborted.'
            STOP
          ENDIF
        ENDDO
      ENDDO
C
C____________________________________________________________________C
C Parameters are ok so let's zero COEFM
C
      I = 0
      CALL MATOP( COEFM, DZERO, NCFM, NCFM, I )
C
C (nder+1) is the number of the matrix column being filled in.
C inode is the number of the matrix row being filled in.
C
      DO NDER = 0, NNDS - 1
        ICOL = NDER + 1
        DO INODE = 1, NNDS
          IROW = INODE
          IF ( NDER.EQ.0 ) THEN
            COEFM( IROW, ICOL )  = 1.0d0
          ELSE
            FAC = WORK( INODE )/DBLE( NDER )
            COEFM( IROW, ICOL ) = COEFM( IROW, ICOL-1 )*FAC
          ENDIF
        ENDDO
      ENDDO
C
C Ok - this matrix is now ready for inversion -
C For this we use the LAPACK routines DGETRF and DGETRI
C First perform LU decomposition
C
      CALL DGETRF( NNDS, NNDS, COEFM, NCFM, IPCM, INFO )
C
C     . Check that LU decomposition has gone without
C     . problem.
C     .
C
      IF ( INFO.NE.0 ) THEN
         PRINT *,' Subroutine GFDCFD.'
         PRINT *,' The LAPACK subroutine DGETRF has'
         PRINT *,' returned ',INFO,' as a value of '
         PRINT *,' INFO in LU decomposition of COEFM matrix.'
         PRINT *,' Program aborted.'
         STOP
      ENDIF
C     .
C     . Now compute the inverse with the LAPACK routine
C     . DGETRI.
C     .
      CALL DGETRI( NNDS, COEFM, NCFM, IPCM, WORK, NCFM, INFO )
C     .
C     . Check that inversion has gone without problem.
C     .
      IF ( INFO.NE.0 ) THEN
         PRINT *,' Subroutine GFDCFD.'
         PRINT *,' The LAPACK subroutine DGETRI has'
         PRINT *,' returned ',INFO,' as a value of '
         PRINT *,' INFO in inversion of COEFM matrix.'
         PRINT *,' Program aborted.'
         STOP
      ENDIF
C     .
      RETURN
      END
C*********************************************************************
C*********************************************************************
C subroutine Adapted Solution Vector DeRivative **********************
C            -       -        -      - -        **********************
C Steve Gibbons Mon Oct 25 15:17:59 BST 1999                         C
C____________________________________________________________________C
C                                                                    C
C V is a solution vector with NH (= inarr(3) ) harmonic radial       C
C functions and NR (= inarr(2) ) grid nodes in each function.        C
C The position of the j^{th} node of the i^{th} harmonic is given by C
C INDFUN( j, i, INARR) and the radial value is given by              C
C  XARR( j ) - an array which is passed to the routine FDCMBD        C
C in order to calculate the finite difference coefficients, FDCM.    C
C XARR itself is not referenced by ASVDR.                            C
C ASVDR returns the radial derivatives 0, 1, ..., IHD of radial      C
C function IH evaluated at node IR.                                  C
C                                                                    C
C IFORMF = INARR(1) should be either 3 or 4 since this is the        C
C arbitrarily spaced mesh version of the code.                       C
C                                                                    C
C NBN is the maximum number of nodes on                              C
C either side which maybe used in central differences.               C
C For instance, if to calculate the derivative of                    C
C f at r_j, you may use the values of f at r = r_{j-2}, r_{j-1},     C
C r_j, r_{j+1} and r_{j+2} then NBN = 2. The value of IHD is checked C
C only for being positive and no greater than NDRVS (the             C
C number of the highest derivative for which coefficients            C
C are stored by the array SVFDC), as SVFDC must be calculated in     C
C advance by a call to SVFDCF which checks NDRVS against             C
C the physical restrictions imposed by the value of NBN.             C
C                                                                    C
C NDRVM restricts the size of NDRVS and is a defining parameter      C
C of the array SVFDCF.                                               C
C                                                                    C
C NBN must be as supplied to SVFDCF.                                 C
C                                                                    C
C IS must be supplied to indicate the finite difference scheme       C
C being employed.                                                    C
C                                                                    C
C ASVDR will use whichever coefficients SVFDC( k, j, nd + 1, K ),    C
C that SVFDCF formed with K = IS.                                    C
C____________________________________________________________________C
C                                                                    C
C Input variables :-                                                 C
C ===============                                                    C
C  Integer                                                           C
C  -------                                                           C
C     IR        : Number of radial grid node.                        C
C     IS        : Number of finite difference scheme used.           C
C     IH        : Number of radial function (harmonic).              C
C     NBN       : Number of bounding nodes. See above.               C
C     IHD       : Highest derivative requested.                      C
C     NFDCM     : Leading dimension of SVFDC. At least (2*NBN+1)     C
C     NR        : Number of radial grid nodes in each function.      C
C     NDRVS     : Number of highest derivative for which             C
C                  coefficients are stored by the array SVFDC.       C
C     NDRVM     : Limit on NDRVS. Array bound for SVFDC.             C
C     INARR     : Integer array dimension ( * ).                     C
C                  Elements may be arbitrary except for              C
C                  INARR( 1 ) = IFORMF - flag for vector format.     C
C                                        See INDFUN                  C
C                  INARR( 2 ) = NRR. Must be consistent with NR.     C
C                  INARR( 3 ) = NH = total number of radial func.s   C
C     NDCS      : Maximum distinct finite difference schemes         C
C                  stored by the array SVFDC.                        C
C                                                                    C
C  Double Precision                                                  C
C  ----------------                                                  C
C                                                                    C
C     V         : Solution vector. Dim ( * ) but length atleast      C
C                  NR*NH.                                            C
C     DERV      : Derivatives. Dim ( * ) but length atleast IHD      C
C                  DERV( i ) is returned containing the value of     C
C                  the i^[th} derivative.                            C
C                                                                    C
C     SVFDC     : Finite difference coefficient matrix.              C
C                  Dimension ( NFDCM, NR, NDRVM+1, NDCS ).           C
C                   Array is generated by the routine svfdcf         C
C                 See documentation for SVFDCF for details.          C
C                                                                    C
C____________________________________________________________________C
C
C*********************************************************************
      SUBROUTINE ASVDR ( V, IR, IS, IH, NBN, IHD, NFDCM, NR, NDRVS,
     1                   NDRVM, DERV, INARR, SVFDC, NDCS )
      IMPLICIT NONE
C____________________________________________________________________C
C Variable declarations - Parameters ................................C
      INTEGER IR, IS, IH, NBN, IHD, NFDCM, NR, NDRVS, NDRVM,
     1        INARR( * ), NDCS
      DOUBLE PRECISION V( * ), DERV( * ), 
     1                 SVFDC( NFDCM, NR, NDRVM+1, NDCS )
C____________________________________________________________________C
C Variable declarations - Working variables .........................C
C
      DOUBLE PRECISION COEF
      INTEGER INODE, ILN, IRN, INDFUN, ID, IND, NRR, IFORMF, IK, ID1
C____________________________________________________________________C
C START OF PROGRAM **************************************************C
C____________________________________________________________________C
C
C Check input parameters ...
C
      IFORMF = INARR( 1 )
      NRR    = INARR( 2 )
C
      IF ( NRR.NE.NR ) THEN
         PRINT *,' Subroutine ASVDR.'
         PRINT *,' INARR( 2 ) = ', NRR
         PRINT *,' NR = ', NR
         PRINT *,' Program aborted.'
         STOP
      ENDIF
C
      IF ( IFORMF.NE.3 .AND. IFORMF.NE.4 ) THEN
         PRINT *,' Subroutine ASVDR.'
         PRINT *,' INARR( 1 ) = ', IFORMF
         PRINT *,' This is an irregular grid routine.'
         PRINT *,' Program aborted.'
         STOP
      ENDIF
C
      IF ( IR.LT.1 .OR. IR.GT.NR ) THEN
        PRINT *,' Subroutine ASVDR.'
        PRINT *,' IR   = ', IR
        PRINT *,' Program aborted.'
        STOP
      ENDIF
C     .
      IF ( IHD.LT.0 .OR. IHD.GT.NDRVS .OR. NDRVS.GT.NDRVM ) THEN
        PRINT *,' Subroutine ASVDR.'
        PRINT *,' IHD   = ', IHD
        PRINT *,' NDRVS = ', NDRVS
        PRINT *,' NDRVM = ', NDRVM
        PRINT *,' Program aborted.'
        STOP
      ENDIF
C     .
      IF ( IS.LT.1 .OR. IS.GT.NDCS ) THEN
        PRINT *,' Subroutine ASVDR.'
        PRINT *,' IS    = ', IS
        PRINT *,' NDCS  = ', NDCS
        PRINT *,' Program aborted.'
        STOP
      ENDIF
C     .
C     . Calculate furthest left and furthest right mode
C     . to be used to form derivative.
C     .
      ILN = MAX(  1, IR - NBN )
      IRN = MIN( NR, IR + NBN )
C
      DO ID = 0, IHD
        ID1 = ID + 1
        DERV( ID1 ) = 0.0d0
        DO INODE = ILN, IRN
          IK = INODE - IR + NBN + 1
          COEF = SVFDC( IK, IR, ID1, IS )
          IND = INDFUN( INODE, IH, INARR )
          DERV( ID1 ) = DERV( ID1 ) + COEF*V( IND )
        ENDDO
      ENDDO
C     .
      RETURN
      END
C*********************************************************************

C*********************************************************************
C subroutine VECtor OPeration ****************************************
C Steve Gibbons 22.4.97 Fills vector with a constant, multiplies a   C
C                       vector by a constant or adds a constant.     C
C____________________________________________________________________C
C Input variables :-                                                 C
C ===============                                                    C
C  Integer                                                           C
C  -------                                                           C
C     IOP	: Type of operation required.                        C
C                  IOP=0  -->  Each element of the vector = CONST    C
C                  IOP=1  -->  Each el. is multiplied by CONST       C
C                  IOP=2  -->  Each el. is added to CONST            C
C     N		: Length of the vector.                              C
C  Double Precision                                                  C
C  ----------------                                                  C
C     VEC	: Vector - dimension ( N )                           C
C     CONST     : Double precision constant.                         C
C____________________________________________________________________C
C
C*********************************************************************
      SUBROUTINE VECOP ( VEC, CONST, N, IOP )
      IMPLICIT NONE
C____________________________________________________________________C
C Variable declarations - Parameters ................................C
      INTEGER N, IOP
      DOUBLE PRECISION VEC( N ), CONST
C____________________________________________________________________C
C Variable declarations - Working variables .........................C
      INTEGER I
C____________________________________________________________________C
C START OF PROGRAM **************************************************C
C____________________________________________________________________C
C
C First do the case of making an element constant.
      IF ( IOP.EQ.0 ) THEN
         DO I = 1, N
            VEC ( I ) = CONST
         ENDDO
         RETURN
      ENDIF
C Now do multiplying a vector
      IF ( IOP.EQ.1 ) THEN
         DO I = 1, N
            VEC ( I ) = VEC( I )*CONST
         ENDDO
         RETURN
      ENDIF
C Now do adding a vector
      IF ( IOP.EQ.2 ) THEN
         DO I = 1, N
            VEC ( I ) = VEC( I ) + CONST
         ENDDO
         RETURN
      ENDIF
C____________________________________________________________________C

      PRINT *,' Subroutine VECOP. IOP must be 0, 1 or 2.'
      PRINT *,'Program aborted.'
      STOP
      END
C*********************************************************************


C*********************************************************************
C subroutine Fast Fourier Transform ReaL Version *********************
C            -    -       -         -  - -       *********************
C Steve Gibbons 22.4.97 (Based on FOUR1 from Numerical Recipes       C
C      I've made slight alterations for double precision, error      C
C      checking at the start and rescaling the coefficients after    C
C      the forward transform.)                                       C
C____________________________________________________________________C
C Input variables :-                                                 C
C ===============                                                    C
C  Integer                                                           C
C  -------                                                           C
C     NN	: Number of data to be transformed.		     C
C                  NN MUST be a power of two. This is checked for by C
C                  the subroutine POWTWO.                            C
C                                                                    C
C  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -       C
C                                                                    C
C    The real function f( PHI ) must be periodic, satisfying         C
C    f( PHI ) = f( PHI + 2*PI )                                      C
C                                                                    C
C    f( PHI ) has the formulation                                    C
C                                                                    C
C    f( PHI ) = sum_{m=0}^M [ c_m cos(m PHI) + s_m sin( m PHI) ]     C
C                                                                    C
C    FFTRLV will transform between the coefficients {c_m, s_m}       C
C    and discretised values of f at NN equally spaced points PHI_i   C
C    with                                                            C
C                                                                    C
C    PHI_i = (i-1)*DELTA_phi          and                            C
C                                                                    C
C    DELTA_phi = 2*PI/dble( NN )                                     C
C                                                                    C
C  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -       C
C                                                                    C
C     ISIGN	:                                                    C
C                                                                    C
C    If ISIGN = 1, on entry DATA( 2*I - 1 ) must contain the value   C
C    of f at the point PHI = PHI_i as defined above.                 C
C                                                                    C
C    On exit, DATA( 2*m + 1 ) will contain the coeff. c_m and        C
C             DATA( 2*m + 2 ) will contain the coeff. s_m and        C
C                                                                    C
C    for m = 0, M.                                                   C
C      For accuracy, M **MUST** be less than NN/2                    C
C    (This is the Nyquist criterion - see Numerical Recipes ).       C
C                                                                    C
C    If ISIGN = -1, DATA must be zero except for the values          C
C    c_m contained in DATA( 2*m + 1 )  and                           C
C    s_m contained in DATA( 2*m + 2 ) on input.                      C
C                                                                    C
C    On output, f( PHI_i ) will be contained in DATA( 2*I - 1 ).     C
C                                                                    C
C    Note that this routine SCALES the coefficients after            C
C    transforming, which the Numerical Recipes version doesn't.      C
C                                                                    C
C  Double Precision                                                  C
C  ----------------                                                  C
C     DATA	: Array of dimension (2*NN). See above ...           C
C____________________________________________________________________C
C
C*********************************************************************
      SUBROUTINE FFTRLV ( DATA, NN, ISIGN )
      IMPLICIT NONE
C____________________________________________________________________C
C Variable declarations - Parameters ................................C
      INTEGER NN, ISIGN
      DOUBLE PRECISION DATA( 2 * NN )
C____________________________________________________________________C
C Variable declarations - Working variables .........................C
      DOUBLE PRECISION TEMPR,TEMPI,THETA,WPR,WPI,WR,WI,WTEMP,FAC
      INTEGER I,J,N,M,MMAX,ISTEP
c     LOGICAL POWT
      DOUBLE PRECISION PI, HTHETA
      PARAMETER (PI=3.14159265358979312D0)
C____________________________________________________________________C
C START OF PROGRAM **************************************************C
C____________________________________________________________________C
C 
C  SJG Mon Feb  5 18:18:33 WET 2001
C  I have commented out the power of two checking
C  for speed: reinstate if necessary:
C 
C First check that the number NN supplied is infact a power of 2.
c     CALL POWTWO ( NN, POWT)
c     IF ( .NOT. POWT) THEN
c        PRINT *,' Subroutine FFTRLV. NN is not a power of 2.'
c        PRINT *,'Program aborted.'
c        STOP
c     ENDIF
C Now check that ISIGN = 1 or -1.
      IF ( ISIGN*ISIGN .NE. 1) THEN
         PRINT *,' Subroutine FFTRLV. ISIGN must be 1 or -1.'
         PRINT *,'Program aborted.'
         STOP
      ENDIF
C____________________________________________________________________C
C Here we will zero the imaginary part of the function
C just incase it contains non-zero data which may harm the
C outcome.
C
      IF ( ISIGN.EQ.1 ) THEN
        DO I = 1, NN
          DATA( 2*I ) = 0.0d0
        ENDDO
      ENDIF
C
C____________________________________________________________________C
C
      N = 2 * NN
      J = 1
      DO I = 1, N, 2
         IF ( J.GT.I ) THEN
            TEMPR = DATA( J )
            TEMPI = DATA( J+1 )
            DATA( J ) = DATA( I )
            DATA( J+1 ) = DATA( I+1 )
            DATA( I ) = TEMPR
            DATA( I+1 ) = TEMPI
         ENDIF
C ... now calculate inverse bit map for next value of I.
         M = N/2
 500     IF ( M.GE.2 .AND. J.GT.M ) THEN
            J = J-M
            M = M/2
            GOTO 500
         ENDIF
         J = J + M
      ENDDO
C ............. now do the real part of transform.
      MMAX = 2
 502  IF ( N.GT.MMAX ) THEN
         ISTEP = 2*MMAX
         HTHETA = PI/DBLE(ISIGN*MMAX)
c        THETA = 2.0d0*PI/DBLE(ISIGN*MMAX)
C set THETA temporarily
         THETA = DSIN(HTHETA)
         WPR = -2.0d0*THETA*THETA
c        WPR = -2.0d0*DSIN(HTHETA)**2
         THETA = 2.0d0*HTHETA
         WPI = DSIN(THETA)
         WR = 1.0d0
         WI = 0.0d0
         DO M = 1, MMAX, 2
            DO I = M, N, ISTEP
               J = I + MMAX
               TEMPR = WR*DATA( J ) - WI*DATA( J + 1 )
               TEMPI = WR*DATA( J + 1 ) + WI*DATA( J )
               DATA( J ) = DATA( I ) - TEMPR
               DATA( J + 1 ) = DATA( I + 1 ) - TEMPI
               DATA( I ) = DATA( I ) + TEMPR
               DATA( I + 1 ) = DATA( I + 1 ) + TEMPI
            ENDDO
            WTEMP = WR
            WR = WR*WPR - WI*WPI + WR
            WI = WI*WPR + WTEMP*WPI + WI
         ENDDO
         MMAX = ISTEP
         GOTO 502
      ENDIF
C
      IF ( ISIGN.EQ.1 ) THEN
        FAC = 2.0d0/DBLE( NN )
        DO I = 1, NN
          DATA( I ) = FAC*DATA( I )
        ENDDO
        DO I = NN+1, 2*NN
          DATA( I ) = 0.0d0
        ENDDO
      ENDIF
C
      RETURN
      END
C*********************************************************************

C*********************************************************************
C subroutine MATrix OPeration ****************************************
C Steve Gibbons 23.4.97 Does operation on a two-dimensional array.   C
C                                                Can set equal to a  C
C                       constant; multiply by a constant or have a   C
C                       constant added to it.                        C
C____________________________________________________________________C
C Input variables :-                                                 C
C ===============                                                    C
C  Integer                                                           C
C  -------                                                           C
C     IOP	: Type of operation to be done.                      C
C                    WHOLE MATRIX OPERATIONS                         C
C                  IOP=0  -->  Each element of the matrix = CONST    C
C                  IOP=1  -->  Each el. is multiplied by CONST       C
C                  IOP=2  -->  Each el. is added to CONST            C
C     NDIM1     : First dimension of the matrix.	             C
C     NDIM2     : Second dimension of the matrix.	             C
C  Double Precision                                                  C
C  ----------------                                                  C
C     MAT	: Matrix with dimension ( NDIM1, NDIM2 )             C
C     CONST     : Double precision constant.                         C
C____________________________________________________________________C
C
C*********************************************************************
      SUBROUTINE MATOP ( MAT, CONST, NDIM1, NDIM2, IOP)
      IMPLICIT NONE
C____________________________________________________________________C
C Variable declarations - Parameters ................................C
      INTEGER NDIM1, NDIM2, IOP
      DOUBLE PRECISION MAT ( NDIM1, NDIM2 ), CONST
C____________________________________________________________________C
C Variable declarations - Working variables .........................C
      INTEGER I, J
C____________________________________________________________________C
C START OF PROGRAM **************************************************C
C____________________________________________________________________C
C
C First do case IOP=0 
      IF ( IOP.EQ.0 ) THEN
         DO J = 1, NDIM2
            DO I = 1, NDIM1
               MAT ( I, J) = CONST
            ENDDO
         ENDDO
         RETURN
      ENDIF
C Now do case IOP=1
      IF ( IOP.EQ.1 ) THEN
         DO J = 1, NDIM2
            DO I = 1, NDIM1
               MAT ( I, J) = MAT ( I, J)*CONST
            ENDDO
         ENDDO
         RETURN
      ENDIF
C Now do case IOP=2
      IF ( IOP.EQ.2 ) THEN
         DO J = 1, NDIM2
            DO I = 1, NDIM1
               MAT ( I, J) = MAT ( I, J) + CONST
            ENDDO
         ENDDO
         RETURN
      ENDIF
      PRINT *,' Subroutine MATOP. IOP must be 0,1 or 2.'
      STOP
      END
C*********************************************************************
C*********************************************************************
C subroutine POWer of TWO checker ************************************
C            ---      ---         ************************************
C Steve Gibbons 11.4.97        					     C
C____________________________________________________________________C
C Checks an integer INPUT and returns a .TRUE. logical variable POT  C
C if and only if INPUT is a power of 2.				     C
C____________________________________________________________________C
C Input variables :-                                                 C
C ===============                                                    C
C  Integer                                                           C
C  -------                                                           C
C     INPUT	: Integer to be tested				     C
C____________________________________________________________________C
C Output :-                                                          C
C ======                                                             C
C  Logical 							     C
C  ------- 							     C
C     POT	: Power Of Two ? 				     C
C____________________________________________________________________C
C*********************************************************************
      SUBROUTINE POWTWO(INPUT,POT)
      IMPLICIT NONE
C____________________________________________________________________C
C Variable declarations - Parameters ................................C
      INTEGER INPUT
      LOGICAL POT
C____________________________________________________________________C
C Variable declarations - Working variables .........................C
      INTEGER ITWO,N,IREM
      PARAMETER (ITWO=2)
C____________________________________________________________________C
C START OF PROGRAM **************************************************C
C____________________________________________________________________C
C First of all put N = INPUT so that INPUT won't be altered
      N=INPUT
C Then check that N is atleast equal to 2 - otherwise
C it is obviously not a power of 2.
      IF ( N.LT.2 ) THEN
         POT=.FALSE.
         RETURN
      ENDIF
C Begin our iteration ....
 500  CONTINUE
      IREM = MOD ( N , ITWO )
      IF ( IREM.NE.0 ) THEN
         POT=.FALSE.
         RETURN
      ENDIF
      IF ( N.EQ.2 ) THEN
         POT=.TRUE.
         RETURN
      ELSE
         N=N/2
         GOTO 500
      ENDIF
      END
C*********************************************************************
C*********************************************************************
C function Element of Matrix MULTiplication evaluate *****************
C          -          -      ----                    *****************
C Steve Gibbons Fri Oct 22 17:23:46 BST 1999                         C
C____________________________________________________________________C
C                                                                    C
C Let the matrices B and C be respectively 'm by k' and 'k by n'     C
C double precision matrices. Then A = B C is an 'm by n' matrix      C
C with a_{ij} = \sum_{l=1,k} b_{il} c_{lj}                           C
C                                                                    C
C EMMULT returns the value of A( I, J)                               C
C____________________________________________________________________C
C                                                                    C
C Input variables :-                                                 C
C ===============                                                    C
C                                                                    C
C  Integer                                                           C
C  -------                                                           C
C     I         : Row containing desired element of A.               C
C     J         : Column containing desired element of A.            C
C     LDB       : Leading dimension of matrix B.                     C
C     LDC       : Leading dimension of matrix C.                     C
C     K         : Number of columns of matrix B and                  C
C                 Number of rows of matrix C.                        C
C                                                                    C
C  Double Precision                                                  C
C  ----------------                                                  C
C                                                                    C
C     B         : First matrix. Dim ( LDB, * )                       C
C     C         : First matrix. Dim ( LDC, * )                       C
C                                                                    C
C____________________________________________________________________C
C
C*********************************************************************
      FUNCTION EMMULT( I, J, LDB, LDC, K, B, C )
      IMPLICIT NONE
C____________________________________________________________________C
C Variable declarations - Parameters ................................C
      INTEGER I, J, LDB, LDC, K
      DOUBLE PRECISION EMMULT, B( LDB, * ), C( LDB, * )
C____________________________________________________________________C
C Variable declarations - Working variables .........................C
      INTEGER L
C____________________________________________________________________C
C START OF PROGRAM **************************************************C
C____________________________________________________________________C
C
C First check the validity of the integer inputs.
C K must not exceed LDC
C I must not exceed LDB
C     .
      IF ( K.GT.LDC .OR. I.GT.LDB ) THEN
        PRINT *,' Function EMMULT.'
        PRINT *,' K    = ', K,' I    = ', I
        PRINT *,' LDB  = ', LDB,' LDC  = ', LDC
        PRINT *,' Program aborted.'
        STOP
      ENDIF
C     .
      EMMULT = 0.0d0
C     . 
      DO L = 1, K
        EMMULT = EMMULT + B( I, L )*C( L, J )
      ENDDO
C     .
      RETURN
      END
C*********************************************************************

C*********************************************************************
C function PMM *******************************************************
C Steve Gibbons 16.4.97                                              C
C____________________________________________________________________C
C Gives the Schmidt Normalised Legendre Function P_m^m ( X )         C
C from eqn. 175 in my notes ie                                       C
C                                                                    C
C                   ( 2m - 1)!! (1- XX)^(m/2) * SQRT (2.0d0 )        C
C   P_m^m( X ) =  ---------------------------------------------      C
C                        SQRT ( (2m)! )                              C
C                                                                    C
C       for m non-zero and                                           C
C                                                                    C
C   P_0^0( X ) = 1.0d0                                               C
C                                                                    C
C N.B. The double factorial sign means the product of all ODD        C
C integers between 1 and ( 2m - 1 ).                                 C
C Best to use the form in eq 184 to calculate PMM.                   C
C____________________________________________________________________C
C Input Variables :-                                                 C
C ===============                                                    C
C  Integer                                                           C
C  -------                                                           C
C     M         : Well it's M isn't it!                              C
C  Double Precision                                                  C
C  ----------------                                                  C
C     SINE      : Sin (theta)                                        C
C____________________________________________________________________C
C
C*********************************************************************
      FUNCTION PMM ( M, S )
      IMPLICIT NONE
      DOUBLE PRECISION PMM
C____________________________________________________________________C
C Variable declarations - Parameters ................................C
      INTEGER M
      DOUBLE PRECISION S
C____________________________________________________________________C
C Variable declarations - Working variables .........................C
      INTEGER I
      DOUBLE PRECISION TWOI
C____________________________________________________________________C
C START OF FUNCTION *************************************************C
C____________________________________________________________________C
C
      IF ( M.LT.0 ) THEN
         PRINT *,' Function PMM. M < 0 error. Stopped.'
         STOP
      ENDIF
      IF ( M.EQ.0 ) THEN
         PMM = 1.0d0
         RETURN
      ENDIF
C ................................. so M is greater than 0
      PMM = DSQRT ( 2.0d0 )
      DO I = 1, M
         TWOI = 2.0d0*DBLE(I)
         PMM = PMM * S * DSQRT ( (TWOI - 1.0d0 )/TWOI )
      ENDDO
      RETURN
      END
C*********************************************************************
C*********************************************************************
C function PMM1 ******************************************************
C Steve Gibbons 16.4.97                                              C
C____________________________________________________________________C
C Evaluates the Schmidt Normalised Legendre Function P_(m+1)^m (X)   C
C according to equation 179 in my notes ; i.e.                       C
C                                                                    C
C    P_(m+1)^m (X) = SQRT( 2m+1 ).X.P^m_m(X)                         C
C                                                                    C
C____________________________________________________________________C
C Input Variables :-                                                 C
C ===============                                                    C
C  Integer                                                           C
C  -------                                                           C
C     M         : Well it's M isn't it!                              C
C  Double Precision                                                  C
C  ----------------                                                  C
C     X         : Cos(theta)                                         C
C     PMM0      : P_m^m(X) as evaluated by Function PMM.             C
C____________________________________________________________________C
C
C*********************************************************************
      FUNCTION PMM1 ( M, X, PMM0 )
      IMPLICIT NONE
      DOUBLE PRECISION PMM1
C____________________________________________________________________C
C Variable declarations - Parameters ................................C
      INTEGER M
      DOUBLE PRECISION X, PMM0
C____________________________________________________________________C
C Variable declarations - Working variables .........................C
      DOUBLE PRECISION RM
C____________________________________________________________________C
C START OF FUNCTION *************************************************C
C____________________________________________________________________C
C
      IF ( M.LT.0 ) THEN
         PRINT *,' Function PMM1. M < 0 error. Stopped.'
         STOP
      ENDIF
      RM = DBLE( M )
      PMM1 = X*PMM0*DSQRT( 2.0d0*RM+1.0d0 )
      RETURN
      END
C*********************************************************************
C*********************************************************************
C function PLM *******************************************************
C Steve Gibbons 16.4.97                                              C
C____________________________________________________________________C
C Calculates the Schmidt Normalised Legendre Function P_l^m (x)      C
C given P_(l-1)^m and P_(l-2)^m according to equation 183 in my notesC
C i.e.                                                               C
C   P_l^m( X ) = { A * P_(l-1)^m - B * P_(l-2)^m }/C                 C
C                                                                    C
C where A = (2*l - 1)*X ,                                            C
C                                                                    C
C B = SQRT( (L+M-1)*(L-M-1) ) and C = SQRT( (L+M)*(L-M) )            C
C                                                                    C
C____________________________________________________________________C
C Input Variables :-                                                 C
C ===============                                                    C
C  Integer                                                           C
C  -------                                                           C
C     M         : Well it's M isn't it!                              C
C     L         : Well it's L isn't it!                              C
C  Double Precision                                                  C
C  ----------------                                                  C
C     X         : Cos(theta)                                         C
C     PLMIN1    : P_(l-1)^m ( X )                                    C
C     PLMIN2    : P_(l-2)^m ( X )                                    C
C____________________________________________________________________C
C Local Variables :-                                                 C
C ===============                                                    C
C  Integer                                                           C
C  -------                                                           C
C     RL        : DBLE  ( L )                                        C
C     RM        : DBLE  ( M )                                        C
C____________________________________________________________________C
C
C*********************************************************************
      FUNCTION PLM ( L, M, X, PLMIN1, PLMIN2 )
      IMPLICIT NONE
      DOUBLE PRECISION PLM
C____________________________________________________________________C
C Variable declarations - Parameters ................................C
      INTEGER M,L
      DOUBLE PRECISION X, PLMIN1, PLMIN2
C____________________________________________________________________C
C Variable declarations - Working variables .........................C
      DOUBLE PRECISION RM,RL
C____________________________________________________________________C
C START OF FUNCTION *************************************************C
C____________________________________________________________________C
C
      IF ( L.LT.2 ) THEN
         PRINT *,' You are trying to run function PLM with L < 2.'
         PRINT *,' Program aborted.'
         STOP
      ENDIF
      IF ( M.LT.0 .OR. M.GT.L ) THEN
         PRINT *,' M is out of range in function PLM.'
         PRINT *,' Program aborted.'
         STOP
      ENDIF
      IF ( L.EQ.M ) THEN
         PRINT *,' PLM function called with L = M .'
         PRINT *,' Division by zero would follow.'
         PRINT *,' Program aborted.'
         STOP
      ENDIF
      RM = DBLE( M )
      RL = DBLE( L )
      PLM = ( 2.0d0*RL - 1.0d0 )*X*PLMIN1
      PLM = PLM - PLMIN2*DSQRT( ( RL+RM-1.0d0 )*( RL-RM-1.0d0 ) )
      PLM = PLM/DSQRT( (RL+RM)*(RL-RM) )
      RETURN
      END
C*********************************************************************
C*********************************************************************
C function DPMM ******************************************************
C Steve Gibbons 16.4.97                                              C
C____________________________________________________________________C
C Calculates the derivative of P_m^m(theta) according to equation    C
C 185 in my notes i.e.                                               C
C                                                                    C
C d P_m^m(theta)/ d(theta) = sqrt(2.0)*M*cos(theta)/sin(theta) * A   C
C  with A = product_(i=1)^m ( SQRT( (2i-1)/2i ) * sin(theta) )       C
C____________________________________________________________________C
C Input Variables :-                                                 C
C ===============                                                    C
C  Integer                                                           C
C  -------                                                           C
C     M         : Well it's M isn't it!                              C
C  Double Precision                                                  C
C  ----------------                                                  C
C     SINE      : Sin (theta)                                        C
C____________________________________________________________________C
C
C*********************************************************************
      FUNCTION DPMM ( M , C , S )
      IMPLICIT NONE
      DOUBLE PRECISION DPMM
C____________________________________________________________________C
C Variable declarations - Parameters ................................C
      INTEGER M
      DOUBLE PRECISION C, S
C____________________________________________________________________C
C Variable declarations - Working variables .........................C
      INTEGER I
      DOUBLE PRECISION TWOI
C____________________________________________________________________C
C START OF FUNCTION *************************************************C
C____________________________________________________________________C
C
      IF ( M.LT.0 ) THEN
         PRINT *,' Function DPMM. M < 0 error. Stopped.'
         STOP
      ENDIF
      IF ( M.EQ.0 ) THEN
         DPMM = 0.0d0
         RETURN
      ENDIF
C ................................. so M is greater than 0
      DPMM = DSQRT ( 2.0d0 )*DBLE(M)*C/S
      DO I = 1, M
         TWOI = 2.0d0*DBLE(I)
         DPMM = DPMM * S * DSQRT ( (TWOI - 1.0d0 )/TWOI )
      ENDDO
      RETURN
      END
C*********************************************************************
C*********************************************************************
C function DPMM1 *****************************************************
C Steve Gibbons 18.4.97                                              C
C____________________________________________________________________C
C Calculates d(P_(m+1)^m/d(theta) by equation 185 in my notes        C
C____________________________________________________________________C
C Input Variables :-                                                 C
C ===============                                                    C
C  Integer                                                           C
C  -------                                                           C
C     M         : Well it's M isn't it!                              C
C  Double Precision                                                  C
C  ----------------                                                  C
C     X         : Cos(theta)                                         C
C     S         : Sin(theta)                                         C
C     PMM0      : P_m^m(X) as evaluated by Function PMM.             C
C     DPMM0     : d(P_m^m(X))/d(theta) as evaluated by Function DPMM C
C____________________________________________________________________C
C
C*********************************************************************
      FUNCTION DPMM1 ( M, X, S, PMM0, DPMM0)
      IMPLICIT NONE
      DOUBLE PRECISION DPMM1
C____________________________________________________________________C
C Variable declarations - Parameters ................................C
      INTEGER M
      DOUBLE PRECISION X, S, PMM0, DPMM0
C____________________________________________________________________C
C Variable declarations - Working variables .........................C
      DOUBLE PRECISION RM
C____________________________________________________________________C
C START OF FUNCTION *************************************************C
C____________________________________________________________________C
C
      IF ( M.LT.0 ) THEN
         PRINT *,' Function DPMM1. M < 0 error. Stopped.'
         STOP
      ENDIF
      RM = DBLE( M )
      DPMM1 = DSQRT( 2.0d0*RM+1.0d0 )*( X*DPMM0 - S*PMM0 )
      RETURN
      END
C*********************************************************************
C*********************************************************************
C function DPLM *******************************************************
C Steve Gibbons 18.4.97                                              C
C____________________________________________________________________C
C                                                                    C
C Calculates general P_l^m derivative from eqn 187 in my notes       C
C____________________________________________________________________C
C Input Variables :-                                                 C
C ===============                                                    C
C  Integer                                                           C
C  -------                                                           C
C     M         : Well it's M isn't it!                              C
C     L         : Well it's L isn't it!                              C
C  Double Precision                                                  C
C  ----------------                                                  C
C     X         : Cos(theta)                                         C
C     S         : Sin(theta)                                         C
C     PLMIN1    : P_(l-1)^m ( X )                                    C
C     DPLMN1   : P_(l-1)^m ( X ) derivative                         C
C     DPLMN2   : P_(l-2)^m ( X ) derivative                         C
C____________________________________________________________________C
C Local Variables :-                                                 C
C ===============                                                    C
C  Integer                                                           C
C  -------                                                           C
C     RL        : DBLE  ( L )                                        C
C     RM        : DBLE  ( M )                                        C
C____________________________________________________________________C
C
C*********************************************************************
      FUNCTION DPLM ( L, M, X, S, PLMIN1, DPLMN1, DPLMN2 )
      IMPLICIT NONE
      DOUBLE PRECISION DPLM
C____________________________________________________________________C
C Variable declarations - Parameters ................................C
      INTEGER M,L
      DOUBLE PRECISION X, S, PLMIN1, DPLMN1, DPLMN2
C____________________________________________________________________C
C Variable declarations - Working variables .........................C
      DOUBLE PRECISION RM,RL
C____________________________________________________________________C
C START OF FUNCTION *************************************************C
C____________________________________________________________________C
C
      IF ( L.LT.2 ) THEN
         PRINT *,' You are trying to run function DPLM with L < 2.'
         PRINT *,' Program aborted.'
         STOP
      ENDIF
      IF ( M.LT.0 .OR. M.GT.L ) THEN
         PRINT *,' M is out of range in function DPLM.'
         PRINT *,' Program aborted.'
         STOP
      ENDIF
      IF ( L.EQ.M ) THEN
         PRINT *,' DPLM function called with L = M .'
         PRINT *,' Division by zero would follow.'
         PRINT *,' Program aborted.'
         STOP
      ENDIF
      RM = DBLE( M )
      RL = DBLE( L )
      DPLM = ( 2.0d0*RL - 1.0d0 )*(X*DPLMN1-S*PLMIN1)
      DPLM = DPLM - DPLMN2*DSQRT( ( RL+RM-1.0d0 )*( RL-RM-1.0d0 ) )
      DPLM = DPLM/DSQRT( (RL+RM)*(RL-RM) )
      RETURN
      END
C*********************************************************************
C*********************************************************************
C function SQuare Root of L*(L+1) ************************************
C          --     -       -  - -  ************************************
C Steve Gibbons 25.4.97                                              C
C____________________________________________________________________C
C
      FUNCTION SQRLL1 ( L )
      IMPLICIT NONE
      DOUBLE PRECISION SQRLL1,Q
      INTEGER L
   
      IF ( L.LT.0 ) THEN
         PRINT *,' Function SQRLL1. L less than 0. Program aborted.'
         STOP
      ENDIF
      Q = DBLE( L )
      SQRLL1 = DSQRT( Q*Q + Q )
      RETURN
      END
C*********************************************************************
C*********************************************************************
C integer function INDex FUNction ************************************
C                  ---   ---      ************************************
C Steve Gibbons Thu Sep 16 11:13:38 BST 1999                         C
C____________________________________________________________________C
C                                                                    C
C Returns the location in the solution vector of the value of the    C
C IR^{th} grid node of the IH^{th} harmonic.                         C
C                                                                    C
C                                                                    C
C                                                                    C
C____________________________________________________________________C
C Input variables :-                                                 C
C ===============                                                    C
C  Integer                                                           C
C  -------                                                           C
C     IR        : Number of radial grid node.                        C
C     IH        : Number of harmonic.                                C
C     INARR     : Array of integers. Dimension ( * )                 C
C                 At present the elements of INARR are as follows.   C
C                                                                    C
C                 INARR( 1 ) = IFORMF. Format flag.                  C
C                  This flag decides how the elements are arranged   C
C                  in the solution vector.                           C
C                                                                    C
C                  The current options are:-                         C
C                                                                    C
C                   IFORMF = 1/3. INDFUN = ( IR - 1 )*NH + IH        C
C                   IFORMF = 2/4. INDFUN = ( IH - 1 )*NR + IR        C
C                                                                    C
C  Here, NH is the TOTAL number of harmonics in the solution         C
C  vector - or atleast the part of which is visible to that part     C
C  of the program. NR is the number of radial grid nodes             C
C  corresponding to that harmonic (which for cases IFORMF = 1 and    C
C  IFORMF = 2 is identical for all harmonics - more complicated      C
C  options (for example magnetic fields which have to be resolved    C
C  beyond the region of fluid flow) may be added later and this      C
C  routine should be flexible to all possibilities with extra        C
C  constraints being added in other elements of INARR.               C
C                                                                    C
C  IFORMF = 1/3 is the option likely to be used for the purpose of   C
C  solution as it allows the banded formation of a matrix.           C
C                                                                    C
C  IFORMF = 2/4 is the option likely to be used for the purpose of   C
C  displaying solutions as it stores adjacent nodes for each         C
C  harmonic together.                                                C
C                                                                    C
C                 INARR( 2 ) = NR. See above.                        C
C                 INARR( 3 ) = NH. See above.                        C
C                                                                    C
C____________________________________________________________________C
C
C*********************************************************************
      FUNCTION INDFUN ( IR, IH, INARR )
      IMPLICIT NONE
C____________________________________________________________________C
C Variable declarations - Parameters ................................C
      INTEGER INDFUN, IR, IH, INARR( * )
C____________________________________________________________________C
C Variable declarations - Working variables .........................C
      INTEGER IFORMF, NR, NH
C____________________________________________________________________C
C START OF PROGRAM **************************************************C
C____________________________________________________________________C
C
      IFORMF = INARR( 1 )
      NR     = INARR( 2 )
      NH     = INARR( 3 )
C
      IF ( IR.LT.1 .OR. IR.GT.NR ) THEN
        PRINT *,' Function INDFUN. IR = ', IR
        PRINT *,' NR = ', NR,'. Program aborted.'
        STOP
      ENDIF
C
      IF ( IH.LT.1 .OR. IH.GT.NH ) THEN
        PRINT *,' Function INDFUN. IH = ', IH
        PRINT *,' NH = ', NH,'. Program aborted.'
        STOP
      ENDIF
C
      IF ( IFORMF.EQ.1 .OR. IFORMF.EQ.3 ) THEN
        INDFUN = ( IR - 1 )*NH + IH
        RETURN
      ENDIF
C
      IF ( IFORMF.EQ.2 .OR. IFORMF.EQ.4 ) THEN
        INDFUN = ( IH - 1 )*NR + IR
        RETURN
      ENDIF
C
      PRINT *,' Function INDFUN. IFORMF = ', IFORMF
      PRINT *,' Not current option. Program aborted.'
      STOP
      END
C*********************************************************************
      integer function idamax(n,dx,incx)
c
c     finds the index of element having max. absolute value.
c     jack dongarra, linpack, 3/11/78.
c     modified 3/93 to return if incx .le. 0.
c     modified 12/3/93, array(1) declarations changed to array(*)
c
      double precision dx(*),dmax
      integer i,incx,ix,n
c
      idamax = 0
      if( n.lt.1 .or. incx.le.0 ) return
      idamax = 1
      if(n.eq.1)return
      if(incx.eq.1)go to 20
c
c        code for increment not equal to 1
c
      ix = 1
      dmax = dabs(dx(1))
      ix = ix + incx
      do 10 i = 2,n
         if(dabs(dx(ix)).le.dmax) go to 5
         idamax = i
         dmax = dabs(dx(ix))
    5    ix = ix + incx
   10 continue
      return
c
c        code for increment equal to 1
c
   20 dmax = dabs(dx(1))
      do 30 i = 2,n
         if(dabs(dx(i)).le.dmax) go to 30
         idamax = i
         dmax = dabs(dx(i))
   30 continue
      return
      end
      SUBROUTINE DGEMM ( TRANSA, TRANSB, M, N, K, ALPHA, A, LDA, B, LDB,
     $                   BETA, C, LDC )
*     .. Scalar Arguments ..
      CHARACTER*1        TRANSA, TRANSB
      INTEGER            M, N, K, LDA, LDB, LDC
      DOUBLE PRECISION   ALPHA, BETA
*     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), B( LDB, * ), C( LDC, * )
*     ..
*
*  Purpose
*  =======
*
*  DGEMM  performs one of the matrix-matrix operations
*
*     C := alpha*op( A )*op( B ) + beta*C,
*
*  where  op( X ) is one of
*
*     op( X ) = X   or   op( X ) = X',
*
*  alpha and beta are scalars, and A, B and C are matrices, with op( A )
*  an m by k matrix,  op( B )  a  k by n matrix and  C an m by n matrix.
*
*  Parameters
*  ==========
*
*  TRANSA - CHARACTER*1.
*           On entry, TRANSA specifies the form of op( A ) to be used in
*           the matrix multiplication as follows:
*
*              TRANSA = 'N' or 'n',  op( A ) = A.
*
*              TRANSA = 'T' or 't',  op( A ) = A'.
*
*              TRANSA = 'C' or 'c',  op( A ) = A'.
*
*           Unchanged on exit.
*
*  TRANSB - CHARACTER*1.
*           On entry, TRANSB specifies the form of op( B ) to be used in
*           the matrix multiplication as follows:
*
*              TRANSB = 'N' or 'n',  op( B ) = B.
*
*              TRANSB = 'T' or 't',  op( B ) = B'.
*
*              TRANSB = 'C' or 'c',  op( B ) = B'.
*
*           Unchanged on exit.
*
*  M      - INTEGER.
*           On entry,  M  specifies  the number  of rows  of the  matrix
*           op( A )  and of the  matrix  C.  M  must  be at least  zero.
*           Unchanged on exit.
*
*  N      - INTEGER.
*           On entry,  N  specifies the number  of columns of the matrix
*           op( B ) and the number of columns of the matrix C. N must be
*           at least zero.
*           Unchanged on exit.
*
*  K      - INTEGER.
*           On entry,  K  specifies  the number of columns of the matrix
*           op( A ) and the number of rows of the matrix op( B ). K must
*           be at least  zero.
*           Unchanged on exit.
*
*  ALPHA  - DOUBLE PRECISION.
*           On entry, ALPHA specifies the scalar alpha.
*           Unchanged on exit.
*
*  A      - DOUBLE PRECISION array of DIMENSION ( LDA, ka ), where ka is
*           k  when  TRANSA = 'N' or 'n',  and is  m  otherwise.
*           Before entry with  TRANSA = 'N' or 'n',  the leading  m by k
*           part of the array  A  must contain the matrix  A,  otherwise
*           the leading  k by m  part of the array  A  must contain  the
*           matrix A.
*           Unchanged on exit.
*
*  LDA    - INTEGER.
*           On entry, LDA specifies the first dimension of A as declared
*           in the calling (sub) program. When  TRANSA = 'N' or 'n' then
*           LDA must be at least  max( 1, m ), otherwise  LDA must be at
*           least  max( 1, k ).
*           Unchanged on exit.
*
*  B      - DOUBLE PRECISION array of DIMENSION ( LDB, kb ), where kb is
*           n  when  TRANSB = 'N' or 'n',  and is  k  otherwise.
*           Before entry with  TRANSB = 'N' or 'n',  the leading  k by n
*           part of the array  B  must contain the matrix  B,  otherwise
*           the leading  n by k  part of the array  B  must contain  the
*           matrix B.
*           Unchanged on exit.
*
*  LDB    - INTEGER.
*           On entry, LDB specifies the first dimension of B as declared
*           in the calling (sub) program. When  TRANSB = 'N' or 'n' then
*           LDB must be at least  max( 1, k ), otherwise  LDB must be at
*           least  max( 1, n ).
*           Unchanged on exit.
*
*  BETA   - DOUBLE PRECISION.
*           On entry,  BETA  specifies the scalar  beta.  When  BETA  is
*           supplied as zero then C need not be set on input.
*           Unchanged on exit.
*
*  C      - DOUBLE PRECISION array of DIMENSION ( LDC, n ).
*           Before entry, the leading  m by n  part of the array  C must
*           contain the matrix  C,  except when  beta  is zero, in which
*           case C need not be set on entry.
*           On exit, the array  C  is overwritten by the  m by n  matrix
*           ( alpha*op( A )*op( B ) + beta*C ).
*
*  LDC    - INTEGER.
*           On entry, LDC specifies the first dimension of C as declared
*           in  the  calling  (sub)  program.   LDC  must  be  at  least
*           max( 1, m ).
*           Unchanged on exit.
*
*
*  Level 3 Blas routine.
*
*  -- Written on 8-February-1989.
*     Jack Dongarra, Argonne National Laboratory.
*     Iain Duff, AERE Harwell.
*     Jeremy Du Croz, Numerical Algorithms Group Ltd.
*     Sven Hammarling, Numerical Algorithms Group Ltd.
*
*
*     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
*     .. External Subroutines ..
      EXTERNAL           XERBLA
*     .. Intrinsic Functions ..
      INTRINSIC          MAX
*     .. Local Scalars ..
      LOGICAL            NOTA, NOTB
      INTEGER            I, INFO, J, L, NCOLA, NROWA, NROWB
      DOUBLE PRECISION   TEMP
*     .. Parameters ..
      DOUBLE PRECISION   ONE         , ZERO
      PARAMETER        ( ONE = 1.0D+0, ZERO = 0.0D+0 )
*     ..
*     .. Executable Statements ..
*
*     Set  NOTA  and  NOTB  as  true if  A  and  B  respectively are not
*     transposed and set  NROWA, NCOLA and  NROWB  as the number of rows
*     and  columns of  A  and the  number of  rows  of  B  respectively.
*
      NOTA  = LSAME( TRANSA, 'N' )
      NOTB  = LSAME( TRANSB, 'N' )
      IF( NOTA )THEN
         NROWA = M
         NCOLA = K
      ELSE
         NROWA = K
         NCOLA = M
      END IF
      IF( NOTB )THEN
         NROWB = K
      ELSE
         NROWB = N
      END IF
*
*     Test the input parameters.
*
      INFO = 0
      IF(      ( .NOT.NOTA                 ).AND.
     $         ( .NOT.LSAME( TRANSA, 'C' ) ).AND.
     $         ( .NOT.LSAME( TRANSA, 'T' ) )      )THEN
         INFO = 1
      ELSE IF( ( .NOT.NOTB                 ).AND.
     $         ( .NOT.LSAME( TRANSB, 'C' ) ).AND.
     $         ( .NOT.LSAME( TRANSB, 'T' ) )      )THEN
         INFO = 2
      ELSE IF( M  .LT.0               )THEN
         INFO = 3
      ELSE IF( N  .LT.0               )THEN
         INFO = 4
      ELSE IF( K  .LT.0               )THEN
         INFO = 5
      ELSE IF( LDA.LT.MAX( 1, NROWA ) )THEN
         INFO = 8
      ELSE IF( LDB.LT.MAX( 1, NROWB ) )THEN
         INFO = 10
      ELSE IF( LDC.LT.MAX( 1, M     ) )THEN
         INFO = 13
      END IF
      IF( INFO.NE.0 )THEN
         CALL XERBLA( 'DGEMM ', INFO )
         RETURN
      END IF
*
*     Quick return if possible.
*
      IF( ( M.EQ.0 ).OR.( N.EQ.0 ).OR.
     $    ( ( ( ALPHA.EQ.ZERO ).OR.( K.EQ.0 ) ).AND.( BETA.EQ.ONE ) ) )
     $   RETURN
*
*     And if  alpha.eq.zero.
*
      IF( ALPHA.EQ.ZERO )THEN
         IF( BETA.EQ.ZERO )THEN
            DO 20, J = 1, N
               DO 10, I = 1, M
                  C( I, J ) = ZERO
   10          CONTINUE
   20       CONTINUE
         ELSE
            DO 40, J = 1, N
               DO 30, I = 1, M
                  C( I, J ) = BETA*C( I, J )
   30          CONTINUE
   40       CONTINUE
         END IF
         RETURN
      END IF
*
*     Start the operations.
*
      IF( NOTB )THEN
         IF( NOTA )THEN
*
*           Form  C := alpha*A*B + beta*C.
*
            DO 90, J = 1, N
               IF( BETA.EQ.ZERO )THEN
                  DO 50, I = 1, M
                     C( I, J ) = ZERO
   50             CONTINUE
               ELSE IF( BETA.NE.ONE )THEN
                  DO 60, I = 1, M
                     C( I, J ) = BETA*C( I, J )
   60             CONTINUE
               END IF
               DO 80, L = 1, K
                  IF( B( L, J ).NE.ZERO )THEN
                     TEMP = ALPHA*B( L, J )
                     DO 70, I = 1, M
                        C( I, J ) = C( I, J ) + TEMP*A( I, L )
   70                CONTINUE
                  END IF
   80          CONTINUE
   90       CONTINUE
         ELSE
*
*           Form  C := alpha*A'*B + beta*C
*
            DO 120, J = 1, N
               DO 110, I = 1, M
                  TEMP = ZERO
                  DO 100, L = 1, K
                     TEMP = TEMP + A( L, I )*B( L, J )
  100             CONTINUE
                  IF( BETA.EQ.ZERO )THEN
                     C( I, J ) = ALPHA*TEMP
                  ELSE
                     C( I, J ) = ALPHA*TEMP + BETA*C( I, J )
                  END IF
  110          CONTINUE
  120       CONTINUE
         END IF
      ELSE
         IF( NOTA )THEN
*
*           Form  C := alpha*A*B' + beta*C
*
            DO 170, J = 1, N
               IF( BETA.EQ.ZERO )THEN
                  DO 130, I = 1, M
                     C( I, J ) = ZERO
  130             CONTINUE
               ELSE IF( BETA.NE.ONE )THEN
                  DO 140, I = 1, M
                     C( I, J ) = BETA*C( I, J )
  140             CONTINUE
               END IF
               DO 160, L = 1, K
                  IF( B( J, L ).NE.ZERO )THEN
                     TEMP = ALPHA*B( J, L )
                     DO 150, I = 1, M
                        C( I, J ) = C( I, J ) + TEMP*A( I, L )
  150                CONTINUE
                  END IF
  160          CONTINUE
  170       CONTINUE
         ELSE
*
*           Form  C := alpha*A'*B' + beta*C
*
            DO 200, J = 1, N
               DO 190, I = 1, M
                  TEMP = ZERO
                  DO 180, L = 1, K
                     TEMP = TEMP + A( L, I )*B( J, L )
  180             CONTINUE
                  IF( BETA.EQ.ZERO )THEN
                     C( I, J ) = ALPHA*TEMP
                  ELSE
                     C( I, J ) = ALPHA*TEMP + BETA*C( I, J )
                  END IF
  190          CONTINUE
  200       CONTINUE
         END IF
      END IF
*
      RETURN
*
*     End of DGEMM .
*
      END
      SUBROUTINE DTRSM ( SIDE, UPLO, TRANSA, DIAG, M, N, ALPHA, A, LDA,
     $                   B, LDB )
*     .. Scalar Arguments ..
      CHARACTER*1        SIDE, UPLO, TRANSA, DIAG
      INTEGER            M, N, LDA, LDB
      DOUBLE PRECISION   ALPHA
*     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), B( LDB, * )
*     ..
*
*  Purpose
*  =======
*
*  DTRSM  solves one of the matrix equations
*
*     op( A )*X = alpha*B,   or   X*op( A ) = alpha*B,
*
*  where alpha is a scalar, X and B are m by n matrices, A is a unit, or
*  non-unit,  upper or lower triangular matrix  and  op( A )  is one  of
*
*     op( A ) = A   or   op( A ) = A'.
*
*  The matrix X is overwritten on B.
*
*  Parameters
*  ==========
*
*  SIDE   - CHARACTER*1.
*           On entry, SIDE specifies whether op( A ) appears on the left
*           or right of X as follows:
*
*              SIDE = 'L' or 'l'   op( A )*X = alpha*B.
*
*              SIDE = 'R' or 'r'   X*op( A ) = alpha*B.
*
*           Unchanged on exit.
*
*  UPLO   - CHARACTER*1.
*           On entry, UPLO specifies whether the matrix A is an upper or
*           lower triangular matrix as follows:
*
*              UPLO = 'U' or 'u'   A is an upper triangular matrix.
*
*              UPLO = 'L' or 'l'   A is a lower triangular matrix.
*
*           Unchanged on exit.
*
*  TRANSA - CHARACTER*1.
*           On entry, TRANSA specifies the form of op( A ) to be used in
*           the matrix multiplication as follows:
*
*              TRANSA = 'N' or 'n'   op( A ) = A.
*
*              TRANSA = 'T' or 't'   op( A ) = A'.
*
*              TRANSA = 'C' or 'c'   op( A ) = A'.
*
*           Unchanged on exit.
*
*  DIAG   - CHARACTER*1.
*           On entry, DIAG specifies whether or not A is unit triangular
*           as follows:
*
*              DIAG = 'U' or 'u'   A is assumed to be unit triangular.
*
*              DIAG = 'N' or 'n'   A is not assumed to be unit
*                                  triangular.
*
*           Unchanged on exit.
*
*  M      - INTEGER.
*           On entry, M specifies the number of rows of B. M must be at
*           least zero.
*           Unchanged on exit.
*
*  N      - INTEGER.
*           On entry, N specifies the number of columns of B.  N must be
*           at least zero.
*           Unchanged on exit.
*
*  ALPHA  - DOUBLE PRECISION.
*           On entry,  ALPHA specifies the scalar  alpha. When  alpha is
*           zero then  A is not referenced and  B need not be set before
*           entry.
*           Unchanged on exit.
*
*  A      - DOUBLE PRECISION array of DIMENSION ( LDA, k ), where k is m
*           when  SIDE = 'L' or 'l'  and is  n  when  SIDE = 'R' or 'r'.
*           Before entry  with  UPLO = 'U' or 'u',  the  leading  k by k
*           upper triangular part of the array  A must contain the upper
*           triangular matrix  and the strictly lower triangular part of
*           A is not referenced.
*           Before entry  with  UPLO = 'L' or 'l',  the  leading  k by k
*           lower triangular part of the array  A must contain the lower
*           triangular matrix  and the strictly upper triangular part of
*           A is not referenced.
*           Note that when  DIAG = 'U' or 'u',  the diagonal elements of
*           A  are not referenced either,  but are assumed to be  unity.
*           Unchanged on exit.
*
*  LDA    - INTEGER.
*           On entry, LDA specifies the first dimension of A as declared
*           in the calling (sub) program.  When  SIDE = 'L' or 'l'  then
*           LDA  must be at least  max( 1, m ),  when  SIDE = 'R' or 'r'
*           then LDA must be at least max( 1, n ).
*           Unchanged on exit.
*
*  B      - DOUBLE PRECISION array of DIMENSION ( LDB, n ).
*           Before entry,  the leading  m by n part of the array  B must
*           contain  the  right-hand  side  matrix  B,  and  on exit  is
*           overwritten by the solution matrix  X.
*
*  LDB    - INTEGER.
*           On entry, LDB specifies the first dimension of B as declared
*           in  the  calling  (sub)  program.   LDB  must  be  at  least
*           max( 1, m ).
*           Unchanged on exit.
*
*
*  Level 3 Blas routine.
*
*
*  -- Written on 8-February-1989.
*     Jack Dongarra, Argonne National Laboratory.
*     Iain Duff, AERE Harwell.
*     Jeremy Du Croz, Numerical Algorithms Group Ltd.
*     Sven Hammarling, Numerical Algorithms Group Ltd.
*
*
*     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
*     .. External Subroutines ..
      EXTERNAL           XERBLA
*     .. Intrinsic Functions ..
      INTRINSIC          MAX
*     .. Local Scalars ..
      LOGICAL            LSIDE, NOUNIT, UPPER
      INTEGER            I, INFO, J, K, NROWA
      DOUBLE PRECISION   TEMP
*     .. Parameters ..
      DOUBLE PRECISION   ONE         , ZERO
      PARAMETER        ( ONE = 1.0D+0, ZERO = 0.0D+0 )
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters.
*
      LSIDE  = LSAME( SIDE  , 'L' )
      IF( LSIDE )THEN
         NROWA = M
      ELSE
         NROWA = N
      END IF
      NOUNIT = LSAME( DIAG  , 'N' )
      UPPER  = LSAME( UPLO  , 'U' )
*
      INFO   = 0
      IF(      ( .NOT.LSIDE                ).AND.
     $         ( .NOT.LSAME( SIDE  , 'R' ) )      )THEN
         INFO = 1
      ELSE IF( ( .NOT.UPPER                ).AND.
     $         ( .NOT.LSAME( UPLO  , 'L' ) )      )THEN
         INFO = 2
      ELSE IF( ( .NOT.LSAME( TRANSA, 'N' ) ).AND.
     $         ( .NOT.LSAME( TRANSA, 'T' ) ).AND.
     $         ( .NOT.LSAME( TRANSA, 'C' ) )      )THEN
         INFO = 3
      ELSE IF( ( .NOT.LSAME( DIAG  , 'U' ) ).AND.
     $         ( .NOT.LSAME( DIAG  , 'N' ) )      )THEN
         INFO = 4
      ELSE IF( M  .LT.0               )THEN
         INFO = 5
      ELSE IF( N  .LT.0               )THEN
         INFO = 6
      ELSE IF( LDA.LT.MAX( 1, NROWA ) )THEN
         INFO = 9
      ELSE IF( LDB.LT.MAX( 1, M     ) )THEN
         INFO = 11
      END IF
      IF( INFO.NE.0 )THEN
         CALL XERBLA( 'DTRSM ', INFO )
         RETURN
      END IF
*
*     Quick return if possible.
*
      IF( N.EQ.0 )
     $   RETURN
*
*     And when  alpha.eq.zero.
*
      IF( ALPHA.EQ.ZERO )THEN
         DO 20, J = 1, N
            DO 10, I = 1, M
               B( I, J ) = ZERO
   10       CONTINUE
   20    CONTINUE
         RETURN
      END IF
*
*     Start the operations.
*
      IF( LSIDE )THEN
         IF( LSAME( TRANSA, 'N' ) )THEN
*
*           Form  B := alpha*inv( A )*B.
*
            IF( UPPER )THEN
               DO 60, J = 1, N
                  IF( ALPHA.NE.ONE )THEN
                     DO 30, I = 1, M
                        B( I, J ) = ALPHA*B( I, J )
   30                CONTINUE
                  END IF
                  DO 50, K = M, 1, -1
                     IF( B( K, J ).NE.ZERO )THEN
                        IF( NOUNIT )
     $                     B( K, J ) = B( K, J )/A( K, K )
                        DO 40, I = 1, K - 1
                           B( I, J ) = B( I, J ) - B( K, J )*A( I, K )
   40                   CONTINUE
                     END IF
   50             CONTINUE
   60          CONTINUE
            ELSE
               DO 100, J = 1, N
                  IF( ALPHA.NE.ONE )THEN
                     DO 70, I = 1, M
                        B( I, J ) = ALPHA*B( I, J )
   70                CONTINUE
                  END IF
                  DO 90 K = 1, M
                     IF( B( K, J ).NE.ZERO )THEN
                        IF( NOUNIT )
     $                     B( K, J ) = B( K, J )/A( K, K )
                        DO 80, I = K + 1, M
                           B( I, J ) = B( I, J ) - B( K, J )*A( I, K )
   80                   CONTINUE
                     END IF
   90             CONTINUE
  100          CONTINUE
            END IF
         ELSE
*
*           Form  B := alpha*inv( A' )*B.
*
            IF( UPPER )THEN
               DO 130, J = 1, N
                  DO 120, I = 1, M
                     TEMP = ALPHA*B( I, J )
                     DO 110, K = 1, I - 1
                        TEMP = TEMP - A( K, I )*B( K, J )
  110                CONTINUE
                     IF( NOUNIT )
     $                  TEMP = TEMP/A( I, I )
                     B( I, J ) = TEMP
  120             CONTINUE
  130          CONTINUE
            ELSE
               DO 160, J = 1, N
                  DO 150, I = M, 1, -1
                     TEMP = ALPHA*B( I, J )
                     DO 140, K = I + 1, M
                        TEMP = TEMP - A( K, I )*B( K, J )
  140                CONTINUE
                     IF( NOUNIT )
     $                  TEMP = TEMP/A( I, I )
                     B( I, J ) = TEMP
  150             CONTINUE
  160          CONTINUE
            END IF
         END IF
      ELSE
         IF( LSAME( TRANSA, 'N' ) )THEN
*
*           Form  B := alpha*B*inv( A ).
*
            IF( UPPER )THEN
               DO 210, J = 1, N
                  IF( ALPHA.NE.ONE )THEN
                     DO 170, I = 1, M
                        B( I, J ) = ALPHA*B( I, J )
  170                CONTINUE
                  END IF
                  DO 190, K = 1, J - 1
                     IF( A( K, J ).NE.ZERO )THEN
                        DO 180, I = 1, M
                           B( I, J ) = B( I, J ) - A( K, J )*B( I, K )
  180                   CONTINUE
                     END IF
  190             CONTINUE
                  IF( NOUNIT )THEN
                     TEMP = ONE/A( J, J )
                     DO 200, I = 1, M
                        B( I, J ) = TEMP*B( I, J )
  200                CONTINUE
                  END IF
  210          CONTINUE
            ELSE
               DO 260, J = N, 1, -1
                  IF( ALPHA.NE.ONE )THEN
                     DO 220, I = 1, M
                        B( I, J ) = ALPHA*B( I, J )
  220                CONTINUE
                  END IF
                  DO 240, K = J + 1, N
                     IF( A( K, J ).NE.ZERO )THEN
                        DO 230, I = 1, M
                           B( I, J ) = B( I, J ) - A( K, J )*B( I, K )
  230                   CONTINUE
                     END IF
  240             CONTINUE
                  IF( NOUNIT )THEN
                     TEMP = ONE/A( J, J )
                     DO 250, I = 1, M
                       B( I, J ) = TEMP*B( I, J )
  250                CONTINUE
                  END IF
  260          CONTINUE
            END IF
         ELSE
*
*           Form  B := alpha*B*inv( A' ).
*
            IF( UPPER )THEN
               DO 310, K = N, 1, -1
                  IF( NOUNIT )THEN
                     TEMP = ONE/A( K, K )
                     DO 270, I = 1, M
                        B( I, K ) = TEMP*B( I, K )
  270                CONTINUE
                  END IF
                  DO 290, J = 1, K - 1
                     IF( A( J, K ).NE.ZERO )THEN
                        TEMP = A( J, K )
                        DO 280, I = 1, M
                           B( I, J ) = B( I, J ) - TEMP*B( I, K )
  280                   CONTINUE
                     END IF
  290             CONTINUE
                  IF( ALPHA.NE.ONE )THEN
                     DO 300, I = 1, M
                        B( I, K ) = ALPHA*B( I, K )
  300                CONTINUE
                  END IF
  310          CONTINUE
            ELSE
               DO 360, K = 1, N
                  IF( NOUNIT )THEN
                     TEMP = ONE/A( K, K )
                     DO 320, I = 1, M
                        B( I, K ) = TEMP*B( I, K )
  320                CONTINUE
                  END IF
                  DO 340, J = K + 1, N
                     IF( A( J, K ).NE.ZERO )THEN
                        TEMP = A( J, K )
                        DO 330, I = 1, M
                           B( I, J ) = B( I, J ) - TEMP*B( I, K )
  330                   CONTINUE
                     END IF
  340             CONTINUE
                  IF( ALPHA.NE.ONE )THEN
                     DO 350, I = 1, M
                        B( I, K ) = ALPHA*B( I, K )
  350                CONTINUE
                  END IF
  360          CONTINUE
            END IF
         END IF
      END IF
*
      RETURN
*
*     End of DTRSM .
*
      END
      SUBROUTINE DGEMV ( TRANS, M, N, ALPHA, A, LDA, X, INCX,
     $                   BETA, Y, INCY )
*     .. Scalar Arguments ..
      DOUBLE PRECISION   ALPHA, BETA
      INTEGER            INCX, INCY, LDA, M, N
      CHARACTER*1        TRANS
*     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), X( * ), Y( * )
*     ..
*
*  Purpose
*  =======
*
*  DGEMV  performs one of the matrix-vector operations
*
*     y := alpha*A*x + beta*y,   or   y := alpha*A'*x + beta*y,
*
*  where alpha and beta are scalars, x and y are vectors and A is an
*  m by n matrix.
*
*  Parameters
*  ==========
*
*  TRANS  - CHARACTER*1.
*           On entry, TRANS specifies the operation to be performed as
*           follows:
*
*              TRANS = 'N' or 'n'   y := alpha*A*x + beta*y.
*
*              TRANS = 'T' or 't'   y := alpha*A'*x + beta*y.
*
*              TRANS = 'C' or 'c'   y := alpha*A'*x + beta*y.
*
*           Unchanged on exit.
*
*  M      - INTEGER.
*           On entry, M specifies the number of rows of the matrix A.
*           M must be at least zero.
*           Unchanged on exit.
*
*  N      - INTEGER.
*           On entry, N specifies the number of columns of the matrix A.
*           N must be at least zero.
*           Unchanged on exit.
*
*  ALPHA  - DOUBLE PRECISION.
*           On entry, ALPHA specifies the scalar alpha.
*           Unchanged on exit.
*
*  A      - DOUBLE PRECISION array of DIMENSION ( LDA, n ).
*           Before entry, the leading m by n part of the array A must
*           contain the matrix of coefficients.
*           Unchanged on exit.
*
*  LDA    - INTEGER.
*           On entry, LDA specifies the first dimension of A as declared
*           in the calling (sub) program. LDA must be at least
*           max( 1, m ).
*           Unchanged on exit.
*
*  X      - DOUBLE PRECISION array of DIMENSION at least
*           ( 1 + ( n - 1 )*abs( INCX ) ) when TRANS = 'N' or 'n'
*           and at least
*           ( 1 + ( m - 1 )*abs( INCX ) ) otherwise.
*           Before entry, the incremented array X must contain the
*           vector x.
*           Unchanged on exit.
*
*  INCX   - INTEGER.
*           On entry, INCX specifies the increment for the elements of
*           X. INCX must not be zero.
*           Unchanged on exit.
*
*  BETA   - DOUBLE PRECISION.
*           On entry, BETA specifies the scalar beta. When BETA is
*           supplied as zero then Y need not be set on input.
*           Unchanged on exit.
*
*  Y      - DOUBLE PRECISION array of DIMENSION at least
*           ( 1 + ( m - 1 )*abs( INCY ) ) when TRANS = 'N' or 'n'
*           and at least
*           ( 1 + ( n - 1 )*abs( INCY ) ) otherwise.
*           Before entry with BETA non-zero, the incremented array Y
*           must contain the vector y. On exit, Y is overwritten by the
*           updated vector y.
*
*  INCY   - INTEGER.
*           On entry, INCY specifies the increment for the elements of
*           Y. INCY must not be zero.
*           Unchanged on exit.
*
*
*  Level 2 Blas routine.
*
*  -- Written on 22-October-1986.
*     Jack Dongarra, Argonne National Lab.
*     Jeremy Du Croz, Nag Central Office.
*     Sven Hammarling, Nag Central Office.
*     Richard Hanson, Sandia National Labs.
*
*
*     .. Parameters ..
      DOUBLE PRECISION   ONE         , ZERO
      PARAMETER        ( ONE = 1.0D+0, ZERO = 0.0D+0 )
*     .. Local Scalars ..
      DOUBLE PRECISION   TEMP
      INTEGER            I, INFO, IX, IY, J, JX, JY, KX, KY, LENX, LENY
*     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
*     .. External Subroutines ..
      EXTERNAL           XERBLA
*     .. Intrinsic Functions ..
      INTRINSIC          MAX
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters.
*
      INFO = 0
      IF     ( .NOT.LSAME( TRANS, 'N' ).AND.
     $         .NOT.LSAME( TRANS, 'T' ).AND.
     $         .NOT.LSAME( TRANS, 'C' )      )THEN
         INFO = 1
      ELSE IF( M.LT.0 )THEN
         INFO = 2
      ELSE IF( N.LT.0 )THEN
         INFO = 3
      ELSE IF( LDA.LT.MAX( 1, M ) )THEN
         INFO = 6
      ELSE IF( INCX.EQ.0 )THEN
         INFO = 8
      ELSE IF( INCY.EQ.0 )THEN
         INFO = 11
      END IF
      IF( INFO.NE.0 )THEN
         CALL XERBLA( 'DGEMV ', INFO )
         RETURN
      END IF
*
*     Quick return if possible.
*
      IF( ( M.EQ.0 ).OR.( N.EQ.0 ).OR.
     $    ( ( ALPHA.EQ.ZERO ).AND.( BETA.EQ.ONE ) ) )
     $   RETURN
*
*     Set  LENX  and  LENY, the lengths of the vectors x and y, and set
*     up the start points in  X  and  Y.
*
      IF( LSAME( TRANS, 'N' ) )THEN
         LENX = N
         LENY = M
      ELSE
         LENX = M
         LENY = N
      END IF
      IF( INCX.GT.0 )THEN
         KX = 1
      ELSE
         KX = 1 - ( LENX - 1 )*INCX
      END IF
      IF( INCY.GT.0 )THEN
         KY = 1
      ELSE
         KY = 1 - ( LENY - 1 )*INCY
      END IF
*
*     Start the operations. In this version the elements of A are
*     accessed sequentially with one pass through A.
*
*     First form  y := beta*y.
*
      IF( BETA.NE.ONE )THEN
         IF( INCY.EQ.1 )THEN
            IF( BETA.EQ.ZERO )THEN
               DO 10, I = 1, LENY
                  Y( I ) = ZERO
   10          CONTINUE
            ELSE
               DO 20, I = 1, LENY
                  Y( I ) = BETA*Y( I )
   20          CONTINUE
            END IF
         ELSE
            IY = KY
            IF( BETA.EQ.ZERO )THEN
               DO 30, I = 1, LENY
                  Y( IY ) = ZERO
                  IY      = IY   + INCY
   30          CONTINUE
            ELSE
               DO 40, I = 1, LENY
                  Y( IY ) = BETA*Y( IY )
                  IY      = IY           + INCY
   40          CONTINUE
            END IF
         END IF
      END IF
      IF( ALPHA.EQ.ZERO )
     $   RETURN
      IF( LSAME( TRANS, 'N' ) )THEN
*
*        Form  y := alpha*A*x + y.
*
         JX = KX
         IF( INCY.EQ.1 )THEN
            DO 60, J = 1, N
               IF( X( JX ).NE.ZERO )THEN
                  TEMP = ALPHA*X( JX )
                  DO 50, I = 1, M
                     Y( I ) = Y( I ) + TEMP*A( I, J )
   50             CONTINUE
               END IF
               JX = JX + INCX
   60       CONTINUE
         ELSE
            DO 80, J = 1, N
               IF( X( JX ).NE.ZERO )THEN
                  TEMP = ALPHA*X( JX )
                  IY   = KY
                  DO 70, I = 1, M
                     Y( IY ) = Y( IY ) + TEMP*A( I, J )
                     IY      = IY      + INCY
   70             CONTINUE
               END IF
               JX = JX + INCX
   80       CONTINUE
         END IF
      ELSE
*
*        Form  y := alpha*A'*x + y.
*
         JY = KY
         IF( INCX.EQ.1 )THEN
            DO 100, J = 1, N
               TEMP = ZERO
               DO 90, I = 1, M
                  TEMP = TEMP + A( I, J )*X( I )
   90          CONTINUE
               Y( JY ) = Y( JY ) + ALPHA*TEMP
               JY      = JY      + INCY
  100       CONTINUE
         ELSE
            DO 120, J = 1, N
               TEMP = ZERO
               IX   = KX
               DO 110, I = 1, M
                  TEMP = TEMP + A( I, J )*X( IX )
                  IX   = IX   + INCX
  110          CONTINUE
               Y( JY ) = Y( JY ) + ALPHA*TEMP
               JY      = JY      + INCY
  120       CONTINUE
         END IF
      END IF
*
      RETURN
*
*     End of DGEMV .
*
      END
      subroutine  dswap (n,dx,incx,dy,incy)
c
c     interchanges two vectors.
c     uses unrolled loops for increments equal one.
c     jack dongarra, linpack, 3/11/78.
c     modified 12/3/93, array(1) declarations changed to array(*)
c
      double precision dx(*),dy(*),dtemp
      integer i,incx,incy,ix,iy,m,mp1,n
c
      if(n.le.0)return
      if(incx.eq.1.and.incy.eq.1)go to 20
c
c       code for unequal increments or equal increments not equal
c         to 1
c
      ix = 1
      iy = 1
      if(incx.lt.0)ix = (-n+1)*incx + 1
      if(incy.lt.0)iy = (-n+1)*incy + 1
      do 10 i = 1,n
        dtemp = dx(ix)
        dx(ix) = dy(iy)
        dy(iy) = dtemp
        ix = ix + incx
        iy = iy + incy
   10 continue
      return
c
c       code for both increments equal to 1
c
c
c       clean-up loop
c
   20 m = mod(n,3)
      if( m .eq. 0 ) go to 40
      do 30 i = 1,m
        dtemp = dx(i)
        dx(i) = dy(i)
        dy(i) = dtemp
   30 continue
      if( n .lt. 3 ) return
   40 mp1 = m + 1
      do 50 i = mp1,n,3
        dtemp = dx(i)
        dx(i) = dy(i)
        dy(i) = dtemp
        dtemp = dx(i + 1)
        dx(i + 1) = dy(i + 1)
        dy(i + 1) = dtemp
        dtemp = dx(i + 2)
        dx(i + 2) = dy(i + 2)
        dy(i + 2) = dtemp
   50 continue
      return
      end
      SUBROUTINE DGER  ( M, N, ALPHA, X, INCX, Y, INCY, A, LDA )
*     .. Scalar Arguments ..
      DOUBLE PRECISION   ALPHA
      INTEGER            INCX, INCY, LDA, M, N
*     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), X( * ), Y( * )
*     ..
*
*  Purpose
*  =======
*
*  DGER   performs the rank 1 operation
*
*     A := alpha*x*y' + A,
*
*  where alpha is a scalar, x is an m element vector, y is an n element
*  vector and A is an m by n matrix.
*
*  Parameters
*  ==========
*
*  M      - INTEGER.
*           On entry, M specifies the number of rows of the matrix A.
*           M must be at least zero.
*           Unchanged on exit.
*
*  N      - INTEGER.
*           On entry, N specifies the number of columns of the matrix A.
*           N must be at least zero.
*           Unchanged on exit.
*
*  ALPHA  - DOUBLE PRECISION.
*           On entry, ALPHA specifies the scalar alpha.
*           Unchanged on exit.
*
*  X      - DOUBLE PRECISION array of dimension at least
*           ( 1 + ( m - 1 )*abs( INCX ) ).
*           Before entry, the incremented array X must contain the m
*           element vector x.
*           Unchanged on exit.
*
*  INCX   - INTEGER.
*           On entry, INCX specifies the increment for the elements of
*           X. INCX must not be zero.
*           Unchanged on exit.
*
*  Y      - DOUBLE PRECISION array of dimension at least
*           ( 1 + ( n - 1 )*abs( INCY ) ).
*           Before entry, the incremented array Y must contain the n
*           element vector y.
*           Unchanged on exit.
*
*  INCY   - INTEGER.
*           On entry, INCY specifies the increment for the elements of
*           Y. INCY must not be zero.
*           Unchanged on exit.
*
*  A      - DOUBLE PRECISION array of DIMENSION ( LDA, n ).
*           Before entry, the leading m by n part of the array A must
*           contain the matrix of coefficients. On exit, A is
*           overwritten by the updated matrix.
*
*  LDA    - INTEGER.
*           On entry, LDA specifies the first dimension of A as declared
*           in the calling (sub) program. LDA must be at least
*           max( 1, m ).
*           Unchanged on exit.
*
*
*  Level 2 Blas routine.
*
*  -- Written on 22-October-1986.
*     Jack Dongarra, Argonne National Lab.
*     Jeremy Du Croz, Nag Central Office.
*     Sven Hammarling, Nag Central Office.
*     Richard Hanson, Sandia National Labs.
*
*
*     .. Parameters ..
      DOUBLE PRECISION   ZERO
      PARAMETER        ( ZERO = 0.0D+0 )
*     .. Local Scalars ..
      DOUBLE PRECISION   TEMP
      INTEGER            I, INFO, IX, J, JY, KX
*     .. External Subroutines ..
      EXTERNAL           XERBLA
*     .. Intrinsic Functions ..
      INTRINSIC          MAX
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters.
*
      INFO = 0
      IF     ( M.LT.0 )THEN
         INFO = 1
      ELSE IF( N.LT.0 )THEN
         INFO = 2
      ELSE IF( INCX.EQ.0 )THEN
         INFO = 5
      ELSE IF( INCY.EQ.0 )THEN
         INFO = 7
      ELSE IF( LDA.LT.MAX( 1, M ) )THEN
         INFO = 9
      END IF
      IF( INFO.NE.0 )THEN
         CALL XERBLA( 'DGER  ', INFO )
         RETURN
      END IF
*
*     Quick return if possible.
*
      IF( ( M.EQ.0 ).OR.( N.EQ.0 ).OR.( ALPHA.EQ.ZERO ) )
     $   RETURN
*
*     Start the operations. In this version the elements of A are
*     accessed sequentially with one pass through A.
*
      IF( INCY.GT.0 )THEN
         JY = 1
      ELSE
         JY = 1 - ( N - 1 )*INCY
      END IF
      IF( INCX.EQ.1 )THEN
         DO 20, J = 1, N
            IF( Y( JY ).NE.ZERO )THEN
               TEMP = ALPHA*Y( JY )
               DO 10, I = 1, M
                  A( I, J ) = A( I, J ) + X( I )*TEMP
   10          CONTINUE
            END IF
            JY = JY + INCY
   20    CONTINUE
      ELSE
         IF( INCX.GT.0 )THEN
            KX = 1
         ELSE
            KX = 1 - ( M - 1 )*INCX
         END IF
         DO 40, J = 1, N
            IF( Y( JY ).NE.ZERO )THEN
               TEMP = ALPHA*Y( JY )
               IX   = KX
               DO 30, I = 1, M
                  A( I, J ) = A( I, J ) + X( IX )*TEMP
                  IX        = IX        + INCX
   30          CONTINUE
            END IF
            JY = JY + INCY
   40    CONTINUE
      END IF
*
      RETURN
*
*     End of DGER  .
*
      END
      subroutine  dscal(n,da,dx,incx)
c
c     scales a vector by a constant.
c     uses unrolled loops for increment equal to one.
c     jack dongarra, linpack, 3/11/78.
c     modified 3/93 to return if incx .le. 0.
c     modified 12/3/93, array(1) declarations changed to array(*)
c
      double precision da,dx(*)
      integer i,incx,m,mp1,n,nincx
c
      if( n.le.0 .or. incx.le.0 )return
      if(incx.eq.1)go to 20
c
c        code for increment not equal to 1
c
      nincx = n*incx
      do 10 i = 1,nincx,incx
        dx(i) = da*dx(i)
   10 continue
      return
c
c        code for increment equal to 1
c
c
c        clean-up loop
c
   20 m = mod(n,5)
      if( m .eq. 0 ) go to 40
      do 30 i = 1,m
        dx(i) = da*dx(i)
   30 continue
      if( n .lt. 5 ) return
   40 mp1 = m + 1
      do 50 i = mp1,n,5
        dx(i) = da*dx(i)
        dx(i + 1) = da*dx(i + 1)
        dx(i + 2) = da*dx(i + 2)
        dx(i + 3) = da*dx(i + 3)
        dx(i + 4) = da*dx(i + 4)
   50 continue
      return
      end
      SUBROUTINE DTRMM ( SIDE, UPLO, TRANSA, DIAG, M, N, ALPHA, A, LDA,
     $                   B, LDB )
*     .. Scalar Arguments ..
      CHARACTER*1        SIDE, UPLO, TRANSA, DIAG
      INTEGER            M, N, LDA, LDB
      DOUBLE PRECISION   ALPHA
*     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), B( LDB, * )
*     ..
*
*  Purpose
*  =======
*
*  DTRMM  performs one of the matrix-matrix operations
*
*     B := alpha*op( A )*B,   or   B := alpha*B*op( A ),
*
*  where  alpha  is a scalar,  B  is an m by n matrix,  A  is a unit, or
*  non-unit,  upper or lower triangular matrix  and  op( A )  is one  of
*
*     op( A ) = A   or   op( A ) = A'.
*
*  Parameters
*  ==========
*
*  SIDE   - CHARACTER*1.
*           On entry,  SIDE specifies whether  op( A ) multiplies B from
*           the left or right as follows:
*
*              SIDE = 'L' or 'l'   B := alpha*op( A )*B.
*
*              SIDE = 'R' or 'r'   B := alpha*B*op( A ).
*
*           Unchanged on exit.
*
*  UPLO   - CHARACTER*1.
*           On entry, UPLO specifies whether the matrix A is an upper or
*           lower triangular matrix as follows:
*
*              UPLO = 'U' or 'u'   A is an upper triangular matrix.
*
*              UPLO = 'L' or 'l'   A is a lower triangular matrix.
*
*           Unchanged on exit.
*
*  TRANSA - CHARACTER*1.
*           On entry, TRANSA specifies the form of op( A ) to be used in
*           the matrix multiplication as follows:
*
*              TRANSA = 'N' or 'n'   op( A ) = A.
*
*              TRANSA = 'T' or 't'   op( A ) = A'.
*
*              TRANSA = 'C' or 'c'   op( A ) = A'.
*
*           Unchanged on exit.
*
*  DIAG   - CHARACTER*1.
*           On entry, DIAG specifies whether or not A is unit triangular
*           as follows:
*
*              DIAG = 'U' or 'u'   A is assumed to be unit triangular.
*
*              DIAG = 'N' or 'n'   A is not assumed to be unit
*                                  triangular.
*
*           Unchanged on exit.
*
*  M      - INTEGER.
*           On entry, M specifies the number of rows of B. M must be at
*           least zero.
*           Unchanged on exit.
*
*  N      - INTEGER.
*           On entry, N specifies the number of columns of B.  N must be
*           at least zero.
*           Unchanged on exit.
*
*  ALPHA  - DOUBLE PRECISION.
*           On entry,  ALPHA specifies the scalar  alpha. When  alpha is
*           zero then  A is not referenced and  B need not be set before
*           entry.
*           Unchanged on exit.
*
*  A      - DOUBLE PRECISION array of DIMENSION ( LDA, k ), where k is m
*           when  SIDE = 'L' or 'l'  and is  n  when  SIDE = 'R' or 'r'.
*           Before entry  with  UPLO = 'U' or 'u',  the  leading  k by k
*           upper triangular part of the array  A must contain the upper
*           triangular matrix  and the strictly lower triangular part of
*           A is not referenced.
*           Before entry  with  UPLO = 'L' or 'l',  the  leading  k by k
*           lower triangular part of the array  A must contain the lower
*           triangular matrix  and the strictly upper triangular part of
*           A is not referenced.
*           Note that when  DIAG = 'U' or 'u',  the diagonal elements of
*           A  are not referenced either,  but are assumed to be  unity.
*           Unchanged on exit.
*
*  LDA    - INTEGER.
*           On entry, LDA specifies the first dimension of A as declared
*           in the calling (sub) program.  When  SIDE = 'L' or 'l'  then
*           LDA  must be at least  max( 1, m ),  when  SIDE = 'R' or 'r'
*           then LDA must be at least max( 1, n ).
*           Unchanged on exit.
*
*  B      - DOUBLE PRECISION array of DIMENSION ( LDB, n ).
*           Before entry,  the leading  m by n part of the array  B must
*           contain the matrix  B,  and  on exit  is overwritten  by the
*           transformed matrix.
*
*  LDB    - INTEGER.
*           On entry, LDB specifies the first dimension of B as declared
*           in  the  calling  (sub)  program.   LDB  must  be  at  least
*           max( 1, m ).
*           Unchanged on exit.
*
*
*  Level 3 Blas routine.
*
*  -- Written on 8-February-1989.
*     Jack Dongarra, Argonne National Laboratory.
*     Iain Duff, AERE Harwell.
*     Jeremy Du Croz, Numerical Algorithms Group Ltd.
*     Sven Hammarling, Numerical Algorithms Group Ltd.
*
*
*     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
*     .. External Subroutines ..
      EXTERNAL           XERBLA
*     .. Intrinsic Functions ..
      INTRINSIC          MAX
*     .. Local Scalars ..
      LOGICAL            LSIDE, NOUNIT, UPPER
      INTEGER            I, INFO, J, K, NROWA
      DOUBLE PRECISION   TEMP
*     .. Parameters ..
      DOUBLE PRECISION   ONE         , ZERO
      PARAMETER        ( ONE = 1.0D+0, ZERO = 0.0D+0 )
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters.
*
      LSIDE  = LSAME( SIDE  , 'L' )
      IF( LSIDE )THEN
         NROWA = M
      ELSE
         NROWA = N
      END IF
      NOUNIT = LSAME( DIAG  , 'N' )
      UPPER  = LSAME( UPLO  , 'U' )
*
      INFO   = 0
      IF(      ( .NOT.LSIDE                ).AND.
     $         ( .NOT.LSAME( SIDE  , 'R' ) )      )THEN
         INFO = 1
      ELSE IF( ( .NOT.UPPER                ).AND.
     $         ( .NOT.LSAME( UPLO  , 'L' ) )      )THEN
         INFO = 2
      ELSE IF( ( .NOT.LSAME( TRANSA, 'N' ) ).AND.
     $         ( .NOT.LSAME( TRANSA, 'T' ) ).AND.
     $         ( .NOT.LSAME( TRANSA, 'C' ) )      )THEN
         INFO = 3
      ELSE IF( ( .NOT.LSAME( DIAG  , 'U' ) ).AND.
     $         ( .NOT.LSAME( DIAG  , 'N' ) )      )THEN
         INFO = 4
      ELSE IF( M  .LT.0               )THEN
         INFO = 5
      ELSE IF( N  .LT.0               )THEN
         INFO = 6
      ELSE IF( LDA.LT.MAX( 1, NROWA ) )THEN
         INFO = 9
      ELSE IF( LDB.LT.MAX( 1, M     ) )THEN
         INFO = 11
      END IF
      IF( INFO.NE.0 )THEN
         CALL XERBLA( 'DTRMM ', INFO )
         RETURN
      END IF
*
*     Quick return if possible.
*
      IF( N.EQ.0 )
     $   RETURN
*
*     And when  alpha.eq.zero.
*
      IF( ALPHA.EQ.ZERO )THEN
         DO 20, J = 1, N
            DO 10, I = 1, M
               B( I, J ) = ZERO
   10       CONTINUE
   20    CONTINUE
         RETURN
      END IF
*
*     Start the operations.
*
      IF( LSIDE )THEN
         IF( LSAME( TRANSA, 'N' ) )THEN
*
*           Form  B := alpha*A*B.
*
            IF( UPPER )THEN
               DO 50, J = 1, N
                  DO 40, K = 1, M
                     IF( B( K, J ).NE.ZERO )THEN
                        TEMP = ALPHA*B( K, J )
                        DO 30, I = 1, K - 1
                           B( I, J ) = B( I, J ) + TEMP*A( I, K )
   30                   CONTINUE
                        IF( NOUNIT )
     $                     TEMP = TEMP*A( K, K )
                        B( K, J ) = TEMP
                     END IF
   40             CONTINUE
   50          CONTINUE
            ELSE
               DO 80, J = 1, N
                  DO 70 K = M, 1, -1
                     IF( B( K, J ).NE.ZERO )THEN
                        TEMP      = ALPHA*B( K, J )
                        B( K, J ) = TEMP
                        IF( NOUNIT )
     $                     B( K, J ) = B( K, J )*A( K, K )
                        DO 60, I = K + 1, M
                           B( I, J ) = B( I, J ) + TEMP*A( I, K )
   60                   CONTINUE
                     END IF
   70             CONTINUE
   80          CONTINUE
            END IF
         ELSE
*
*           Form  B := alpha*B*A'.
*
            IF( UPPER )THEN
               DO 110, J = 1, N
                  DO 100, I = M, 1, -1
                     TEMP = B( I, J )
                     IF( NOUNIT )
     $                  TEMP = TEMP*A( I, I )
                     DO 90, K = 1, I - 1
                        TEMP = TEMP + A( K, I )*B( K, J )
   90                CONTINUE
                     B( I, J ) = ALPHA*TEMP
  100             CONTINUE
  110          CONTINUE
            ELSE
               DO 140, J = 1, N
                  DO 130, I = 1, M
                     TEMP = B( I, J )
                     IF( NOUNIT )
     $                  TEMP = TEMP*A( I, I )
                     DO 120, K = I + 1, M
                        TEMP = TEMP + A( K, I )*B( K, J )
  120                CONTINUE
                     B( I, J ) = ALPHA*TEMP
  130             CONTINUE
  140          CONTINUE
            END IF
         END IF
      ELSE
         IF( LSAME( TRANSA, 'N' ) )THEN
*
*           Form  B := alpha*B*A.
*
            IF( UPPER )THEN
               DO 180, J = N, 1, -1
                  TEMP = ALPHA
                  IF( NOUNIT )
     $               TEMP = TEMP*A( J, J )
                  DO 150, I = 1, M
                     B( I, J ) = TEMP*B( I, J )
  150             CONTINUE
                  DO 170, K = 1, J - 1
                     IF( A( K, J ).NE.ZERO )THEN
                        TEMP = ALPHA*A( K, J )
                        DO 160, I = 1, M
                           B( I, J ) = B( I, J ) + TEMP*B( I, K )
  160                   CONTINUE
                     END IF
  170             CONTINUE
  180          CONTINUE
            ELSE
               DO 220, J = 1, N
                  TEMP = ALPHA
                  IF( NOUNIT )
     $               TEMP = TEMP*A( J, J )
                  DO 190, I = 1, M
                     B( I, J ) = TEMP*B( I, J )
  190             CONTINUE
                  DO 210, K = J + 1, N
                     IF( A( K, J ).NE.ZERO )THEN
                        TEMP = ALPHA*A( K, J )
                        DO 200, I = 1, M
                           B( I, J ) = B( I, J ) + TEMP*B( I, K )
  200                   CONTINUE
                     END IF
  210             CONTINUE
  220          CONTINUE
            END IF
         ELSE
*
*           Form  B := alpha*B*A'.
*
            IF( UPPER )THEN
               DO 260, K = 1, N
                  DO 240, J = 1, K - 1
                     IF( A( J, K ).NE.ZERO )THEN
                        TEMP = ALPHA*A( J, K )
                        DO 230, I = 1, M
                           B( I, J ) = B( I, J ) + TEMP*B( I, K )
  230                   CONTINUE
                     END IF
  240             CONTINUE
                  TEMP = ALPHA
                  IF( NOUNIT )
     $               TEMP = TEMP*A( K, K )
                  IF( TEMP.NE.ONE )THEN
                     DO 250, I = 1, M
                        B( I, K ) = TEMP*B( I, K )
  250                CONTINUE
                  END IF
  260          CONTINUE
            ELSE
               DO 300, K = N, 1, -1
                  DO 280, J = K + 1, N
                     IF( A( J, K ).NE.ZERO )THEN
                        TEMP = ALPHA*A( J, K )
                        DO 270, I = 1, M
                           B( I, J ) = B( I, J ) + TEMP*B( I, K )
  270                   CONTINUE
                     END IF
  280             CONTINUE
                  TEMP = ALPHA
                  IF( NOUNIT )
     $               TEMP = TEMP*A( K, K )
                  IF( TEMP.NE.ONE )THEN
                     DO 290, I = 1, M
                        B( I, K ) = TEMP*B( I, K )
  290                CONTINUE
                  END IF
  300          CONTINUE
            END IF
         END IF
      END IF
*
      RETURN
*
*     End of DTRMM .
*
      END
      SUBROUTINE DTRMV ( UPLO, TRANS, DIAG, N, A, LDA, X, INCX )
*     .. Scalar Arguments ..
      INTEGER            INCX, LDA, N
      CHARACTER*1        DIAG, TRANS, UPLO
*     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), X( * )
*     ..
*
*  Purpose
*  =======
*
*  DTRMV  performs one of the matrix-vector operations
*
*     x := A*x,   or   x := A'*x,
*
*  where x is an n element vector and  A is an n by n unit, or non-unit,
*  upper or lower triangular matrix.
*
*  Parameters
*  ==========
*
*  UPLO   - CHARACTER*1.
*           On entry, UPLO specifies whether the matrix is an upper or
*           lower triangular matrix as follows:
*
*              UPLO = 'U' or 'u'   A is an upper triangular matrix.
*
*              UPLO = 'L' or 'l'   A is a lower triangular matrix.
*
*           Unchanged on exit.
*
*  TRANS  - CHARACTER*1.
*           On entry, TRANS specifies the operation to be performed as
*           follows:
*
*              TRANS = 'N' or 'n'   x := A*x.
*
*              TRANS = 'T' or 't'   x := A'*x.
*
*              TRANS = 'C' or 'c'   x := A'*x.
*
*           Unchanged on exit.
*
*  DIAG   - CHARACTER*1.
*           On entry, DIAG specifies whether or not A is unit
*           triangular as follows:
*
*              DIAG = 'U' or 'u'   A is assumed to be unit triangular.
*
*              DIAG = 'N' or 'n'   A is not assumed to be unit
*                                  triangular.
*
*           Unchanged on exit.
*
*  N      - INTEGER.
*           On entry, N specifies the order of the matrix A.
*           N must be at least zero.
*           Unchanged on exit.
*
*  A      - DOUBLE PRECISION array of DIMENSION ( LDA, n ).
*           Before entry with  UPLO = 'U' or 'u', the leading n by n
*           upper triangular part of the array A must contain the upper
*           triangular matrix and the strictly lower triangular part of
*           A is not referenced.
*           Before entry with UPLO = 'L' or 'l', the leading n by n
*           lower triangular part of the array A must contain the lower
*           triangular matrix and the strictly upper triangular part of
*           A is not referenced.
*           Note that when  DIAG = 'U' or 'u', the diagonal elements of
*           A are not referenced either, but are assumed to be unity.
*           Unchanged on exit.
*
*  LDA    - INTEGER.
*           On entry, LDA specifies the first dimension of A as declared
*           in the calling (sub) program. LDA must be at least
*           max( 1, n ).
*           Unchanged on exit.
*
*  X      - DOUBLE PRECISION array of dimension at least
*           ( 1 + ( n - 1 )*abs( INCX ) ).
*           Before entry, the incremented array X must contain the n
*           element vector x. On exit, X is overwritten with the
*           tranformed vector x.
*
*  INCX   - INTEGER.
*           On entry, INCX specifies the increment for the elements of
*           X. INCX must not be zero.
*           Unchanged on exit.
*
*
*  Level 2 Blas routine.
*
*  -- Written on 22-October-1986.
*     Jack Dongarra, Argonne National Lab.
*     Jeremy Du Croz, Nag Central Office.
*     Sven Hammarling, Nag Central Office.
*     Richard Hanson, Sandia National Labs.
*
*
*     .. Parameters ..
      DOUBLE PRECISION   ZERO
      PARAMETER        ( ZERO = 0.0D+0 )
*     .. Local Scalars ..
      DOUBLE PRECISION   TEMP
      INTEGER            I, INFO, IX, J, JX, KX
      LOGICAL            NOUNIT
*     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
*     .. External Subroutines ..
      EXTERNAL           XERBLA
*     .. Intrinsic Functions ..
      INTRINSIC          MAX
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters.
*
      INFO = 0
      IF     ( .NOT.LSAME( UPLO , 'U' ).AND.
     $         .NOT.LSAME( UPLO , 'L' )      )THEN
         INFO = 1
      ELSE IF( .NOT.LSAME( TRANS, 'N' ).AND.
     $         .NOT.LSAME( TRANS, 'T' ).AND.
     $         .NOT.LSAME( TRANS, 'C' )      )THEN
         INFO = 2
      ELSE IF( .NOT.LSAME( DIAG , 'U' ).AND.
     $         .NOT.LSAME( DIAG , 'N' )      )THEN
         INFO = 3
      ELSE IF( N.LT.0 )THEN
         INFO = 4
      ELSE IF( LDA.LT.MAX( 1, N ) )THEN
         INFO = 6
      ELSE IF( INCX.EQ.0 )THEN
         INFO = 8
      END IF
      IF( INFO.NE.0 )THEN
         CALL XERBLA( 'DTRMV ', INFO )
         RETURN
      END IF
*
*     Quick return if possible.
*
      IF( N.EQ.0 )
     $   RETURN
*
      NOUNIT = LSAME( DIAG, 'N' )
*
*     Set up the start point in X if the increment is not unity. This
*     will be  ( N - 1 )*INCX  too small for descending loops.
*
      IF( INCX.LE.0 )THEN
         KX = 1 - ( N - 1 )*INCX
      ELSE IF( INCX.NE.1 )THEN
         KX = 1
      END IF
*
*     Start the operations. In this version the elements of A are
*     accessed sequentially with one pass through A.
*
      IF( LSAME( TRANS, 'N' ) )THEN
*
*        Form  x := A*x.
*
         IF( LSAME( UPLO, 'U' ) )THEN
            IF( INCX.EQ.1 )THEN
               DO 20, J = 1, N
                  IF( X( J ).NE.ZERO )THEN
                     TEMP = X( J )
                     DO 10, I = 1, J - 1
                        X( I ) = X( I ) + TEMP*A( I, J )
   10                CONTINUE
                     IF( NOUNIT )
     $                  X( J ) = X( J )*A( J, J )
                  END IF
   20          CONTINUE
            ELSE
               JX = KX
               DO 40, J = 1, N
                  IF( X( JX ).NE.ZERO )THEN
                     TEMP = X( JX )
                     IX   = KX
                     DO 30, I = 1, J - 1
                        X( IX ) = X( IX ) + TEMP*A( I, J )
                        IX      = IX      + INCX
   30                CONTINUE
                     IF( NOUNIT )
     $                  X( JX ) = X( JX )*A( J, J )
                  END IF
                  JX = JX + INCX
   40          CONTINUE
            END IF
         ELSE
            IF( INCX.EQ.1 )THEN
               DO 60, J = N, 1, -1
                  IF( X( J ).NE.ZERO )THEN
                     TEMP = X( J )
                     DO 50, I = N, J + 1, -1
                        X( I ) = X( I ) + TEMP*A( I, J )
   50                CONTINUE
                     IF( NOUNIT )
     $                  X( J ) = X( J )*A( J, J )
                  END IF
   60          CONTINUE
            ELSE
               KX = KX + ( N - 1 )*INCX
               JX = KX
               DO 80, J = N, 1, -1
                  IF( X( JX ).NE.ZERO )THEN
                     TEMP = X( JX )
                     IX   = KX
                     DO 70, I = N, J + 1, -1
                        X( IX ) = X( IX ) + TEMP*A( I, J )
                        IX      = IX      - INCX
   70                CONTINUE
                     IF( NOUNIT )
     $                  X( JX ) = X( JX )*A( J, J )
                  END IF
                  JX = JX - INCX
   80          CONTINUE
            END IF
         END IF
      ELSE
*
*        Form  x := A'*x.
*
         IF( LSAME( UPLO, 'U' ) )THEN
            IF( INCX.EQ.1 )THEN
               DO 100, J = N, 1, -1
                  TEMP = X( J )
                  IF( NOUNIT )
     $               TEMP = TEMP*A( J, J )
                  DO 90, I = J - 1, 1, -1
                     TEMP = TEMP + A( I, J )*X( I )
   90             CONTINUE
                  X( J ) = TEMP
  100          CONTINUE
            ELSE
               JX = KX + ( N - 1 )*INCX
               DO 120, J = N, 1, -1
                  TEMP = X( JX )
                  IX   = JX
                  IF( NOUNIT )
     $               TEMP = TEMP*A( J, J )
                  DO 110, I = J - 1, 1, -1
                     IX   = IX   - INCX
                     TEMP = TEMP + A( I, J )*X( IX )
  110             CONTINUE
                  X( JX ) = TEMP
                  JX      = JX   - INCX
  120          CONTINUE
            END IF
         ELSE
            IF( INCX.EQ.1 )THEN
               DO 140, J = 1, N
                  TEMP = X( J )
                  IF( NOUNIT )
     $               TEMP = TEMP*A( J, J )
                  DO 130, I = J + 1, N
                     TEMP = TEMP + A( I, J )*X( I )
  130             CONTINUE
                  X( J ) = TEMP
  140          CONTINUE
            ELSE
               JX = KX
               DO 160, J = 1, N
                  TEMP = X( JX )
                  IX   = JX
                  IF( NOUNIT )
     $               TEMP = TEMP*A( J, J )
                  DO 150, I = J + 1, N
                     IX   = IX   + INCX
                     TEMP = TEMP + A( I, J )*X( IX )
  150             CONTINUE
                  X( JX ) = TEMP
                  JX      = JX   + INCX
  160          CONTINUE
            END IF
         END IF
      END IF
*
      RETURN
*
*     End of DTRMV .
*
      END
      SUBROUTINE DGETRF( M, N, A, LDA, IPIV, INFO )
*
*  -- LAPACK routine (version 2.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     March 31, 1993
*
*     .. Scalar Arguments ..
      INTEGER            INFO, LDA, M, N
*     ..
*     .. Array Arguments ..
      INTEGER            IPIV( * )
      DOUBLE PRECISION   A( LDA, * )
*     ..
*
*  Purpose
*  =======
*
*  DGETRF computes an LU factorization of a general M-by-N matrix A
*  using partial pivoting with row interchanges.
*
*  The factorization has the form
*     A = P * L * U
*  where P is a permutation matrix, L is lower triangular with unit
*  diagonal elements (lower trapezoidal if m > n), and U is upper
*  triangular (upper trapezoidal if m < n).
*
*  This is the right-looking Level 3 BLAS version of the algorithm.
*
*  Arguments
*  =========
*
*  M       (input) INTEGER
*          The number of rows of the matrix A.  M >= 0.
*
*  N       (input) INTEGER
*          The number of columns of the matrix A.  N >= 0.
*
*  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
*          On entry, the M-by-N matrix to be factored.
*          On exit, the factors L and U from the factorization
*          A = P*L*U; the unit diagonal elements of L are not stored.
*
*  LDA     (input) INTEGER
*          The leading dimension of the array A.  LDA >= max(1,M).
*
*  IPIV    (output) INTEGER array, dimension (min(M,N))
*          The pivot indices; for 1 <= i <= min(M,N), row i of the
*          matrix was interchanged with row IPIV(i).
*
*  INFO    (output) INTEGER
*          = 0:  successful exit
*          < 0:  if INFO = -i, the i-th argument had an illegal value
*          > 0:  if INFO = i, U(i,i) is exactly zero. The factorization
*                has been completed, but the factor U is exactly
*                singular, and division by zero will occur if it is used
*                to solve a system of equations.
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ONE
      PARAMETER          ( ONE = 1.0D+0 )
*     ..
*     .. Local Scalars ..
      INTEGER            I, IINFO, J, JB, NB
*     ..
*     .. External Subroutines ..
      EXTERNAL           DGEMM, DGETF2, DLASWP, DTRSM, XERBLA
*     ..
*     .. External Functions ..
      INTEGER            ILAENV
      EXTERNAL           ILAENV
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters.
*
      INFO = 0
      IF( M.LT.0 ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( LDA.LT.MAX( 1, M ) ) THEN
         INFO = -4
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DGETRF', -INFO )
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( M.EQ.0 .OR. N.EQ.0 )
     $   RETURN
*
*     Determine the block size for this environment.
*
      NB = ILAENV( 1, 'DGETRF', ' ', M, N, -1, -1 )
      IF( NB.LE.1 .OR. NB.GE.MIN( M, N ) ) THEN
*
*        Use unblocked code.
*
         CALL DGETF2( M, N, A, LDA, IPIV, INFO )
      ELSE
*
*        Use blocked code.
*
         DO 20 J = 1, MIN( M, N ), NB
            JB = MIN( MIN( M, N )-J+1, NB )
*
*           Factor diagonal and subdiagonal blocks and test for exact
*           singularity.
*
            CALL DGETF2( M-J+1, JB, A( J, J ), LDA, IPIV( J ), IINFO )
*
*           Adjust INFO and the pivot indices.
*
            IF( INFO.EQ.0 .AND. IINFO.GT.0 )
     $         INFO = IINFO + J - 1
            DO 10 I = J, MIN( M, J+JB-1 )
               IPIV( I ) = J - 1 + IPIV( I )
   10       CONTINUE
*
*           Apply interchanges to columns 1:J-1.
*
            CALL DLASWP( J-1, A, LDA, J, J+JB-1, IPIV, 1 )
*
            IF( J+JB.LE.N ) THEN
*
*              Apply interchanges to columns J+JB:N.
*
               CALL DLASWP( N-J-JB+1, A( 1, J+JB ), LDA, J, J+JB-1,
     $                      IPIV, 1 )
*
*              Compute block row of U.
*
               CALL DTRSM( 'Left', 'Lower', 'No transpose', 'Unit', JB,
     $                     N-J-JB+1, ONE, A( J, J ), LDA, A( J, J+JB ),
     $                     LDA )
               IF( J+JB.LE.M ) THEN
*
*                 Update trailing submatrix.
*
                  CALL DGEMM( 'No transpose', 'No transpose', M-J-JB+1,
     $                        N-J-JB+1, JB, -ONE, A( J+JB, J ), LDA,
     $                        A( J, J+JB ), LDA, ONE, A( J+JB, J+JB ),
     $                        LDA )
               END IF
            END IF
   20    CONTINUE
      END IF
      RETURN
*
*     End of DGETRF
*
      END
      SUBROUTINE DGETRI( N, A, LDA, IPIV, WORK, LWORK, INFO )
*
*  -- LAPACK routine (version 2.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     September 30, 1994
*
*     .. Scalar Arguments ..
      INTEGER            INFO, LDA, LWORK, N
*     ..
*     .. Array Arguments ..
      INTEGER            IPIV( * )
      DOUBLE PRECISION   A( LDA, * ), WORK( LWORK )
*     ..
*
*  Purpose
*  =======
*
*  DGETRI computes the inverse of a matrix using the LU factorization
*  computed by DGETRF.
*
*  This method inverts U and then computes inv(A) by solving the system
*  inv(A)*L = inv(U) for inv(A).
*
*  Arguments
*  =========
*
*  N       (input) INTEGER
*          The order of the matrix A.  N >= 0.
*
*  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
*          On entry, the factors L and U from the factorization
*          A = P*L*U as computed by DGETRF.
*          On exit, if INFO = 0, the inverse of the original matrix A.
*
*  LDA     (input) INTEGER
*          The leading dimension of the array A.  LDA >= max(1,N).
*
*  IPIV    (input) INTEGER array, dimension (N)
*          The pivot indices from DGETRF; for 1<=i<=N, row i of the
*          matrix was interchanged with row IPIV(i).
*
*  WORK    (workspace/output) DOUBLE PRECISION array, dimension (LWORK)
*          On exit, if INFO=0, then WORK(1) returns the optimal LWORK.
*
*  LWORK   (input) INTEGER
*          The dimension of the array WORK.  LWORK >= max(1,N).
*          For optimal performance LWORK >= N*NB, where NB is
*          the optimal blocksize returned by ILAENV.
*
*  INFO    (output) INTEGER
*          = 0:  successful exit
*          < 0:  if INFO = -i, the i-th argument had an illegal value
*          > 0:  if INFO = i, U(i,i) is exactly zero; the matrix is
*                singular and its inverse could not be computed.
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE
      PARAMETER          ( ZERO = 0.0D+0, ONE = 1.0D+0 )
*     ..
*     .. Local Scalars ..
      INTEGER            I, IWS, J, JB, JJ, JP, LDWORK, NB, NBMIN, NN
*     ..
*     .. External Functions ..
      INTEGER            ILAENV
      EXTERNAL           ILAENV
*     ..
*     .. External Subroutines ..
      EXTERNAL           DGEMM, DGEMV, DSWAP, DTRSM, DTRTRI, XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters.
*
      INFO = 0
      WORK( 1 ) = MAX( N, 1 )
      IF( N.LT.0 ) THEN
         INFO = -1
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -3
      ELSE IF( LWORK.LT.MAX( 1, N ) ) THEN
         INFO = -6
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DGETRI', -INFO )
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( N.EQ.0 )
     $   RETURN
*
*     Form inv(U).  If INFO > 0 from DTRTRI, then U is singular,
*     and the inverse is not computed.
*
      CALL DTRTRI( 'Upper', 'Non-unit', N, A, LDA, INFO )
      IF( INFO.GT.0 )
     $   RETURN
*
*     Determine the block size for this environment.
*
      NB = ILAENV( 1, 'DGETRI', ' ', N, -1, -1, -1 )
      NBMIN = 2
      LDWORK = N
      IF( NB.GT.1 .AND. NB.LT.N ) THEN
         IWS = MAX( LDWORK*NB, 1 )
         IF( LWORK.LT.IWS ) THEN
            NB = LWORK / LDWORK
            NBMIN = MAX( 2, ILAENV( 2, 'DGETRI', ' ', N, -1, -1, -1 ) )
         END IF
      ELSE
         IWS = N
      END IF
*
*     Solve the equation inv(A)*L = inv(U) for inv(A).
*
      IF( NB.LT.NBMIN .OR. NB.GE.N ) THEN
*
*        Use unblocked code.
*
         DO 20 J = N, 1, -1
*
*           Copy current column of L to WORK and replace with zeros.
*
            DO 10 I = J + 1, N
               WORK( I ) = A( I, J )
               A( I, J ) = ZERO
   10       CONTINUE
*
*           Compute current column of inv(A).
*
            IF( J.LT.N )
     $         CALL DGEMV( 'No transpose', N, N-J, -ONE, A( 1, J+1 ),
     $                     LDA, WORK( J+1 ), 1, ONE, A( 1, J ), 1 )
   20    CONTINUE
      ELSE
*
*        Use blocked code.
*
         NN = ( ( N-1 ) / NB )*NB + 1
         DO 50 J = NN, 1, -NB
            JB = MIN( NB, N-J+1 )
*
*           Copy current block column of L to WORK and replace with
*           zeros.
*
            DO 40 JJ = J, J + JB - 1
               DO 30 I = JJ + 1, N
                  WORK( I+( JJ-J )*LDWORK ) = A( I, JJ )
                  A( I, JJ ) = ZERO
   30          CONTINUE
   40       CONTINUE
*
*           Compute current block column of inv(A).
*
            IF( J+JB.LE.N )
     $         CALL DGEMM( 'No transpose', 'No transpose', N, JB,
     $                     N-J-JB+1, -ONE, A( 1, J+JB ), LDA,
     $                     WORK( J+JB ), LDWORK, ONE, A( 1, J ), LDA )
            CALL DTRSM( 'Right', 'Lower', 'No transpose', 'Unit', N, JB,
     $                  ONE, WORK( J ), LDWORK, A( 1, J ), LDA )
   50    CONTINUE
      END IF
*
*     Apply column interchanges.
*
      DO 60 J = N - 1, 1, -1
         JP = IPIV( J )
         IF( JP.NE.J )
     $      CALL DSWAP( N, A( 1, J ), 1, A( 1, JP ), 1 )
   60 CONTINUE
*
      WORK( 1 ) = IWS
      RETURN
*
*     End of DGETRI
*
      END
      SUBROUTINE DGETF2( M, N, A, LDA, IPIV, INFO )
*
*  -- LAPACK routine (version 2.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     June 30, 1992
*
*     .. Scalar Arguments ..
      INTEGER            INFO, LDA, M, N
*     ..
*     .. Array Arguments ..
      INTEGER            IPIV( * )
      DOUBLE PRECISION   A( LDA, * )
*     ..
*
*  Purpose
*  =======
*
*  DGETF2 computes an LU factorization of a general m-by-n matrix A
*  using partial pivoting with row interchanges.
*
*  The factorization has the form
*     A = P * L * U
*  where P is a permutation matrix, L is lower triangular with unit
*  diagonal elements (lower trapezoidal if m > n), and U is upper
*  triangular (upper trapezoidal if m < n).
*
*  This is the right-looking Level 2 BLAS version of the algorithm.
*
*  Arguments
*  =========
*
*  M       (input) INTEGER
*          The number of rows of the matrix A.  M >= 0.
*
*  N       (input) INTEGER
*          The number of columns of the matrix A.  N >= 0.
*
*  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
*          On entry, the m by n matrix to be factored.
*          On exit, the factors L and U from the factorization
*          A = P*L*U; the unit diagonal elements of L are not stored.
*
*  LDA     (input) INTEGER
*          The leading dimension of the array A.  LDA >= max(1,M).
*
*  IPIV    (output) INTEGER array, dimension (min(M,N))
*          The pivot indices; for 1 <= i <= min(M,N), row i of the
*          matrix was interchanged with row IPIV(i).
*
*  INFO    (output) INTEGER
*          = 0: successful exit
*          < 0: if INFO = -k, the k-th argument had an illegal value
*          > 0: if INFO = k, U(k,k) is exactly zero. The factorization
*               has been completed, but the factor U is exactly
*               singular, and division by zero will occur if it is used
*               to solve a system of equations.
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ONE, ZERO
      PARAMETER          ( ONE = 1.0D+0, ZERO = 0.0D+0 )
*     ..
*     .. Local Scalars ..
      INTEGER            J, JP
*     ..
*     .. External Functions ..
      INTEGER            IDAMAX
      EXTERNAL           IDAMAX
*     ..
*     .. External Subroutines ..
      EXTERNAL           DGER, DSCAL, DSWAP, XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters.
*
      INFO = 0
      IF( M.LT.0 ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( LDA.LT.MAX( 1, M ) ) THEN
         INFO = -4
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DGETF2', -INFO )
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( M.EQ.0 .OR. N.EQ.0 )
     $   RETURN
*
      DO 10 J = 1, MIN( M, N )
*
*        Find pivot and test for singularity.
*
         JP = J - 1 + IDAMAX( M-J+1, A( J, J ), 1 )
         IPIV( J ) = JP
         IF( A( JP, J ).NE.ZERO ) THEN
*
*           Apply the interchange to columns 1:N.
*
            IF( JP.NE.J )
     $         CALL DSWAP( N, A( J, 1 ), LDA, A( JP, 1 ), LDA )
*
*           Compute elements J+1:M of J-th column.
*
            IF( J.LT.M )
     $         CALL DSCAL( M-J, ONE / A( J, J ), A( J+1, J ), 1 )
*
         ELSE IF( INFO.EQ.0 ) THEN
*
            INFO = J
         END IF
*
         IF( J.LT.MIN( M, N ) ) THEN
*
*           Update trailing submatrix.
*
            CALL DGER( M-J, N-J, -ONE, A( J+1, J ), 1, A( J, J+1 ), LDA,
     $                 A( J+1, J+1 ), LDA )
         END IF
   10 CONTINUE
      RETURN
*
*     End of DGETF2
*
      END
      SUBROUTINE DLASWP( N, A, LDA, K1, K2, IPIV, INCX )
*
*  -- LAPACK auxiliary routine (version 2.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     October 31, 1992
*
*     .. Scalar Arguments ..
      INTEGER            INCX, K1, K2, LDA, N
*     ..
*     .. Array Arguments ..
      INTEGER            IPIV( * )
      DOUBLE PRECISION   A( LDA, * )
*     ..
*
*  Purpose
*  =======
*
*  DLASWP performs a series of row interchanges on the matrix A.
*  One row interchange is initiated for each of rows K1 through K2 of A.
*
*  Arguments
*  =========
*
*  N       (input) INTEGER
*          The number of columns of the matrix A.
*
*  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
*          On entry, the matrix of column dimension N to which the row
*          interchanges will be applied.
*          On exit, the permuted matrix.
*
*  LDA     (input) INTEGER
*          The leading dimension of the array A.
*
*  K1      (input) INTEGER
*          The first element of IPIV for which a row interchange will
*          be done.
*
*  K2      (input) INTEGER
*          The last element of IPIV for which a row interchange will
*          be done.
*
*  IPIV    (input) INTEGER array, dimension (M*abs(INCX))
*          The vector of pivot indices.  Only the elements in positions
*          K1 through K2 of IPIV are accessed.
*          IPIV(K) = L implies rows K and L are to be interchanged.
*
*  INCX    (input) INTEGER
*          The increment between successive values of IPIV.  If IPIV
*          is negative, the pivots are applied in reverse order.
*
* =====================================================================
*
*     .. Local Scalars ..
      INTEGER            I, IP, IX
*     ..
*     .. External Subroutines ..
      EXTERNAL           DSWAP
*     ..
*     .. Executable Statements ..
*
*     Interchange row I with row IPIV(I) for each of rows K1 through K2.
*
      IF( INCX.EQ.0 )
     $   RETURN
      IF( INCX.GT.0 ) THEN
         IX = K1
      ELSE
         IX = 1 + ( 1-K2 )*INCX
      END IF
      IF( INCX.EQ.1 ) THEN
         DO 10 I = K1, K2
            IP = IPIV( I )
            IF( IP.NE.I )
     $         CALL DSWAP( N, A( I, 1 ), LDA, A( IP, 1 ), LDA )
   10    CONTINUE
      ELSE IF( INCX.GT.1 ) THEN
         DO 20 I = K1, K2
            IP = IPIV( IX )
            IF( IP.NE.I )
     $         CALL DSWAP( N, A( I, 1 ), LDA, A( IP, 1 ), LDA )
            IX = IX + INCX
   20    CONTINUE
      ELSE IF( INCX.LT.0 ) THEN
         DO 30 I = K2, K1, -1
            IP = IPIV( IX )
            IF( IP.NE.I )
     $         CALL DSWAP( N, A( I, 1 ), LDA, A( IP, 1 ), LDA )
            IX = IX + INCX
   30    CONTINUE
      END IF
*
      RETURN
*
*     End of DLASWP
*
      END
      SUBROUTINE XERBLA( SRNAME, INFO )
*
*  -- LAPACK auxiliary routine (preliminary version) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     February 29, 1992
*
*     .. Scalar Arguments ..
      CHARACTER*6        SRNAME
      INTEGER            INFO
*     ..
*
*  Purpose
*  =======
*
*  XERBLA  is an error handler for the LAPACK routines.
*  It is called by an LAPACK routine if an input parameter has an
*  invalid value.  A message is printed and execution stops.
*
*  Installers may consider modifying the STOP statement in order to
*  call system-specific exception-handling facilities.
*
*  Arguments
*  =========
*
*  SRNAME  (input) CHARACTER*6
*          The name of the routine which called XERBLA.
*
*  INFO    (input) INTEGER
*          The position of the invalid parameter in the parameter list
*          of the calling routine.
*
*
      WRITE( *, FMT = 9999 )SRNAME, INFO
*
      STOP
*
 9999 FORMAT( ' ** On entry to ', A6, ' parameter number ', I2, ' had ',
     $      'an illegal value' )
*
*     End of XERBLA
*
      END
      SUBROUTINE DTRTRI( UPLO, DIAG, N, A, LDA, INFO )
*
*  -- LAPACK routine (version 2.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     March 31, 1993
*
*     .. Scalar Arguments ..
      CHARACTER          DIAG, UPLO
      INTEGER            INFO, LDA, N
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * )
*     ..
*
*  Purpose
*  =======
*
*  DTRTRI computes the inverse of a real upper or lower triangular
*  matrix A.
*
*  This is the Level 3 BLAS version of the algorithm.
*
*  Arguments
*  =========
*
*  UPLO    (input) CHARACTER*1
*          = 'U':  A is upper triangular;
*          = 'L':  A is lower triangular.
*
*  DIAG    (input) CHARACTER*1
*          = 'N':  A is non-unit triangular;
*          = 'U':  A is unit triangular.
*
*  N       (input) INTEGER
*          The order of the matrix A.  N >= 0.
*
*  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
*          On entry, the triangular matrix A.  If UPLO = 'U', the
*          leading N-by-N upper triangular part of the array A contains
*          the upper triangular matrix, and the strictly lower
*          triangular part of A is not referenced.  If UPLO = 'L', the
*          leading N-by-N lower triangular part of the array A contains
*          the lower triangular matrix, and the strictly upper
*          triangular part of A is not referenced.  If DIAG = 'U', the
*          diagonal elements of A are also not referenced and are
*          assumed to be 1.
*          On exit, the (triangular) inverse of the original matrix, in
*          the same storage format.
*
*  LDA     (input) INTEGER
*          The leading dimension of the array A.  LDA >= max(1,N).
*
*  INFO    (output) INTEGER
*          = 0: successful exit
*          < 0: if INFO = -i, the i-th argument had an illegal value
*          > 0: if INFO = i, A(i,i) is exactly zero.  The triangular
*               matrix is singular and its inverse can not be computed.
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ONE, ZERO
      PARAMETER          ( ONE = 1.0D+0, ZERO = 0.0D+0 )
*     ..
*     .. Local Scalars ..
      LOGICAL            NOUNIT, UPPER
      INTEGER            J, JB, NB, NN
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      INTEGER            ILAENV
      EXTERNAL           LSAME, ILAENV
*     ..
*     .. External Subroutines ..
      EXTERNAL           DTRMM, DTRSM, DTRTI2, XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters.
*
      INFO = 0
      UPPER = LSAME( UPLO, 'U' )
      NOUNIT = LSAME( DIAG, 'N' )
      IF( .NOT.UPPER .AND. .NOT.LSAME( UPLO, 'L' ) ) THEN
         INFO = -1
      ELSE IF( .NOT.NOUNIT .AND. .NOT.LSAME( DIAG, 'U' ) ) THEN
         INFO = -2
      ELSE IF( N.LT.0 ) THEN
         INFO = -3
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -5
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DTRTRI', -INFO )
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( N.EQ.0 )
     $   RETURN
*
*     Check for singularity if non-unit.
*
      IF( NOUNIT ) THEN
         DO 10 INFO = 1, N
            IF( A( INFO, INFO ).EQ.ZERO )
     $         RETURN
   10    CONTINUE
         INFO = 0
      END IF
*
*     Determine the block size for this environment.
*
      NB = ILAENV( 1, 'DTRTRI', UPLO // DIAG, N, -1, -1, -1 )
      IF( NB.LE.1 .OR. NB.GE.N ) THEN
*
*        Use unblocked code
*
         CALL DTRTI2( UPLO, DIAG, N, A, LDA, INFO )
      ELSE
*
*        Use blocked code
*
         IF( UPPER ) THEN
*
*           Compute inverse of upper triangular matrix
*
            DO 20 J = 1, N, NB
               JB = MIN( NB, N-J+1 )
*
*              Compute rows 1:j-1 of current block column
*
               CALL DTRMM( 'Left', 'Upper', 'No transpose', DIAG, J-1,
     $                     JB, ONE, A, LDA, A( 1, J ), LDA )
               CALL DTRSM( 'Right', 'Upper', 'No transpose', DIAG, J-1,
     $                     JB, -ONE, A( J, J ), LDA, A( 1, J ), LDA )
*
*              Compute inverse of current diagonal block
*
               CALL DTRTI2( 'Upper', DIAG, JB, A( J, J ), LDA, INFO )
   20       CONTINUE
         ELSE
*
*           Compute inverse of lower triangular matrix
*
            NN = ( ( N-1 ) / NB )*NB + 1
            DO 30 J = NN, 1, -NB
               JB = MIN( NB, N-J+1 )
               IF( J+JB.LE.N ) THEN
*
*                 Compute rows j+jb:n of current block column
*
                  CALL DTRMM( 'Left', 'Lower', 'No transpose', DIAG,
     $                        N-J-JB+1, JB, ONE, A( J+JB, J+JB ), LDA,
     $                        A( J+JB, J ), LDA )
                  CALL DTRSM( 'Right', 'Lower', 'No transpose', DIAG,
     $                        N-J-JB+1, JB, -ONE, A( J, J ), LDA,
     $                        A( J+JB, J ), LDA )
               END IF
*
*              Compute inverse of current diagonal block
*
               CALL DTRTI2( 'Lower', DIAG, JB, A( J, J ), LDA, INFO )
   30       CONTINUE
         END IF
      END IF
*
      RETURN
*
*     End of DTRTRI
*
      END
      SUBROUTINE DTRTI2( UPLO, DIAG, N, A, LDA, INFO )
*
*  -- LAPACK routine (version 2.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     February 29, 1992
*
*     .. Scalar Arguments ..
      CHARACTER          DIAG, UPLO
      INTEGER            INFO, LDA, N
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * )
*     ..
*
*  Purpose
*  =======
*
*  DTRTI2 computes the inverse of a real upper or lower triangular
*  matrix.
*
*  This is the Level 2 BLAS version of the algorithm.
*
*  Arguments
*  =========
*
*  UPLO    (input) CHARACTER*1
*          Specifies whether the matrix A is upper or lower triangular.
*          = 'U':  Upper triangular
*          = 'L':  Lower triangular
*
*  DIAG    (input) CHARACTER*1
*          Specifies whether or not the matrix A is unit triangular.
*          = 'N':  Non-unit triangular
*          = 'U':  Unit triangular
*
*  N       (input) INTEGER
*          The order of the matrix A.  N >= 0.
*
*  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
*          On entry, the triangular matrix A.  If UPLO = 'U', the
*          leading n by n upper triangular part of the array A contains
*          the upper triangular matrix, and the strictly lower
*          triangular part of A is not referenced.  If UPLO = 'L', the
*          leading n by n lower triangular part of the array A contains
*          the lower triangular matrix, and the strictly upper
*          triangular part of A is not referenced.  If DIAG = 'U', the
*          diagonal elements of A are also not referenced and are
*          assumed to be 1.
*
*          On exit, the (triangular) inverse of the original matrix, in
*          the same storage format.
*
*  LDA     (input) INTEGER
*          The leading dimension of the array A.  LDA >= max(1,N).
*
*  INFO    (output) INTEGER
*          = 0: successful exit
*          < 0: if INFO = -k, the k-th argument had an illegal value
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ONE
      PARAMETER          ( ONE = 1.0D+0 )
*     ..
*     .. Local Scalars ..
      LOGICAL            NOUNIT, UPPER
      INTEGER            J
      DOUBLE PRECISION   AJJ
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
*     ..
*     .. External Subroutines ..
      EXTERNAL           DSCAL, DTRMV, XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters.
*
      INFO = 0
      UPPER = LSAME( UPLO, 'U' )
      NOUNIT = LSAME( DIAG, 'N' )
      IF( .NOT.UPPER .AND. .NOT.LSAME( UPLO, 'L' ) ) THEN
         INFO = -1
      ELSE IF( .NOT.NOUNIT .AND. .NOT.LSAME( DIAG, 'U' ) ) THEN
         INFO = -2
      ELSE IF( N.LT.0 ) THEN
         INFO = -3
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -5
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DTRTI2', -INFO )
         RETURN
      END IF
*
      IF( UPPER ) THEN
*
*        Compute inverse of upper triangular matrix.
*
         DO 10 J = 1, N
            IF( NOUNIT ) THEN
               A( J, J ) = ONE / A( J, J )
               AJJ = -A( J, J )
            ELSE
               AJJ = -ONE
            END IF
*
*           Compute elements 1:j-1 of j-th column.
*
            CALL DTRMV( 'Upper', 'No transpose', DIAG, J-1, A, LDA,
     $                  A( 1, J ), 1 )
            CALL DSCAL( J-1, AJJ, A( 1, J ), 1 )
   10    CONTINUE
      ELSE
*
*        Compute inverse of lower triangular matrix.
*
         DO 20 J = N, 1, -1
            IF( NOUNIT ) THEN
               A( J, J ) = ONE / A( J, J )
               AJJ = -A( J, J )
            ELSE
               AJJ = -ONE
            END IF
            IF( J.LT.N ) THEN
*
*              Compute elements j+1:n of j-th column.
*
               CALL DTRMV( 'Lower', 'No transpose', DIAG, N-J,
     $                     A( J+1, J+1 ), LDA, A( J+1, J ), 1 )
               CALL DSCAL( N-J, AJJ, A( J+1, J ), 1 )
            END IF
   20    CONTINUE
      END IF
*
      RETURN
*
*     End of DTRTI2
*
      END
      INTEGER          FUNCTION ILAENV( ISPEC, NAME, OPTS, N1, N2, N3,
     $                 N4 )
*
*  -- LAPACK auxiliary routine (version 2.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     September 30, 1994
*
*     .. Scalar Arguments ..
      CHARACTER*( * )    NAME, OPTS
      INTEGER            ISPEC, N1, N2, N3, N4
*     ..
*
*  Purpose
*  =======
*
*  ILAENV is called from the LAPACK routines to choose problem-dependent
*  parameters for the local environment.  See ISPEC for a description of
*  the parameters.
*
*  This version provides a set of parameters which should give good,
*  but not optimal, performance on many of the currently available
*  computers.  Users are encouraged to modify this subroutine to set
*  the tuning parameters for their particular machine using the option
*  and problem size information in the arguments.
*
*  This routine will not function correctly if it is converted to all
*  lower case.  Converting it to all upper case is allowed.
*
*  Arguments
*  =========
*
*  ISPEC   (input) INTEGER
*          Specifies the parameter to be returned as the value of
*          ILAENV.
*          = 1: the optimal blocksize; if this value is 1, an unblocked
*               algorithm will give the best performance.
*          = 2: the minimum block size for which the block routine
*               should be used; if the usable block size is less than
*               this value, an unblocked routine should be used.
*          = 3: the crossover point (in a block routine, for N less
*               than this value, an unblocked routine should be used)
*          = 4: the number of shifts, used in the nonsymmetric
*               eigenvalue routines
*          = 5: the minimum column dimension for blocking to be used;
*               rectangular blocks must have dimension at least k by m,
*               where k is given by ILAENV(2,...) and m by ILAENV(5,...)
*          = 6: the crossover point for the SVD (when reducing an m by n
*               matrix to bidiagonal form, if max(m,n)/min(m,n) exceeds
*               this value, a QR factorization is used first to reduce
*               the matrix to a triangular form.)
*          = 7: the number of processors
*          = 8: the crossover point for the multishift QR and QZ methods
*               for nonsymmetric eigenvalue problems.
*
*  NAME    (input) CHARACTER*(*)
*          The name of the calling subroutine, in either upper case or
*          lower case.
*
*  OPTS    (input) CHARACTER*(*)
*          The character options to the subroutine NAME, concatenated
*          into a single character string.  For example, UPLO = 'U',
*          TRANS = 'T', and DIAG = 'N' for a triangular routine would
*          be specified as OPTS = 'UTN'.
*
*  N1      (input) INTEGER
*  N2      (input) INTEGER
*  N3      (input) INTEGER
*  N4      (input) INTEGER
*          Problem dimensions for the subroutine NAME; these may not all
*          be required.
*
* (ILAENV) (output) INTEGER
*          >= 0: the value of the parameter specified by ISPEC
*          < 0:  if ILAENV = -k, the k-th argument had an illegal value.
*
*  Further Details
*  ===============
*
*  The following conventions have been used when calling ILAENV from the
*  LAPACK routines:
*  1)  OPTS is a concatenation of all of the character options to
*      subroutine NAME, in the same order that they appear in the
*      argument list for NAME, even if they are not used in determining
*      the value of the parameter specified by ISPEC.
*  2)  The problem dimensions N1, N2, N3, N4 are specified in the order
*      that they appear in the argument list for NAME.  N1 is used
*      first, N2 second, and so on, and unused problem dimensions are
*      passed a value of -1.
*  3)  The parameter value returned by ILAENV is checked for validity in
*      the calling subroutine.  For example, ILAENV is used to retrieve
*      the optimal blocksize for STRTRI as follows:
*
*      NB = ILAENV( 1, 'STRTRI', UPLO // DIAG, N, -1, -1, -1 )
*      IF( NB.LE.1 ) NB = MAX( 1, N )
*
*  =====================================================================
*
*     .. Local Scalars ..
      LOGICAL            CNAME, SNAME
      CHARACTER*1        C1
      CHARACTER*2        C2, C4
      CHARACTER*3        C3
      CHARACTER*6        SUBNAM
      INTEGER            I, IC, IZ, NB, NBMIN, NX
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          CHAR, ICHAR, INT, MIN, REAL
*     ..
*     .. Executable Statements ..
*
      GO TO ( 100, 100, 100, 400, 500, 600, 700, 800 ) ISPEC
*
*     Invalid value for ISPEC
*
      ILAENV = -1
      RETURN
*
  100 CONTINUE
*
*     Convert NAME to upper case if the first character is lower case.
*
      ILAENV = 1
      SUBNAM = NAME
      IC = ICHAR( SUBNAM( 1:1 ) )
      IZ = ICHAR( 'Z' )
      IF( IZ.EQ.90 .OR. IZ.EQ.122 ) THEN
*
*        ASCII character set
*
         IF( IC.GE.97 .AND. IC.LE.122 ) THEN
            SUBNAM( 1:1 ) = CHAR( IC-32 )
            DO 10 I = 2, 6
               IC = ICHAR( SUBNAM( I:I ) )
               IF( IC.GE.97 .AND. IC.LE.122 )
     $            SUBNAM( I:I ) = CHAR( IC-32 )
   10       CONTINUE
         END IF
*
      ELSE IF( IZ.EQ.233 .OR. IZ.EQ.169 ) THEN
*
*        EBCDIC character set
*
         IF( ( IC.GE.129 .AND. IC.LE.137 ) .OR.
     $       ( IC.GE.145 .AND. IC.LE.153 ) .OR.
     $       ( IC.GE.162 .AND. IC.LE.169 ) ) THEN
            SUBNAM( 1:1 ) = CHAR( IC+64 )
            DO 20 I = 2, 6
               IC = ICHAR( SUBNAM( I:I ) )
               IF( ( IC.GE.129 .AND. IC.LE.137 ) .OR.
     $             ( IC.GE.145 .AND. IC.LE.153 ) .OR.
     $             ( IC.GE.162 .AND. IC.LE.169 ) )
     $            SUBNAM( I:I ) = CHAR( IC+64 )
   20       CONTINUE
         END IF
*
      ELSE IF( IZ.EQ.218 .OR. IZ.EQ.250 ) THEN
*
*        Prime machines:  ASCII+128
*
         IF( IC.GE.225 .AND. IC.LE.250 ) THEN
            SUBNAM( 1:1 ) = CHAR( IC-32 )
            DO 30 I = 2, 6
               IC = ICHAR( SUBNAM( I:I ) )
               IF( IC.GE.225 .AND. IC.LE.250 )
     $            SUBNAM( I:I ) = CHAR( IC-32 )
   30       CONTINUE
         END IF
      END IF
*
      C1 = SUBNAM( 1:1 )
      SNAME = C1.EQ.'S' .OR. C1.EQ.'D'
      CNAME = C1.EQ.'C' .OR. C1.EQ.'Z'
      IF( .NOT.( CNAME .OR. SNAME ) )
     $   RETURN
      C2 = SUBNAM( 2:3 )
      C3 = SUBNAM( 4:6 )
      C4 = C3( 2:3 )
*
      GO TO ( 110, 200, 300 ) ISPEC
*
  110 CONTINUE
*
*     ISPEC = 1:  block size
*
*     In these examples, separate code is provided for setting NB for
*     real and complex.  We assume that NB will take the same value in
*     single or double precision.
*
      NB = 1
*
      IF( C2.EQ.'GE' ) THEN
         IF( C3.EQ.'TRF' ) THEN
            IF( SNAME ) THEN
               NB = 64
            ELSE
               NB = 64
            END IF
         ELSE IF( C3.EQ.'QRF' .OR. C3.EQ.'RQF' .OR. C3.EQ.'LQF' .OR.
     $            C3.EQ.'QLF' ) THEN
            IF( SNAME ) THEN
               NB = 32
            ELSE
               NB = 32
            END IF
         ELSE IF( C3.EQ.'HRD' ) THEN
            IF( SNAME ) THEN
               NB = 32
            ELSE
               NB = 32
            END IF
         ELSE IF( C3.EQ.'BRD' ) THEN
            IF( SNAME ) THEN
               NB = 32
            ELSE
               NB = 32
            END IF
         ELSE IF( C3.EQ.'TRI' ) THEN
            IF( SNAME ) THEN
               NB = 64
            ELSE
               NB = 64
            END IF
         END IF
      ELSE IF( C2.EQ.'PO' ) THEN
         IF( C3.EQ.'TRF' ) THEN
            IF( SNAME ) THEN
               NB = 64
            ELSE
               NB = 64
            END IF
         END IF
      ELSE IF( C2.EQ.'SY' ) THEN
         IF( C3.EQ.'TRF' ) THEN
            IF( SNAME ) THEN
               NB = 64
            ELSE
               NB = 64
            END IF
         ELSE IF( SNAME .AND. C3.EQ.'TRD' ) THEN
            NB = 1
         ELSE IF( SNAME .AND. C3.EQ.'GST' ) THEN
            NB = 64
         END IF
      ELSE IF( CNAME .AND. C2.EQ.'HE' ) THEN
         IF( C3.EQ.'TRF' ) THEN
            NB = 64
         ELSE IF( C3.EQ.'TRD' ) THEN
            NB = 1
         ELSE IF( C3.EQ.'GST' ) THEN
            NB = 64
         END IF
      ELSE IF( SNAME .AND. C2.EQ.'OR' ) THEN
         IF( C3( 1:1 ).EQ.'G' ) THEN
            IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR.
     $          C4.EQ.'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR.
     $          C4.EQ.'BR' ) THEN
               NB = 32
            END IF
         ELSE IF( C3( 1:1 ).EQ.'M' ) THEN
            IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR.
     $          C4.EQ.'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR.
     $          C4.EQ.'BR' ) THEN
               NB = 32
            END IF
         END IF
      ELSE IF( CNAME .AND. C2.EQ.'UN' ) THEN
         IF( C3( 1:1 ).EQ.'G' ) THEN
            IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR.
     $          C4.EQ.'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR.
     $          C4.EQ.'BR' ) THEN
               NB = 32
            END IF
         ELSE IF( C3( 1:1 ).EQ.'M' ) THEN
            IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR.
     $          C4.EQ.'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR.
     $          C4.EQ.'BR' ) THEN
               NB = 32
            END IF
         END IF
      ELSE IF( C2.EQ.'GB' ) THEN
         IF( C3.EQ.'TRF' ) THEN
            IF( SNAME ) THEN
               IF( N4.LE.64 ) THEN
                  NB = 1
               ELSE
                  NB = 32
               END IF
            ELSE
               IF( N4.LE.64 ) THEN
                  NB = 1
               ELSE
                  NB = 32
               END IF
            END IF
         END IF
      ELSE IF( C2.EQ.'PB' ) THEN
         IF( C3.EQ.'TRF' ) THEN
            IF( SNAME ) THEN
               IF( N2.LE.64 ) THEN
                  NB = 1
               ELSE
                  NB = 32
               END IF
            ELSE
               IF( N2.LE.64 ) THEN
                  NB = 1
               ELSE
                  NB = 32
               END IF
            END IF
         END IF
      ELSE IF( C2.EQ.'TR' ) THEN
         IF( C3.EQ.'TRI' ) THEN
            IF( SNAME ) THEN
               NB = 64
            ELSE
               NB = 64
            END IF
         END IF
      ELSE IF( C2.EQ.'LA' ) THEN
         IF( C3.EQ.'UUM' ) THEN
            IF( SNAME ) THEN
               NB = 64
            ELSE
               NB = 64
            END IF
         END IF
      ELSE IF( SNAME .AND. C2.EQ.'ST' ) THEN
         IF( C3.EQ.'EBZ' ) THEN
            NB = 1
         END IF
      END IF
      ILAENV = NB
      RETURN
*
  200 CONTINUE
*
*     ISPEC = 2:  minimum block size
*
      NBMIN = 2
      IF( C2.EQ.'GE' ) THEN
         IF( C3.EQ.'QRF' .OR. C3.EQ.'RQF' .OR. C3.EQ.'LQF' .OR.
     $       C3.EQ.'QLF' ) THEN
            IF( SNAME ) THEN
               NBMIN = 2
            ELSE
               NBMIN = 2
            END IF
         ELSE IF( C3.EQ.'HRD' ) THEN
            IF( SNAME ) THEN
               NBMIN = 2
            ELSE
               NBMIN = 2
            END IF
         ELSE IF( C3.EQ.'BRD' ) THEN
            IF( SNAME ) THEN
               NBMIN = 2
            ELSE
               NBMIN = 2
            END IF
         ELSE IF( C3.EQ.'TRI' ) THEN
            IF( SNAME ) THEN
               NBMIN = 2
            ELSE
               NBMIN = 2
            END IF
         END IF
      ELSE IF( C2.EQ.'SY' ) THEN
         IF( C3.EQ.'TRF' ) THEN
            IF( SNAME ) THEN
               NBMIN = 8
            ELSE
               NBMIN = 8
            END IF
         ELSE IF( SNAME .AND. C3.EQ.'TRD' ) THEN
            NBMIN = 2
         END IF
      ELSE IF( CNAME .AND. C2.EQ.'HE' ) THEN
         IF( C3.EQ.'TRD' ) THEN
            NBMIN = 2
         END IF
      ELSE IF( SNAME .AND. C2.EQ.'OR' ) THEN
         IF( C3( 1:1 ).EQ.'G' ) THEN
            IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR.
     $          C4.EQ.'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR.
     $          C4.EQ.'BR' ) THEN
               NBMIN = 2
            END IF
         ELSE IF( C3( 1:1 ).EQ.'M' ) THEN
            IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR.
     $          C4.EQ.'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR.
     $          C4.EQ.'BR' ) THEN
               NBMIN = 2
            END IF
         END IF
      ELSE IF( CNAME .AND. C2.EQ.'UN' ) THEN
         IF( C3( 1:1 ).EQ.'G' ) THEN
            IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR.
     $          C4.EQ.'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR.
     $          C4.EQ.'BR' ) THEN
               NBMIN = 2
            END IF
         ELSE IF( C3( 1:1 ).EQ.'M' ) THEN
            IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR.
     $          C4.EQ.'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR.
     $          C4.EQ.'BR' ) THEN
               NBMIN = 2
            END IF
         END IF
      END IF
      ILAENV = NBMIN
      RETURN
*
  300 CONTINUE
*
*     ISPEC = 3:  crossover point
*
      NX = 0
      IF( C2.EQ.'GE' ) THEN
         IF( C3.EQ.'QRF' .OR. C3.EQ.'RQF' .OR. C3.EQ.'LQF' .OR.
     $       C3.EQ.'QLF' ) THEN
            IF( SNAME ) THEN
               NX = 128
            ELSE
               NX = 128
            END IF
         ELSE IF( C3.EQ.'HRD' ) THEN
            IF( SNAME ) THEN
               NX = 128
            ELSE
               NX = 128
            END IF
         ELSE IF( C3.EQ.'BRD' ) THEN
            IF( SNAME ) THEN
               NX = 128
            ELSE
               NX = 128
            END IF
         END IF
      ELSE IF( C2.EQ.'SY' ) THEN
         IF( SNAME .AND. C3.EQ.'TRD' ) THEN
            NX = 1
         END IF
      ELSE IF( CNAME .AND. C2.EQ.'HE' ) THEN
         IF( C3.EQ.'TRD' ) THEN
            NX = 1
         END IF
      ELSE IF( SNAME .AND. C2.EQ.'OR' ) THEN
         IF( C3( 1:1 ).EQ.'G' ) THEN
            IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR.
     $          C4.EQ.'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR.
     $          C4.EQ.'BR' ) THEN
               NX = 128
            END IF
         END IF
      ELSE IF( CNAME .AND. C2.EQ.'UN' ) THEN
         IF( C3( 1:1 ).EQ.'G' ) THEN
            IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR.
     $          C4.EQ.'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR.
     $          C4.EQ.'BR' ) THEN
               NX = 128
            END IF
         END IF
      END IF
      ILAENV = NX
      RETURN
*
  400 CONTINUE
*
*     ISPEC = 4:  number of shifts (used by xHSEQR)
*
      ILAENV = 6
      RETURN
*
  500 CONTINUE
*
*     ISPEC = 5:  minimum column dimension (not used)
*
      ILAENV = 2
      RETURN
*
  600 CONTINUE 
*
*     ISPEC = 6:  crossover point for SVD (used by xGELSS and xGESVD)
*
      ILAENV = INT( REAL( MIN( N1, N2 ) )*1.6E0 )
      RETURN
*
  700 CONTINUE
*
*     ISPEC = 7:  number of processors (not used)
*
      ILAENV = 1
      RETURN
*
  800 CONTINUE
*
*     ISPEC = 8:  crossover point for multishift (used by xHSEQR)
*
      ILAENV = 50
      RETURN
*
*     End of ILAENV
*
      END
      LOGICAL          FUNCTION LSAME( CA, CB )
*
*  -- LAPACK auxiliary routine (version 2.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     September 30, 1994
*
*     .. Scalar Arguments ..
      CHARACTER          CA, CB
*     ..
*
*  Purpose
*  =======
*
*  LSAME returns .TRUE. if CA is the same letter as CB regardless of
*  case.
*
*  Arguments
*  =========
*
*  CA      (input) CHARACTER*1
*  CB      (input) CHARACTER*1
*          CA and CB specify the single characters to be compared.
*
* =====================================================================
*
*     .. Intrinsic Functions ..
      INTRINSIC          ICHAR
*     ..
*     .. Local Scalars ..
      INTEGER            INTA, INTB, ZCODE
*     ..
*     .. Executable Statements ..
*
*     Test if the characters are equal
*
      LSAME = CA.EQ.CB
      IF( LSAME )
     $   RETURN
*
*     Now test for equivalence if both characters are alphabetic.
*
      ZCODE = ICHAR( 'Z' )
*
*     Use 'Z' rather than 'A' so that ASCII can be detected on Prime
*     machines, on which ICHAR returns a value with bit 8 set.
*     ICHAR('A') on Prime machines returns 193 which is the same as
*     ICHAR('A') on an EBCDIC machine.
*
      INTA = ICHAR( CA )
      INTB = ICHAR( CB )
*
      IF( ZCODE.EQ.90 .OR. ZCODE.EQ.122 ) THEN
*
*        ASCII is assumed - ZCODE is the ASCII code of either lower or
*        upper case 'Z'.
*
         IF( INTA.GE.97 .AND. INTA.LE.122 ) INTA = INTA - 32
         IF( INTB.GE.97 .AND. INTB.LE.122 ) INTB = INTB - 32
*
      ELSE IF( ZCODE.EQ.233 .OR. ZCODE.EQ.169 ) THEN
*
*        EBCDIC is assumed - ZCODE is the EBCDIC code of either lower or
*        upper case 'Z'.
*
         IF( INTA.GE.129 .AND. INTA.LE.137 .OR.
     $       INTA.GE.145 .AND. INTA.LE.153 .OR.
     $       INTA.GE.162 .AND. INTA.LE.169 ) INTA = INTA + 64
         IF( INTB.GE.129 .AND. INTB.LE.137 .OR.
     $       INTB.GE.145 .AND. INTB.LE.153 .OR.
     $       INTB.GE.162 .AND. INTB.LE.169 ) INTB = INTB + 64
*
      ELSE IF( ZCODE.EQ.218 .OR. ZCODE.EQ.250 ) THEN
*
*        ASCII is assumed, on Prime machines - ZCODE is the ASCII code
*        plus 128 of either lower or upper case 'Z'.
*
         IF( INTA.GE.225 .AND. INTA.LE.250 ) INTA = INTA - 32
         IF( INTB.GE.225 .AND. INTB.LE.250 ) INTB = INTB - 32
      END IF
      LSAME = INTA.EQ.INTB
*
*     RETURN
*
*     End of LSAME
*
      END
