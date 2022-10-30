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
