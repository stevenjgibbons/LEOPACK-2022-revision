C*********************************************************************
C                                                                    C
C Steve Gibbons -                                                    C
C                                                                    C
C Solution Vector ENergy SPECtrum calculate                          C
C -        -      --     ----                                        C
C Steve Gibbons Mon Mar  6 17:13:05 GMT 2000                         C
C                                                                    C
C Reads in a .xarr, .vecs and .ints file and outputs a file          C
C containing the energy integral spectra and totals.                 C
C                                                                    C
C*********************************************************************
      PROGRAM svenspec
      IMPLICIT NONE
C_____________________________________________________________________
C                                                                    C
C Variable declarations - array defining parameters.                 C
C                                                                    C
C____________________________________________________________________C
C
      INTEGER NRMAX, NHMAX, NDCS, LHMAX, NBN, NCFM, LHLH2M
      PARAMETER ( NRMAX = 200, NHMAX = 5000, LHMAX = 156,
     1            NDCS = LHMAX + 4, NBN = 3, NCFM = 2*NBN + 1,
     2            LHLH2M = LHMAX*(LHMAX+2)  )
C_____________________________________________________________________
C                                                                    C
C Variable declarations - main arrays.                               C
C                                                                    C
C____________________________________________________________________C
C
      INTEGER MHT( NHMAX ), MHL( NHMAX ), MHM( NHMAX ),
     1        MHP( NHMAX ), MHIBC( NDCS ),
     2        MHOBC( NDCS ), LARR( NDCS ), IWORK( NCFM )
      DOUBLE PRECISION XARR( NRMAX ), SV( NRMAX*NHMAX ), 
     1                 SVFDC( NCFM, NRMAX, 2, NDCS ),
     2                 COEFM1( NCFM, NCFM ), COEF1( NCFM ),
     3                 COEFM2( NCFM, NCFM ), COEF2( NCFM )
      DOUBLE PRECISION PVKEA( LHLH2M ), TVKEA( LHLH2M ),
     1                 PFMEA( LHLH2M ), TFMEA( LHLH2M ),
     2                 SPECL( 4, LHMAX ), SPECM( 4, 0:LHMAX )
C_____________________________________________________________________
C                                                                    C
C Variable declarations - scalars and small arrays                   C
C                                                                    C
C____________________________________________________________________C
C
      INTEGER LU, INARR( 3 ), NH, I, IWR, L, M, ICS, IOP,
     1        NR, NCUDS, NDRVS, ILN, IRN, INDSHC, IND, LH
      DOUBLE PRECISION RI, RO, PI, ZERO, PFME, PVKE, TOTKE, TVKE,
     1                 TFME, TOTME, DKE( 2 ), VOLFAC, RTOTKE,
     2                 RTOTME, TOTAL, DLOW, EN1, EN2, EN3, EN4
      CHARACTER *(80) FNAME
C_____________________________________________________________________
C                                                                    C
C Variable declarations - parameters and EXTERNAL declarations.      C
C                                                                    C
C____________________________________________________________________C
C
      PARAMETER ( PI=3.14159265358979312D0, ZERO = 0.0d0, IOP = 0,
     1            DLOW = 1.0d-9 )
C_____________________________________________________________________
C                  **************************************************C
C START OF PROGRAM **************************************************C
C____________________________________________________________________C
C
      I = LHLH2M
      CALL VECOP( PVKEA, ZERO, I, IOP )
      CALL VECOP( TVKEA, ZERO, I, IOP )
      CALL VECOP( PFMEA, ZERO, I, IOP )
      CALL VECOP( TFMEA, ZERO, I, IOP )
C
      I = 4
      CALL MATOP( SPECL, ZERO, I, LHMAX, IOP )
      CALL MATOP( SPECM, ZERO, I, LHMAX+1, IOP )
C
      DO I = 1, 80
        FNAME(I:I) = ' '
      ENDDO
C
 80   FORMAT(A)
      PRINT *,'-------------------------------------------'
      PRINT *,' Enter filename of integer harmonics file.'
      PRINT *,'-------------------------------------------'
      READ ( 5, 80 ) FNAME
      PRINT *,' Integers file = ', FNAME
      LU  = 91
      IWR = 3
c     CALL BIHFRD( NH, NHMAX, MHT, MHL, MHM, LU, FNAME )
      NCUDS = 0
C     (ncuds - number of diff. schemes already in use).
      CALL HMFRD( NH, NHMAX, MHT, MHL, MHM, MHP, NCUDS, NDCS,
     1            MHIBC, MHOBC, LARR, LU, FNAME )
      INARR( 3 ) = NH
C
      PRINT *,'----------------------------------------'
      PRINT *,' Enter filename of radial spacings file.'
      PRINT *,'----------------------------------------'
      READ ( 5, 80 ) FNAME
      PRINT *,' Radial spacings file = ', FNAME
      CALL XARRRD( NR, NRMAX, XARR, LU, FNAME )
C
      PRINT *,'-------------------------------'
      PRINT *,' Enter filename of vector file.'
      PRINT *,'-------------------------------'
      READ ( 5, 80 ) FNAME
      PRINT *,' Vector file = ', FNAME
      CALL SVFRD( INARR, LU, NRMAX, SV, FNAME )
C
      IF ( NR.NE.INARR( 2 ) ) THEN
        PRINT *,' NR         = ', NR
        PRINT *,' INARR( 2 ) = ', INARR( 2 )
        PRINT *,' XARR cannot correspond to SV.'
        PRINT *,' Program aborted.'
        STOP
      ENDIF
C
      RI = XARR( 1  )
      RO = XARR( NR )
C
C----------------------------------------------------------------
C          Calculate array of finite difference coefficients.   |
C----------------------------------------------------------------
C
      NDRVS = 1
      ILN   = 1
      IRN   = NR
      CALL SVFDCF( NR, NDCS, NBN, ILN, IRN, MHIBC, MHOBC,
     1             LARR, NCFM, NCFM, NDRVS, NDRVS, XARR,
     2             IWORK, SVFDC, COEFM1, COEFM2, COEF1, COEF2 )
C
      LH = 0
      DO I = 1, NH
        CALL SHKEER( I, NDCS, NR, INARR, MHT, MHL, MHP, NBN,
     1             NDRVS, NDRVS, NCFM, SV, XARR, DKE, SVFDC )
        L  = MHL( I )
        IF ( L.GT.LHMAX ) THEN
          PRINT *,' Harmonic has been found with L = ', L
          PRINT *,' LHMAX = ', LHMAX
          PRINT *,' Recompile and rerun with larger LHMAX.'
          PRINT *,' Program aborted.'
          STOP
        ENDIF
        IF ( L.GT.LH ) LH = L
        IF ( MHM( I ).LT.0 ) THEN
          M   = -MHM( I )
          ICS = 2
        ELSE
          M   = MHM( I )
          ICS = 1
        ENDIF
        IND = INDSHC( L, M, ICS )
        IF ( MHT( I ).EQ.1 ) THEN
          PVKEA( IND ) = PVKEA( IND ) + DKE( 1 )
          SPECL( 1, L ) = SPECL( 1, L ) + DKE( 1 )
          SPECM( 1, M ) = SPECM( 1, M ) + DKE( 1 )
        ENDIF
        IF ( MHT( I ).EQ.2 ) THEN
          TVKEA( IND ) = TVKEA( IND ) + DKE( 1 )
          SPECL( 2, L ) = SPECL( 2, L ) + DKE( 1 )
          SPECM( 2, M ) = SPECM( 2, M ) + DKE( 1 )
        ENDIF
        IF ( MHT( I ).EQ.4 ) THEN
          PFMEA( IND ) = PFMEA( IND ) + DKE( 1 )
          SPECL( 3, L ) = SPECL( 3, L ) + DKE( 1 )
          SPECM( 3, M ) = SPECM( 3, M ) + DKE( 1 )
        ENDIF
        IF ( MHT( I ).EQ.5 ) THEN
          TFMEA( IND ) = TFMEA( IND ) + DKE( 1 )
          SPECL( 4, L ) = SPECL( 4, L ) + DKE( 1 )
          SPECM( 4, M ) = SPECM( 4, M ) + DKE( 1 )
        ENDIF
      ENDDO
C
C Calculate the totals of kinetic and magnetic energies
C
      PVKE = ZERO
      TVKE = ZERO
      PFME = ZERO
      TFME = ZERO
      DO I = 1, LH*(LH+2)
        PVKE = PVKE + PVKEA( I )
        TVKE = TVKE + TVKEA( I )
        PFME = PFME + PFMEA( I )
        TFME = TFME + TFMEA( I )
      ENDDO
      TOTKE = PVKE + TVKE
      TOTME = PFME + TFME
C
      VOLFAC = 4.0d0*PI*(RO*RO*RO - RI*RI*RI)/3.0d0
C
      RTOTKE = TOTKE/VOLFAC
      RTOTME = TOTME/VOLFAC
C
      PRINT *,'----------------------------------'
      PRINT *,' Enter filename of output file.   '
      PRINT *,'----------------------------------'
      READ ( 5, 80 ) FNAME
      PRINT *,' Output file = ', FNAME
      CALL FOPEN( LU, FNAME, IWR )
      WRITE ( LU, 85 ) PVKE
 85   FORMAT('Poloidal kinetic energy   = ',f20.6)
      WRITE ( LU, 86 ) TVKE
 86   FORMAT('Toroidal kinetic energy   = ',f20.6)
      WRITE ( LU, 87 ) TOTKE
 87   FORMAT('Total kinetic energy      = ',f20.6)
      WRITE ( LU, 88 ) PFME
 88   FORMAT('Poloidal magnetic energy  = ',f20.6)
      WRITE ( LU, 89 ) TFME
 89   FORMAT('Toroidal magnetic energy  = ',f20.6)
      WRITE ( LU, 90 ) TOTME
 90   FORMAT('Total magnetic energy     = ',f20.6)
      WRITE ( LU, 91 ) RTOTKE
 91   FORMAT('Scaled kinetic energy     = ',f20.6)
      WRITE ( LU, 92 ) RTOTME
 92   FORMAT('Scaled magnetic energy    = ',f20.6)
C
      DO L = 1, LH
        TOTAL = ZERO
        EN1   = SPECL( 1, L )
        EN2   = SPECL( 2, L )
        EN3   = SPECL( 3, L )
        EN4   = SPECL( 4, L )
        TOTAL = DABS( EN1 ) + DABS( EN2 ) + DABS( EN3 ) + DABS( EN4 )
        IF ( TOTAL.GE.DLOW ) WRITE ( LU, 117 ) 'L= ', L,
     1                       EN1, EN2, EN3, EN4
      ENDDO
C
      DO M = 0, LH
        TOTAL = ZERO
        EN1   = SPECM( 1, M )
        EN2   = SPECM( 2, M )
        EN3   = SPECM( 3, M )
        EN4   = SPECM( 4, M )
        TOTAL = DABS( EN1 ) + DABS( EN2 ) + DABS( EN3 ) + DABS( EN4 )
        IF ( TOTAL.GE.DLOW ) WRITE ( LU, 117 ) 'M= ', M,
     1                       EN1, EN2, EN3, EN4
      ENDDO
C
 117  FORMAT(A3,I6,1PD16.7,1PD16.7,1PD16.7,1PD16.7)
      CALL FCLOSE( LU, FNAME, 'Error - parameter file' )
      STOP
      END
C*********************************************************************
