C*********************************************************************
C subroutine Spherical Harmonic SELect *******************************
C            -         -        ---    *******************************
C Steve Gibbons 22.3.99                                              C
C____________________________________________________________________C
C Input variables :-                                                 C
C ===============                                                    C
C  Integer                                                           C
C  -------                                                           C
C     LH        : Maximum degree of spherical harmonic               C
C     ISYM      : Equatorial symmetry select flag.                   C
C                  ISYM = 1 - select only equatorial symmetric harms C
C                  ISYM = 2 - select only equatorial anti-sym harms  C
C                  ISYM = 3 - select both EA and ES harmonics        C
C     NMAX      : Maximum number of harmonics permitted.             C
C                 When more than NMAX harmonics have been selected,  C
C                 the requirements are still calculated and the      C
C                 program is aborted on completion of the routine    C
C                 with a message indicating how many harmonics are   C
C                 infact required for the requirements.              C
C     NMODES    : Number of angular modes required.                  C
C     MMODE     : Array dimension ( NMODES ).                        C
C     IFLAG     : Set to 1 to include poloidal velocity, toroidal    C
C                 velocity and temperature harmonics.                C
C                 only option as of 22.3.99                          C
C                                                                    C
C____________________________________________________________________C
C Output variables :-                                                C
C ================                                                   C
C  Integer                                                           C
C  -------                                                           C
C                                                                    C
C     MHT       : Dimension ( NH )                                   C
C     MHL       : Dimension ( NH )                                   C
C     MHM       : Dimension ( NH )                                   C
C     MHC       : Dimension ( NH )                                   C
C                                                                    C
C  mht, mhl, mhm and mhc define the spherical harmonics present      C
C  in the solution vector. For spherical harmonic number I;          C
C                                                                    C
C   MHT( I ) = 1 for a poloidal velocity vector                      C
C   MHT( I ) = 2 for a toroidal velocity vector                      C
C   MHT( I ) = 3 for a temperature / codensity term                  C
C   MHT( I ) = 4 for a poloidal magnetic field vector                C
C   MHT( I ) = 5 for a poloidal magnetic field vector                C
C                                                                    C
C   MHL( I ) = spherical harmonic degree, l                          C
C                                                                    C
C   MHM( I ) = spherical harmonic order, m                           C
C                                                                    C
C   MHC( I ) = 1 for a cosine dependence in phi and                  C
C            = 2  "  "  sine     "        "  "                       C
C                                                                    C
C     NH        : Total number of harmonics used.                    C
C____________________________________________________________________C
C
C*********************************************************************
      SUBROUTINE SHSEL ( LH, ISYM, NMAX, NMODES, MMODE, IFLAG,
     1                   MHT, MHL, MHM, MHC, NH )
      IMPLICIT NONE
C____________________________________________________________________C
C Variable declarations - Parameters ................................C
      INTEGER LH, ISYM, NMAX, NMODES, MMODE( NMODES ), IFLAG,
     2        MHT( * ), MHL( * ), MHM( * ), MHC( * ), NH
C____________________________________________________________________C
C Variable declarations - Working variables .........................C
      INTEGER L, M, ICS, IM, IPTT
      LOGICAL OK, OSEL, OEA, OES, ODDEVN, OE
C____________________________________________________________________C
C START OF PROGRAM **************************************************C
C____________________________________________________________________C
c
      NH = 0
c
      OES = .FALSE.
      OEA = .FALSE.
      IF ( ISYM.EQ.1 .OR. ISYM.EQ.3 ) OES = .TRUE.
      IF ( ISYM.EQ.2 .OR. ISYM.EQ.3 ) OEA = .TRUE.
c
      IF ( NMODES.LT.1 ) RETURN
c
      IF ( IFLAG.NE.1 ) THEN
         PRINT *,' Subroutine SHSEL. IFLAG = ', IFLAG
         PRINT *,' Currently, IFLAG must be 1. '
         STOP
      ENDIF
c
      IF ( ISYM.NE.1 .AND. ISYM.NE.2 .AND. ISYM.NE.3 ) THEN
         PRINT *,' Subroutine SHSEL. ISYM = ', ISYM
         PRINT *,' Currently, ISYM must be 1, 2 or 3.'
         STOP
      ENDIF
c
      OK = .TRUE.
      IF ( IFLAG.EQ.1 ) THEN
         DO IM = 1, NMODES
           M = MMODE( IM )
c          . Handle the monopole temperature term
           IF ( M.EQ.0 .AND. OES ) THEN
             NH = NH + 1
             IF ( NH.GT.NMAX ) OK = .FALSE.
             IF ( OK ) MHT( NH ) = 3
             IF ( OK ) MHL( NH ) = 0
             IF ( OK ) MHM( NH ) = 0
             IF ( OK ) MHC( NH ) = 1
           ENDIF
           DO IPTT = 1, 3
            DO L = M, LH
             OE = ODDEVN ( L-M )
             IF ( L.EQ.0 ) GOTO 501
             DO ICS = 1, 2
               IF ( M.EQ.0 .AND. ICS.EQ.2 ) GOTO 500
               OSEL = .FALSE.
               IF ( OE.AND.OES.AND.(IPTT.EQ.1) ) OSEL = .TRUE.
               IF ( (.NOT. OE).AND.OEA.AND.(IPTT.EQ.1) ) OSEL = .TRUE.
               IF ( OE.AND.OEA.AND.(IPTT.EQ.2) ) OSEL = .TRUE.
               IF ( (.NOT. OE).AND.OES.AND.(IPTT.EQ.2) ) OSEL = .TRUE.
               IF ( OE.AND.OES.AND.(IPTT.EQ.3) ) OSEL = .TRUE.
               IF ( (.NOT. OE).AND.OEA.AND.(IPTT.EQ.3) ) OSEL = .TRUE.
               IF ( OSEL ) NH = NH + 1
               IF ( NH.GT.NMAX ) OK = .FALSE.
               IF ( OSEL .AND. OK ) MHT( NH ) = IPTT
               IF ( OSEL .AND. OK ) MHL( NH ) = L
               IF ( OSEL .AND. OK ) MHM( NH ) = M
               IF ( OSEL .AND. OK ) MHC( NH ) = ICS
 500           CONTINUE
             ENDDO
 501         CONTINUE
            ENDDO
           ENDDO
         ENDDO
      ENDIF
      IF ( OK ) RETURN
      PRINT *,' Subroutine SHSEL. Your specifications'
      PRINT *,' require ',NH,' harmonics. Maximum was set'
      PRINT *,' at ',NMAX,'. Program aborted.'
c
      STOP
      END
C*********************************************************************

