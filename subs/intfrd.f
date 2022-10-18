C*********************************************************************
C subroutine INTeger index File ReaD *********************************
C            ---           -    -  - *********************************
C Steve Gibbons 23.6.99                                              C
C____________________________________________________________________C
C                                                                    C
C  Reads in a file (FNAME) and fills the arrays MHT, MHL, MHM and    C
C MHC. The first line contains LH, NH and the following NH lines     C
C contain MHT( i ), MHL( i ), MHM( i ), MHC( i )                     C
C                                                                    C
C                                                                    C
C____________________________________________________________________C
C Input variables :-                                                 C
C ===============                                                    C
C  Integer                                                           C
C  -------                                                           C
C     NHMAX     : Maximum number of harmonics (all types)            C
C                                                                    C
C  Character                                                         C
C  ---------                                                         C
C     FNAME     : File name (unspecified length)                     C
C                                                                    C
C Output variables :-                                                C
C ================                                                   C
C  Integer                                                           C
C  -------                                                           C
C     NH        : Number of harmonics (all types)                    C
C     LH        : Maximum spherical harmonic degree                  C
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
C                                                                    C
C   MHL( I ) = spherical harmonic degree, l                          C
C                                                                    C
C   MHM( I ) = spherical harmonic order, m                           C
C                                                                    C
C   MHC( I ) = 1 for a cosine dependence in phi and                  C
C            = 2  "  "  sine     "        "  "                       C
C                                                                    C
C____________________________________________________________________C
C
C*********************************************************************
      SUBROUTINE INTFRD ( NHMAX, FNAME, NH, LH, MHT, MHL, MHM, MHC )
      IMPLICIT NONE
C____________________________________________________________________C
C Variable declarations - Parameters ................................C
      INTEGER NHMAX, NH, LH, MHT( * ), MHL( * ), MHM( * ), MHC( * )
      CHARACTER *(*) FNAME
C____________________________________________________________________C
C Variable declarations - Working Variables .........................C
      INTEGER LU, IWR, IH
C____________________________________________________________________C
C START OF PROGRAM **************************************************C
C____________________________________________________________________C
C
C Open file for reading
C
      LU = 999
      IWR = 1
      CALL FOPEN ( LU, FNAME, IWR )
C
      READ ( LU, * ) LH, NH
      IF ( NH.GT.NHMAX ) THEN
         PRINT *,' Subroutine INTFRD. NH = ', NH
         PRINT *,' NHMAX = ', NHMAX
         STOP
      ENDIF
C
      DO IH = 1, NH
         READ ( LU, * ) MHT( IH ), MHL( IH ), MHM( IH ), MHC( IH )
      ENDDO
C
C Close file
C
      CALL FCLOSE ( LU, FNAME, 'error' )
C
      RETURN
      END
C*********************************************************************
