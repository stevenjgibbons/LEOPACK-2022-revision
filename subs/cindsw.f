C*********************************************************************
C subroutine Curl INDex SWitch ***************************************
C            -    ---   --     ***************************************
C Steve Gibbons Sat Sep 25 15:41:55 BST 1999                         C
C____________________________________________________________________C
C                                                                    C
C This routine is for programming conveience rather than scientific  C
C necessity. The integer array MHT has length NH and for each        C
C harmonic, ih, in the solution vector, MHT( ih ) corresponds to     C
C itype - which has the following options.                           C
C                                                                    C
C         ITYPE = 1 for a poloidal velocity harmonic.                C
C         ITYPE = 2 for a toroidal velocity harmonic.                C
C         ITYPE = 3 for a temperature harmonic.                      C
C         ITYPE = 4 for a poloidal magnetic field harmonic.          C
C         ITYPE = 5 for a toroidal magnetic field harmonic.          C
C                                                                    C
C Now when solving the vorticity equation, the curl is taken of      C
C the momentum equation and so equations for toroidal harmonics and  C
C those for poloidal harmonics are interchanged.                     C
C                                                                    C
C CINDSW simply makes a copy of MHT but with the 1s and the 2s       C
C interchanged. This can then be supplied to other routines as the   C
C MHT array for the destination (row) function.                      C
C                                                                    C
C____________________________________________________________________C
C                                                                    C
C Input variables :-                                                 C
C ===============                                                    C
C  Integer                                                           C
C  -------                                                           C
C                                                                    C
C     NH        : Number of spherical harmonics.                     C
C     MHTI      : MHTI( ih ) itype for harmonic 'ih'                 C
C     MHTO      : MHTO( ih ) itype for curl of harmonic 'ih'         C
C                  (only effects velocity harmonics.)                C
C                                                                    C
C                                                                    C
C____________________________________________________________________C
C
C*********************************************************************
      SUBROUTINE CINDSW ( NH, MHTI, MHTO )
      IMPLICIT NONE
C____________________________________________________________________C
C Variable declarations - Parameters ................................C
      INTEGER NH, MHTI( NH ), MHTO( NH )
C____________________________________________________________________C
C Variable declarations - Working variables .........................C
      INTEGER IH
C____________________________________________________________________C
C START OF PROGRAM **************************************************C
C____________________________________________________________________C
c
      DO IH = 1, NH
        IF ( MHTI( IH ).EQ.1 ) MHTO( IH ) = 2
        IF ( MHTI( IH ).EQ.2 ) MHTO( IH ) = 1
        IF ( MHTI( IH ).EQ.3 ) MHTO( IH ) = 3
        IF ( MHTI( IH ).EQ.4 ) MHTO( IH ) = 4
        IF ( MHTI( IH ).EQ.5 ) MHTO( IH ) = 5
      ENDDO
c
      RETURN
      END
C*********************************************************************
