C*********************************************************************
C subroutine Integer Index Array Single Component Extract ************
C            -       -     -     -      -         -       ************
C Steve Gibbons Wed Nov  1 08:51:08 MET 2000                         C
C____________________________________________________________________C
C                                                                    C
C Takes as input the large integer arrays defining spherical         C
C harmonic properties and fills equivalent arrays with indices       C
C corresponding to a single component (1,2,3,4 or 5)                 C
C____________________________________________________________________C
C                                                                    C
C Input variables :-                                                 C
C ===============                                                    C
C  Integer                                                           C
C  -------                                                           C
C                                                                    C
C     NH        : Total number of harmonics in full harm. set.       C
C     ICOMP     : Single component identifier for comparison with    C
C                  MHT array. Must be one of the following:-         C
C                                                                    C
C         ICOMP = 1 for a poloidal velocity harmonic.                C
C         ICOMP = 2 for a toroidal velocity harmonic.                C
C         ICOMP = 3 for a temperature harmonic.                      C
C         ICOMP = 4 for a poloidal magnetic field harmonic.          C
C         ICOMP = 5 for a toroidal magnetic field harmonic.          C
C                                                                    C
C     MHT       : MHT( ih ) contains itype for harmonic 'ih'         C
C                                                                    C
C     MHL       : MHL( ih ) contains degree, l, for harmonic 'ih'    C
C                                                                    C
C     MHM       : MHM( ih ) contains order, m, for harmonic 'ih' if  C
C                  ih has cos m phi dependency and (-m) if ih has    C
C                   sin m phi dependency.                            C
C                                                                    C
C     MHP       : MHP( ih ) contains is, finite difference scheme    C
C                  identifier.                                       C
C                                                                    C
C     NCH       : (Output) number of harmonics in single component   C
C                   set,                                             C
C                                                                    C
C     NCHMAX    : Maximum value for NCH.                             C
C                                                                    C
C     MLC       : Equiv. of MHL for single component set.            C
C     MMC       : Equiv. of MHL for single component set.            C
C     MPC       : Equiv. of MHL for single component set.            C
C                                                                    C
C____________________________________________________________________C
C
C*********************************************************************
      SUBROUTINE IIASCE( NH, ICOMP, MHT, MHL, MHM, MHP, NCH, NCHMAX,
     1                   MLC, MMC, MPC )
      IMPLICIT NONE
C____________________________________________________________________C
C Variable declarations - Parameters ................................C
      INTEGER NH, ICOMP, MHT( * ), MHL( * ), MHM( * ), MHP( * ), NCH,
     1        NCHMAX, MLC( * ), MMC( * ), MPC( * )
C____________________________________________________________________C
C Variable declarations - Working variables .........................C
      INTEGER IH
      LOGICAL OK, KO
C____________________________________________________________________C
C START OF PROGRAM **************************************************C
C____________________________________________________________________C
c
      IF ( ICOMP.LT.1 .OR. ICOMP.GT.5 ) THEN
        PRINT *,' Subroutine IIASCE.'
        PRINT *,' ICOMP = ', ICOMP,' : Invalid value.'
        PRINT *,' Program aborted.'
        STOP
      ENDIF
c
      OK  = .TRUE.
      NCH = 0
      DO IH = 1, NH
        KO = .FALSE.
        IF ( MHT( IH ).EQ.ICOMP ) THEN
          NCH = NCH + 1
          KO  = .TRUE.
        ENDIF
        IF ( NCH.GT.NCHMAX ) OK = .FALSE.
        IF ( KO .AND. OK ) THEN
          MLC( NCH ) = MHL( IH )
          MMC( NCH ) = MHM( IH )
          MPC( NCH ) = MHP( IH )
        ENDIF
      ENDDO
c
      IF ( OK ) RETURN
c
      PRINT *,' Subroutine IIASCE. NCHMAX = ', NCHMAX
      PRINT *,' You require ', NCH,' harmonics for '
      PRINT *,' component ', ICOMP,'. Program aborted.'
      STOP
      END
C*********************************************************************
