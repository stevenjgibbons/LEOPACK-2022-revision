C*********************************************************************
C subroutine Steady Flow Streamline Explicit Time Integration ********
C            -      -    -          -        -    -           ********
C Steve Gibbons Tue May  8 10:17:05 MET DST 2001                     C
C____________________________________________________________________C
C                                                                    C
C A horrifically primitive routine to track with time, the course of C
C a particle in a steady flow given by the vector VEC and its        C
C defining arrays. The format is given by the INARR array            C
C (see INDFUN) harmonics type given by MHT (a harmonic, ih, has      C
C                                                                    C
C      MHT( ih ) = 1 --> poloidal velocity                           C
C      MHT( ih ) = 2 --> toroidal velocity                           C
C      MHT( ih ) = 3 --> temperature                                 C
C      MHT( ih ) = 4 --> poloidal magnetic field                     C
C      MHT( ih ) = 5 --> toroidal magnetic field )                   C
C                                                                    C
C MHL( ih ) is the spherical harmonic degree, l.                     C
C MHM( ih ) is the spherical harmonic oder, m, for cos m phi         C
C           dependence and, -m, for sin m phi dependence.            C
C                                                                    C
C There are NTS time-steps and the time-step size is DELTAT.         C
C                                                                    C
C The double precision array DSPCA has dimensions ( 3, NTS )         C
C with DSPCA( 1, ITS ) containing rad_i, DSPCA( 2, ITS ) containing  C
C the_i and DSPCA( 3, ITS ) containing phi_i.                        C
C                                                                    C
C DSPCA( 1, 1 ), DSPCA( 2, 1 ) and DSPCA( 3, 1 ) must be set on      C
C entry as the initial position.                                     C
C                                                                    C
C____________________________________________________________________C
C                                                                    C
C Input variables :-                                                 C
C ===============                                                    C
C  Integer                                                           C
C  -------                                                           C
C                                                                    C
C     INARR     : Integer array dimension ( * ).                     C
C                  Elements may be arbitrary except for              C
C                  INARR( 1 ) = IFORMF - flag for vector format.     C
C                                        See INDFUN                  C
C                  INARR( 2 ) = NR. Number of radial grid nodes.     C
C                  INARR( 3 ) = NH = total number of radial func.s   C
C     NNDS      : Number of nodes to be used in interpolation.       C
C                 Must be atleast 3.                                 C
C                                                                    C
C     MHT       : Dim (*). See above.                                C
C     MHL       : Dim (*). See above.                                C
C     MHM       : Dim (*). See above.                                C
C                                                                    C
C     NTS       : Number of time-steps to be taken.                  C
C                                                                    C
C  Double Precision                                                  C
C  ----------------                                                  C
C                                                                    C
C     VEC       : Solution vector. Dim ( * ) but length atleast      C
C                  NR*NH.                                            C
C     XARR      : Dim ( * ) but length atleast NR. Location of       C
C                  radial grid nodes.                                C
C                                                                    C
C     DELTAT    : Size of time-step.                                 C
C     DSPCA     : Dim ( 3, NTS )                                     C
C                                                                    C
C____________________________________________________________________C
C
C*********************************************************************
      SUBROUTINE SFSETI( INARR, NNDS, MHT, MHL, MHM, NTS,
     1                   VEC, XARR, DELTAT, DSPCA )
      IMPLICIT NONE
C____________________________________________________________________C
C Variable declarations - Parameters ................................C
      INTEGER          INARR( 3 ), NNDS, MHT( * ), MHL( * ),
     1                 MHM( * ), NTS
      DOUBLE PRECISION VEC( * ), XARR( * ), DELTAT, DSPCA( 3, NTS )
C____________________________________________________________________C
C Variable declarations - Working variables .........................C
      INTEGER          NNDM, ICOMP, ITS
      PARAMETER ( NNDM = 7 )
      INTEGER          IWORK( NNDM )
      DOUBLE PRECISION WORK1( NNDM ), WORK2( NNDM ),
     1                 COEFM( NNDM, NNDM ), RTPFCE,
     2                 RAD, THE, PHI, VRAD, VTHE, VPHI
      EXTERNAL         RTPFCE
C____________________________________________________________________C
C START OF PROGRAM **************************************************C
C____________________________________________________________________C
C
      IF ( NNDS.GT.NNDM ) THEN
        PRINT *,' Subroutine SFSETI.'
        PRINT *,' NNDS = ', NNDS,' NNDM = ', NNDM
        PRINT *,' Program aborted.'
        STOP
      ENDIF
C
      DO ITS = 2, NTS
        RAD = DSPCA( 1, ITS - 1 )
        THE = DSPCA( 2, ITS - 1 )
        PHI = DSPCA( 3, ITS - 1 )
C       .
        ICOMP = 1
        VRAD = RTPFCE( ICOMP, INARR, NNDS, IWORK, MHT, MHL, MHM, RAD,
     1                 THE, PHI, VEC, XARR, WORK1, WORK2, COEFM )
C       .
        ICOMP = 2
        VTHE = RTPFCE( ICOMP, INARR, NNDS, IWORK, MHT, MHL, MHM, RAD,
     1                 THE, PHI, VEC, XARR, WORK1, WORK2, COEFM )
C       .
        ICOMP = 3
        VPHI = RTPFCE( ICOMP, INARR, NNDS, IWORK, MHT, MHL, MHM, RAD,
     1                 THE, PHI, VEC, XARR, WORK1, WORK2, COEFM )
C       .
        DSPCA( 1, ITS ) = RAD + DELTAT*VRAD
        DSPCA( 2, ITS ) = THE + DELTAT*VTHE
        DSPCA( 3, ITS ) = PHI + DELTAT*VPHI
C       .
      ENDDO
C
      RETURN
      END
C*********************************************************************
