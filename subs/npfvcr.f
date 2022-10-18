C*********************************************************************
C subroutine NPlot Format Vector Conversion Routine ******************
C            --    -      -      -          -       ******************
C Steve Gibbons Wed Dec  8 09:16:37 GMT 1999                         C
C____________________________________________________________________C
C                                                                    C
C Writes out the solution vector in the same format as nplot.        C
C This involves reordering the vector and interpolating it onto      C
C an evenly spaced grid.                                             C
C                                                                    C
C____________________________________________________________________C
C                                                                    C
C Input variables :-                                                 C
C ===============                                                    C
C  Integer                                                           C
C  -------                                                           C
C                                                                    C
C     LU        : Logical unit number of file.                       C
C     IWR       : Write flag. = 2 for caution, =3 for overwrite.     C
C     INARR     : Array for format of VEC. (See INDFUN)              C
C                                                                    C
C                 INARR( 1 ) = IFORMF                                C
C                 INARR( 2 ) = NR     See INDFUN for details         C
C                 INARR( 3 ) = NH                                    C
C                                                                    C
C     NROF      : Number of radial nodes required for the output     C
C                  file.                                             C
C                                                                    C
C     NNDS      : Number of nodes allowed for interpolation of func. C
C     IWORK     : Work array. Dim ( NNDS ).                          C
C                                                                    C
C     NOHMAX    : Number of harmonics allowed to be output.          C
C                                                                    C
C     INARRO    : Array for format of VWORK.                         C
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
C  Double Precision                                                  C
C  ----------------                                                  C
C                                                                    C
C     VEC       : Dim ( * ) but atleast NR*NH. Original vector arr.  C
C     XARR      : Dim ( * ) but atleast NR. Radial node values.      C
C                                                                    C
C     VWORK     : Dim ( * ) but atleast NROF*NOHMAX. Work array.     C
C     XWORK     : Dim ( * ) but atleast NROF. Work array.            C
C                                                                    C
C     WORK1     : Dim ( NNDS ). Work array.                          C
C     WORK2     : Dim ( NNDS ). Work array.                          C
C     WORKM     : Dim ( NNDS, NNDS ). Work array.                    C
C                                                                    C
C  Character                                                         C
C  ---------                                                         C
C                                                                    C
C     FNAME     : Name of output file.                               C
C                                                                    C
C     CHVMFF    : Velocity/temperature/magnetic field select.        C
C                 Set to either 'VEL', 'TEM' or 'MAG'                C
C                                                                    C
C                                                                    C
C____________________________________________________________________C
C
C*********************************************************************
      SUBROUTINE NPFVCR( LU, IWR, INARR, NROF, NNDS, IWORK, NOHMAX,
     1                   INARRO, MHT, MHL, MHM, VEC, XARR, VWORK,
     2                   XWORK, WORK1, WORK2, WORKM, FNAME, CHVMFF )
      IMPLICIT NONE
C____________________________________________________________________C
C Variable declarations - Parameters ................................C
      INTEGER LU, IWR, INARR( * ), NROF, NNDS, IWORK( NNDS ), NOHMAX,
     1        INARRO( * ), MHT( * ), MHL( * ), MHM( * )
      DOUBLE PRECISION VEC( * ), XARR( * ), VWORK( * ), XWORK( * ),
     1                 WORK1( NNDS ), WORK2( NNDS ),
     2                 WORKM( NNDS, NNDS )
      CHARACTER *(3) CHVMFF
      CHARACTER *(*) FNAME
C____________________________________________________________________C
C Variable declarations - Working variables .........................C
C
      INTEGER NHT1, NHT2, NHT3, NHT4, NHT5, IFLAG, NOH, IH, NH,
     1        NR, LH, IR, I, ILEN, IIH, IND, INDFUN, IH2, IPOL,
     2        ITOR, NOH2
      DOUBLE PRECISION RI, RO, LOW, ASP, RAD
      PARAMETER ( LOW = 1.0d-6 )
C____________________________________________________________________C
C START OF PROGRAM **************************************************C
C____________________________________________________________________C
C
      IF ( IWR.NE.2 .AND. IWR.NE.3 ) THEN
        PRINT *,' Subroutine NPFVCR.'
        PRINT *,' IWR = ', IWR
        PRINT *,' Program aborted.'
        STOP
      ENDIF
C     .
C     . Let's check validity of CHVMFF
C     .
      IF ( CHVMFF.EQ.'VEL' .OR. CHVMFF.EQ.'vel' .OR.
     1     CHVMFF.EQ.'Vel' ) THEN
        IPOL  = 1
        ITOR  = 2
        IFLAG = 1
      ENDIF
C     .
      IF ( CHVMFF.EQ.'MAG' .OR. CHVMFF.EQ.'mag' .OR.
     1     CHVMFF.EQ.'Mag' ) THEN
        IPOL  = 4
        ITOR  = 5
        IFLAG = 4
      ENDIF
C     .
      IF ( CHVMFF.EQ.'TEM' .OR. CHVMFF.EQ.'tem' .OR.
     1     CHVMFF.EQ.'Tem' ) THEN
        IFLAG = 3
      ENDIF
C     .
      IF ( IFLAG.NE.1 .AND. IFLAG.NE.3 .AND. IFLAG.NE.4 ) THEN
        PRINT *,' Subroutine NPFVCR.'
        PRINT *,' CHVMFF = ', CHVMFF
        PRINT *,' Program aborted.'
        STOP
      ENDIF
C     .
C     . Count up the different kinds of harmonics
C     .
      NH   = INARR( 3 )
      NHT1 = 0
      NHT2 = 0
      NHT3 = 0
      NHT4 = 0
      NHT5 = 0
      DO IH = 1, NH
        IF ( MHT( IH ).EQ.1 ) NHT1 = NHT1 + 1
        IF ( MHT( IH ).EQ.2 ) NHT2 = NHT2 + 1
        IF ( MHT( IH ).EQ.3 ) NHT3 = NHT3 + 1
        IF ( MHT( IH ).EQ.4 ) NHT4 = NHT4 + 1
        IF ( MHT( IH ).EQ.5 ) NHT5 = NHT5 + 1
      ENDDO
C  
      IF ( IFLAG.EQ.1 .AND. NHT1.EQ.0 .AND. NHT2.EQ.0 ) RETURN
      IF ( IFLAG.EQ.3 .AND. NHT3.EQ.0 ) RETURN
      IF ( IFLAG.EQ.4 .AND. NHT4.EQ.0 .AND. NHT5.EQ.0 ) RETURN
C
      IF ( IFLAG.EQ.1 ) NOH = NHT1 + NHT2
C (temperature is a special case for reasons of nplot weirdness)
      IF ( IFLAG.EQ.3 ) NOH = NHT1 + 1
      IF ( IFLAG.EQ.4 ) NOH = NHT4 + NHT5
C
      IF ( NOH.GT.NOHMAX ) THEN
        PRINT *,' Subroutine NPFVCR.'
        PRINT *,' NOH = ', NOH,' NOHMAX = ', NOHMAX
        PRINT *,' Program aborted.'
        STOP
      ENDIF
C
C Set XWORK to NROF evenly spaced nodes
C
      NR = INARR( 2 )
      RI = XARR( 1 )
      RO = XARR( NR )
C
      IF ( RO.LT.LOW ) THEN
        PRINT *,' Subroutine NPFVCR.'
        PRINT *,' RO = ', RO
        PRINT *,' Program aborted.'
        STOP
      ENDIF
      ASP = RI/RO
C
      CALL ESNAAS( NROF, XWORK, RI, RO)
C
C Set parameters in INARRO
C
      INARRO( 1 ) = 4
      INARRO( 2 ) = NROF
      INARRO( 3 ) = NOH
      ILEN = NOH*NROF
      LH   = 0
      NOH2 = 0
C
C Fill elements of VWORK
C First - the case of IFLAG.EQ.1 or IFLAG.EQ.4
C
      IF ( IFLAG.EQ.1 .OR. IFLAG.EQ.4 ) THEN
        DO IH = 1, NH
          IF ( MHT( IH ).NE.IPOL .AND. MHT( IH ).NE.ITOR ) GOTO 60
          NOH2 = NOH2 + 1
          DO IR = 1, NROF
            RAD = XWORK( IR )
C
C interpolate the radial function in VEC to find the correct value
C
            CALL SVRINT( RAD, VEC, XARR, INARR, IH, NNDS, WORK1,
     1                   IWORK, WORK2, WORKM )
C
C calculate location in VWORK
C
            IND = INDFUN( IR, NOH2, INARRO )
C
C nplot radial functions differ by a factor of R
C
            VWORK( IND ) = WORK1( 1 )*RAD
C
          ENDDO
 60     CONTINUE
        ENDDO
      ENDIF
C
C Now case of IFLAG.EQ.3 - a bit trickier due to 
C irregularities in nplot program
C
      IF ( IFLAG.EQ.3 ) THEN
C       .
C       . First we need to add the monopole term
C       .
        NOH2 = 1
        IIH  = 0
        DO IH = 1, NH
          IF ( MHT( IH ).EQ.3 .AND. MHL( IH ).EQ.0 ) IIH = IH
        ENDDO
        IF ( IIH.NE.0 ) THEN
C         .
C         . We have a monopole term in the solution
C         .
          DO IR = 1, NROF
            IND = INDFUN( IR, NOH2, INARRO )
            RAD = XWORK( IR )
            CALL SVRINT( RAD, VEC, XARR, INARR, IIH, NNDS, WORK1,
     1                   IWORK, WORK2, WORKM )
            VWORK( IND ) = WORK1( 1 )*RAD
          ENDDO
C         .
        ELSE
C         .
C         . We do not have a monopole term in the solution
C         . Instead, fill the function with zeros
C         .
          DO IR = 1, NROF
            IND = INDFUN( IR, NOH2, INARRO )
            VWORK( IND ) = 0.0d0
          ENDDO
        ENDIF
C       .
C       . Now do the terms which are not monopole terms
C       .
        DO IH = 1, NH
          IF ( MHT( IH ).EQ.1 ) THEN
C           .
C           . We have found a poloidal harmonic
C           . Now let's loop around the harmonics
C           . to see if we can find the corresponding 
C           . temperature harmonic
C           .
            IIH = 0
            DO IH2 = 1, NH
              IF ( MHT( IH2 ).EQ.3             .AND.
     1             MHL( IH2 ).EQ.MHL( IH )     .AND.
     2             MHM( IH2 ).EQ.MHM( IH )       ) IIH = IH2
            ENDDO
            IF ( IIH.NE.0 ) THEN
C             .
C             . We have an equivalent term in the solution
C             .
              NOH2 = NOH2 + 1
              DO IR = 1, NROF
                IND = INDFUN( IR, NOH2, INARRO )
                RAD = XWORK( IR )
                CALL SVRINT( RAD, VEC, XARR, INARR, IIH, NNDS, WORK1,
     1                       IWORK, WORK2, WORKM )
                VWORK( IND ) = WORK1( 1 )*RAD
              ENDDO
C             .
            ELSE
C             .
C             . We do not have an equivalent term in the solution
C             . Instead, fill the function with zeros
C             .
              NOH2 = NOH2 + 1
              DO IR = 1, NROF
                IND = INDFUN( IR, NOH2, INARRO )
                VWORK( IND ) = 0.0d0
              ENDDO
            ENDIF
C           .
          ENDIF
        ENDDO
C       .
      ENDIF
C
      IF ( NOH2.NE.NOH ) THEN
        PRINT *,' Subroutine NPFVCR.'
        PRINT *,' In principle there should be ',NOH
        PRINT *,' harmonics output. I have counted ',NOH2
        PRINT *,' Program aborted.'
        STOP
      ENDIF
C
C Open output file
C
      CALL FOPEN ( LU, FNAME, IWR )
C
      WRITE ( LU, 100 )
      WRITE ( LU, * ) '  0.0  0.0  0.0  0.0  0.0 '
      WRITE ( LU, 101 ) LH, NOH, NROF, ASP
      WRITE ( LU, 102 ) ( VWORK( I ), I = 1, ILEN )
 100  FORMAT ( ' ' )
 101  FORMAT ( i6,i6,i6,' 0.0  0.0  ',f10.4,' 0.0  0.0 ' )
 102  FORMAT ( 5e15.7 )
C
C Close file.
C
      CALL FCLOSE ( LU, FNAME, 'Error in closing file.')
C
      RETURN
      END
C*********************************************************************
