C*********************************************************************
C subroutine POTential Field Evaluation Routine **********************
C            ---       -     -          -       **********************
C Steve Gibbons Fri Mar 16 16:17:02 MET 2001                         C
C____________________________________________________________________C
C                                                                    C
C For a given outer radius, RO, the spherical polar coordinates      C
C DRAD, DTHE, DPHI, and a magnetic field described by VEC, MHT, MHL, C
C MHM and INARR, POTFER will calculate BRAD, BTHE, BPHI.             C
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
C     MHM       : MHM( ih ) contains order, m, for harmonic 'ih' if  C
C                  ih has cos m phi dependency and (-m) if ih has    C
C                   sin m phi dependency.                            C
C                                                                    C
C  Double Precision                                                  C
C  ----------------                                                  C
C                                                                    C
C     VEC       : Solution vector. Dim ( * ) ( input )               C
C                 Length must be atleast NR*NH                       C
C                                                                    C
C     RO        : Radius of the outer boundary.                      C
C                 (In the scaling applied by dynamo code).           C
C                                                                    C
C     DRAD      : Radius where potential field is to be measured.    C
C                 (In kilometers from the Earth's centre)            C
C     DTHE      : Theta where potential field is to be measured.     C
C     DPHI      : Phi where potential field is to be measured.       C
C                                                                    C
C     BRAD      : Radial component of magnetic field.                C
C     BTHE      : Theta component of magnetic field.                 C
C     BPHI      : Phi component of magnetic field.                   C
C____________________________________________________________________C
C
C*********************************************************************
      SUBROUTINE POTFER( INARR, MHT, MHL, MHM, VEC, RO, DRAD, DTHE,
     1                   DPHI, BRAD, BTHE, BPHI )
      IMPLICIT NONE
C____________________________________________________________________C
C Variable declarations - Parameters ................................C
      INTEGER          INARR( * ), MHT( * ), MHL( * ), MHM( * )
      DOUBLE PRECISION VEC( * ), RO, DRAD, DTHE,
     1                   DPHI, BRAD, BTHE, BPHI
C____________________________________________________________________C
C Variable declarations - Working variables .........................C
      INTEGER          IND, INDFUN, IH, NR, NH, L, M
      DOUBLE PRECISION DPHTRM, PHTRM, DMPHI, P, DP, COSTH,
     1                 SHMPLG, SHDPLG, DPOL, SINTH, DLOW,
     2                 DLP2, DMLFAC, RADCMB, FAC
      EXTERNAL         INDFUN, SHMPLG, SHDPLG
      PARAMETER ( DLOW = 1.0d-9 )
C____________________________________________________________________C
C START OF PROGRAM **************************************************C
C____________________________________________________________________C
C
      BRAD = 0.0d0
      BTHE = 0.0d0
      BPHI = 0.0d0
C     .
      NR    = INARR( 2 )
      NH    = INARR( 3 )
      COSTH = DCOS( DTHE )
      SINTH = DSIN( DTHE )
C     .
C     RADEAR = 6371.0d0
      RADCMB = 3485.0d0
      FAC    = RADCMB/DRAD
C     .
      IF ( DABS( SINTH ).LT.DLOW ) THEN
        PRINT *,' Subroutine POTFER.'
        PRINT *,' SINTH = ', SINTH
        PRINT *,' Program aborted.'
        STOP
      ENDIF
C     .
      DO IH = 1, NH
C       .
        IF ( MHT( IH ).NE.4 ) GOTO 213
        IND    = INDFUN( NR, IH, INARR )
        L      = MHL( IH )
        DLP2   = DBLE( L + 2 )
        DPOL   = VEC( IND )*RO
        DMLFAC = DPOL*DBLE( L )*(FAC**DLP2)
C       .
        IF ( MHM( IH ).LT.0 ) THEN
          M      = -MHM( IH )
          DMPHI  = DBLE( M )*DPHI
          PHTRM  = DSIN( DMPHI )
          DPHTRM = DBLE( M )*DCOS( DMPHI )
        ELSE
          M      = MHM( IH )
          DMPHI  = DBLE( M )*DPHI
          PHTRM  = DCOS( DMPHI )
          DPHTRM = DBLE( M )*DSIN( DMPHI )*(-1.0d0)
        ENDIF
C       .
        P  = SHMPLG( L, M, COSTH )
        DP = SHDPLG( L, M, COSTH )
C       .
        BRAD = BRAD + DMLFAC*DBLE( L + 1 )*P*PHTRM
        BTHE = BTHE - DMLFAC*DP*PHTRM
        BPHI = BPHI - DMLFAC*P*DPHTRM/SINTH
C       .
 213  CONTINUE
      ENDDO
C     .
      RETURN
      END
C*********************************************************************
