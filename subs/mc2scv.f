C*********************************************************************
C subroutine Multiple Component 2 Single Component Vector ************
C            -        -         - -      -         -      ************
C Steve Gibbons Wed Nov  1 11:08:11 MET 2000                         C
C____________________________________________________________________C
C                                                                    C
C   IDIR = 1:                                                        C
C Takes a single component vector VECC defined by MLC, MMC, with     C
C NHC harmonics and fills it up, from grid node ILNR to IRNR, with   C
C values from the corresponding harmonics in the vector VEC which    C
C is defined by the arrays INARR, MHT, MHL and MHM.                  C
C                                                                    C
C   IDIR = -1:                                                       C
C Fills in values of VEC from VECC.                                  C
C                                                                    C
C Note that VECC must be arranged with ind = (ih-1)*nr + ir          C
C                                                                    C
C____________________________________________________________________C
C                                                                    C
C Input variables :-                                                 C
C ===============                                                    C
C  Integer                                                           C
C  -------                                                           C
C                                                                    C
C     INARR     : Int. parameter array corresponding to vectors.     C
C                 Dim ( * ) but the following must be true.          C
C                                                                    C
C                 INARR( 1 ) = IFORMF                                C
C                 INARR( 2 ) = NR      See INDFUN for details        C
C                 INARR( 3 ) = NH                                    C
C                                                                    C
C     MHT       : MHT( ih ) contains itype for harmonic 'ih'         C
C                  (see ICOMP)                                       C
C                                                                    C
C     MHL       : MHL( ih ) contains degree, l, for harmonic 'ih'    C
C                                                                    C
C     MHM       : MHM( ih ) contains order, m, for harmonic 'ih' if  C
C                  ih has cos m phi dependency and (-m) if ih has    C
C                   sin m phi dependency.                            C
C                                                                    C
C N.B. MHT, MHL, MHM and INARR all correspond to the array, VEC.     C
C                                                                    C
C     ILNR      : Lowest radial node to transfer.                    C
C     IRNR      : Highest radial node to transfer.                   C
C                                                                    C
C     ICOMP     : Number of component (1,2,3,4 or 5) -               C
C                 corresponds to MHT( ih ) and                       C
C                                                                    C
C         ICOMP = 1 for a poloidal velocity harmonic.                C
C         ICOMP = 2 for a toroidal velocity harmonic.                C
C         ICOMP = 3 for a temperature harmonic.                      C
C         ICOMP = 4 for a poloidal magnetic field harmonic.          C
C         ICOMP = 5 for a toroidal magnetic field harmonic.          C
C                                                                    C
C     NHC       : Number of harmonic in vector VECC                  C
C                                                                    C
C     MLC       : Dim ( * ). Equivalent of MHL, for VECC.            C
C                                                                    C
C     MMC       : Dim ( * ). Equivalent of MHM, for VECC.            C
C                                                                    C
C     IDIR      : Set to 1 to put values from VEC into VECC          C
C                 Set to -1 to put values from VECC into VEC         C
C                                                                    C
C  Double Precision                                                  C
C  ----------------                                                  C
C                                                                    C
C     VEC       : Dim ( * ). (Multiple component vector).            C
C     VECC      : Dim ( * ). (Single component vector).              C
C                                                                    C
C____________________________________________________________________C
C
C*********************************************************************
      SUBROUTINE MC2SCV( INARR, MHT, MHL, MHM, ILNR, IRNR, ICOMP,
     1                   NHC, MLC, MMC, VEC, VECC, IDIR )
      IMPLICIT NONE
C____________________________________________________________________C
C Variable declarations - Parameters ................................C
      INTEGER INARR( * ), MHT( * ), MHL( * ), MHM( * ), ILNR, IRNR,
     1        ICOMP, NHC, MLC( * ), MMC( * ), IDIR
      DOUBLE PRECISION VEC( * ), VECC( * )
C____________________________________________________________________C
C Variable declarations - Working variables .........................C
      INTEGER IND, INDC, INDFUN, NH, NR, IR, IH, IHC,
     1        IBEGIN
      EXTERNAL INDFUN
C____________________________________________________________________C
C START OF PROGRAM **************************************************C
C____________________________________________________________________C
c
      NR = INARR( 2 )
      NH = INARR( 3 )
c
      IF ( IDIR.NE.1 .AND. IDIR.NE.-1 ) THEN
        PRINT *,' Subroutine MC2SCV.'
        PRINT *,' IDIR = ', IDIR,' : Invalid value.'
        PRINT *,' Program aborted.'
        STOP
      ENDIF
c
      IF ( ICOMP.LT.1 .OR. ICOMP.GT.5 ) THEN
        PRINT *,' Subroutine MC2SCV.'
        PRINT *,' ICOMP = ', ICOMP,' : Invalid value.'
        PRINT *,' Program aborted.'
        STOP
      ENDIF
c
      DO IHC = 1, NHC
        IBEGIN = (IHC-1)*NR 
        DO IH = 1, NH
          IF ( MHT( IH ).EQ.ICOMP   .AND.   MHL( IH ).EQ.MLC( IHC )
     1           .AND.  MHM( IH ).EQ.MMC( IHC )    ) THEN
            DO IR = ILNR, IRNR
              IND  = INDFUN( IR, IH, INARR )
              INDC = IBEGIN + IR
              IF ( IDIR.EQ.1 ) VECC( INDC ) = VEC( IND )
              IF ( IDIR.EQ.-1 ) VEC( IND ) = VECC( INDC )
            ENDDO
          ENDIF
        ENDDO
      ENDDO
c
      RETURN
      END
C*********************************************************************
