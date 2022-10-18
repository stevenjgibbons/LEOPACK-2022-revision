C*********************************************************************
C subroutine Insulator Magnetic Field Matrix Boundary cond. Enforce **
C            -         -        -     -      -              -       **
C Steve Gibbons Mon Oct 18 10:59:54 BST 1999                         C
C____________________________________________________________________C
C                                                                    C
C Applies an insulating boundary condition to all magnetic field     C
C harmonic terms in the matrix. It also zeros all the necessary rows C
C                                                                    C
C Before solving, the right hand side must be treated with the       C
C routine IMGRBE.                                                    C
C                                                                    C
C This routine is only valid when IMF = 1.                           C
C____________________________________________________________________C
C Input variables :-                                                 C
C ===============                                                    C
C  Integer                                                           C
C  -------                                                           C
C     N1        : Leading dimension of the matrix.                   C
C     N2        : Second dimension of the matrix.                    C
C     KL        : Number of lower diagonals in banded matrix.        C
C     KU        : Number of upper diagonals in banded matrix.        C
C     KLE       : Number of extra lower diagonals in banded          C
C                  matrix. This is to make the matrix solvable by    C
C                   LAPACK routines.                                 C
C                                                                    C
C     IMF       : Matrix format flag.                                C
C                                                                    C
C          imf = 1; Matrix is in LAPACK banded format                C
C                   ie element a_{i,j} is stored in                  C
C                   A( kle + ku + 1 + i - j , j )                    C
C                                                                    C
C    ALL other values of IMF are disqualified for this routine.      C
C                                                                    C
C     NR        : Number of radial grid nodes.                       C
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
C                 INARR( 2 ) = NRR. Number of radial grid nodes.     C
C                    (NRR must be consistent with NR above).         C
C                 INARR( 3 ) = NH. Number of harmonics in sol. vect. C
C                                                                    C
C MHTI, MHLI and MHMI all describe the harmonics of the solution     C
C vector (i.e. describe the columns of the matrix)                   C
C                                                                    C
C     MHTI      : Array length ( * ) - atleast length NH             C
C                                                                    C
C         MHTI( IH ) = 1 --> harmonic is poloidal velocity           C
C         MHTI( IH ) = 2 --> harmonic is toroidal velocity           C
C         MHTI( IH ) = 3 --> harmonic is temperature.                C
C         MHTI( IH ) = 4 --> harmonic is poloidal magnetic field     C
C         MHTI( IH ) = 5 --> harmonic is toroidal magnetic field     C
C                                                                    C
C     MHLI      : Array length ( * ) - atleast length NH             C
C                  Sph. harm. degree, l.                             C
C     MHMI      : Array length ( * ) - atleast length NH             C
C                  Sph. harm. order, m, for cos m phi dep.           C
C                 -Sph. harm. order, m, for sin m phi dep.           C
C                                                                    C
C MHTO, MHLO and MHMO all describe the harmonics of the equations    C
C (i.e. describe the rows of the matrix) - they are directly         C
C analogous to MHTI, MHLI and MHMI.                                  C
C                                                                    C
C     NBN       : Number of bounding grid nodes for central          C
C                  differences.                                      C
C                                                                    C
C     ILNC      : Left-most node which may be used to form deriv.s   C
C     IRNC      : Right-most node which may be used to form deriv.s  C
C                                                                    C
C     NDRVS     : Highest derivative stored in FDCM.                 C
C                 Must be atleast 1.                                 C
C                                                                    C
C     NFDCM     : Leading dim of FDCM. See FDCMBD.                   C
C                  (Must be atleast 2*NBN + 1 )                      C
C                                                                    C
C  Double Precision                                                  C
C  ----------------                                                  C
C                                                                    C
C     A         : Matrix. Dim ( N1, N2 )                             C
C     XARR      : Dim ( NR ). i^{th} element gives x_i.              C
C     FDCM      : Dim ( NFDCM, NR, NDRVS ). Finite diff. coeff.s     C
C                Must be formed in advance by a call to FDCMBD.      C
C                                                                    C
C  Character                                                         C
C  ---------                                                         C
C                                                                    C
C     CHBND     : Boundary flag (*). Either 'Inner' or 'Outer'       C
C                                                                    C
C____________________________________________________________________C
C
C*********************************************************************
      SUBROUTINE IMFMBE( N1, N2, KL, KU, KLE, IMF, NR, INARR, MHTI,
     1                   MHLI, MHMI, MHTO, MHLO, MHMO, NBN, ILNC,
     2                   IRNC, NDRVS, NFDCM, A, XARR, FDCM, CHBND )
      IMPLICIT NONE
C____________________________________________________________________C
C Variable declarations - Parameters ................................C
      INTEGER N1, N2, KL, KU, KLE, IMF, NR, INARR( * ), MHTI( * ),
     1        MHLI( * ), MHMI( * ), MHTO( * ), MHLO( * ), MHMO( * ),
     2        NBN, ILNC, IRNC, NDRVS, NFDCM
      DOUBLE PRECISION A( N1, N2 ), XARR( NR ),
     1                 FDCM( NFDCM, NR, NDRVS )
      CHARACTER *(*)   CHBND
C____________________________________________________________________C
C Variable declarations - Working variables .........................C
      INTEGER IFORMF, NRR, NH, IHC, IHR, IHD, INARR2( 6 ), IRR,
     1        IRAD, IFLAG, L, IRC, NRC, INDFUN
      DOUBLE PRECISION JUNK( 1 ), WORK( 2 ), FAC, ZERO
      PARAMETER ( ZERO = 0.0d0 )
      EXTERNAL MHDBCS
      LOGICAL OK
C____________________________________________________________________C
C START OF PROGRAM **************************************************C
C____________________________________________________________________C
C
      JUNK( 1 ) = 0.0d0
C     .
C     . Check value of NR
C     .
      IFORMF = INARR( 1 )
      NRR    = INARR( 2 )
      NH     = INARR( 3 )
C     .
      IF ( NRR.NE.NR ) THEN
        PRINT *,' Subroutine IMFMBE.'
        PRINT *,' NR  = ', NR
        PRINT *,' NRR = ', NRR
        PRINT *,' Program aborted.'
        STOP
      ENDIF
C     .
C     . Check IMF
C     .
      IF ( IMF.NE.1 ) THEN
        PRINT *,' Subroutine IMFMBE.'
        PRINT *,' IMF = ', IMF
        PRINT *,' IMF = 1 is the only permitted option.'
        PRINT *,' Program aborted.'
        STOP
      ENDIF
C     .
C     . Check the validity of CHBND
C     .
      IF ( CHBND(1:1).EQ.'I' .OR. CHBND(1:1).EQ.'i' ) THEN
        IFLAG = 5
        IRR   = 1
        IRAD  = 1
        IF ( ILNC.GT.1 ) THEN
          PRINT *,' Subroutine IMFMBE.'
          PRINT *,' ILNC = ', ILNC
          PRINT *,' Must be atmost 1.'
          PRINT *,' Program aborted.'
          STOP
        ENDIF
        GOTO 50
      ENDIF
C     .
      IF ( CHBND(1:1).EQ.'O' .OR. CHBND(1:1).EQ.'o' ) THEN
        IFLAG = 6
        IRR   = NR
        IRAD  = NR
        IF ( IRNC.LT.NR ) THEN
          PRINT *,' Subroutine IMFMBE.'
          PRINT *,' IRNC = ', IRNC
          PRINT *,' Must be atleast NR.'
          PRINT *,' Program aborted.'
          STOP
        ENDIF
        GOTO 50
      ENDIF
C     .
      PRINT *,' Subroutine IMFMBE. CHBND = ', CHBND
      PRINT *,' Program aborted.'
      STOP
C     .
 50   CONTINUE
      FAC = 1.0d0
C     .
      DO IHC = 1, NH
        IF ( MHTI( IHC ).NE.4 .AND. MHTI( IHC ).NE.5 ) GOTO 60
C       .
C       . Harmonic is either poloidal or toroidal magnetic f.
C       . Let's search for corresponding harmonic in rows.
C       .
        OK = .FALSE.
        L  = MHLI( IHC )
        DO IHR = 1, NH
          IF (   MHTO( IHR ).EQ.MHTI( IHC )   .AND.
     1           MHLO( IHR ).EQ.MHLI( IHC )   .AND.
     2           MHMO( IHR ).EQ.MHMI( IHC )         ) THEN
C           .
            OK = .TRUE.
C           .
C           . Zero the row of the matrix
C           .
            IRC = 1
            NRC = INDFUN( IRR, IHR, INARR )
            CALL BMRCOP( KL, KU, KLE, N2, IRC, NRC, A,
     1                   ZERO, ZERO )
C           .
            INARR2( 1 ) = IFORMF
            INARR2( 2 ) = NR
            INARR2( 3 ) = NH
            INARR2( 6 ) = L
C           .
C           . Case of poloidal magnetic field
C           .
            IF ( MHTI( IHC ).EQ.4 ) THEN
              INARR2( 5 ) = IFLAG
              IHD         = 1
            ENDIF
C           .
C           . Case of toroidal magnetic field
C           .
            IF ( MHTI( IHC ).EQ.5 ) THEN
              INARR2( 5 ) = 1
              IHD         = 0
            ENDIF
C           .
            CALL NGMBCE( N1, N2, KL, KU, KLE, IMF, IHC, IHR, INARR2,
     1                   IHD, NBN, IRR, IRAD, ILNC, IRNC, NFDCM,
     2                   NDRVS, MHDBCS, A, FAC, XARR, FDCM,
     3                   JUNK, JUNK, NR, WORK )
C           .
          ENDIF
        ENDDO
        IF ( .NOT. OK ) THEN
          PRINT *,' Subroutine IMFMBE.'
          PRINT *,' Failiure to find equation corresponding'
          PRINT *,' to harmonic with type = ', MHTI( IHC )
          PRINT *,' l= ',MHLI( IHC ),' and m= ',MHMI( IHC )
          PRINT *,' Program aborted.'
          STOP
        ENDIF
C       .
 60   CONTINUE
      ENDDO
C     .
      RETURN
      END
C*********************************************************************

