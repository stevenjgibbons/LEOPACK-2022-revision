C*********************************************************************
C subroutine Adapted Matrix Surplus Diagonal Element Addition ********
C            -       -      -       -        -       -        ********
C Steve Gibbons Wed Oct 27 09:12:14 GMT 1999                         C
C____________________________________________________________________C
C                                                                    C
C By using finite difference schemes where homogeneous boundary      C
C conditions are implicit in the coefficients used at points close   C
C to the boundary, either the end point or the two end points at     C
C each boundary become obsolete and must be arbitrarily filled       C
C to prevent singularity of the matrix.                              C
C                                                                    C
C AMSDEA loops around each harmonic, checks to see which scheme is   C
C used for this harmonic, looks to see what condition is implied     C
C at the boundary and then adds 0, 1 or 2 diagonal values of         C
C DIAGEL. If we are merely solving a linear system, DIAGEL is        C
C arbitrary. If we solving an eigensystem, DIAGEL will become an     C
C eigenvalue and so should be as highly negative as is necessary     C
C to prevent confusion with the 'real' eigenvalues.                  C
C                                                                    C
C____________________________________________________________________C
C                                                                    C
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
C                 INARR( 2 ) = NR. Number of radial grid nodes.      C
C                 INARR( 3 ) = NH. Number of harmonics in sol. vect. C
C                                                                    C
C     MHP       : Array length ( * ) - atleast length NH             C
C                  Pointer array to finite difference coefficients.  C
C                  MHPI( ih ) = is, which is the 4th index of        C
C                  array SVFDC - indicates f.d. scheme used.         C
C                                                                    C
C     IBCARR    : Dimension ( NDCS ). This must be the array MHIBC   C
C                 submitted to SVFDCF in the case CHBND(1:1) = 'I'   C
C                 and the array MHOBC submitted to SVFDCF in the     C
C                 case CHBND(1:1) = 'O'                              C
C                                                                    C
C     IBCARR( IS ) contains the integer IBC which corresponds        C
C    to the following boundary conditions :-                         C
C                                                                    C
C   IBC       Boundary condition                   No. of nodes      C
C   ---       ------------------                   ------------      C
C                                                                    C
C    1           None                                    0           C
C    2           f = 0                                   1           C
C    3           df/dr = 0                               1           C
C    4           f = 0   and     df/dr = 0               2           C
C    5           f = 0   and   d^2f/dr^2 = 0             2           C
C    6           f/r - df/dr = 0                         1           C
C    7           Insulating magnetic field               1           C
C                                                                    C
C                                                                    C
C                                                                    C
C     NDCS       : Number of distinct differencing coeff.s           C
C                  represented in SVFDC.                             C
C                                                                    C
C  Character                                                         C
C  ---------                                                         C
C                                                                    C
C     CHBND     : Boundary flag (*). Either 'Inner' or 'Outer'       C
C                                                                    C
C  Double Precision                                                  C
C  ----------------                                                  C
C                                                                    C
C     A         : Matrix. Dim ( N1, N2 )                             C
C     DIAGEL    : D.p. const. to be added to diagonals.              C
C                                                                    C
C____________________________________________________________________C
C
C*********************************************************************
      SUBROUTINE AMSDEA( A, N1, N2, KL, KU, KLE, IMF, INARR, 
     1                   MHP, IBCARR, CHBND, DIAGEL, NDCS )
      IMPLICIT NONE
C____________________________________________________________________C
C Variable declarations - Parameters ................................C
      INTEGER N1, N2, KL, KU, KLE, IMF, INARR( * ), MHP( * ), NDCS,
     1        IBCARR( NDCS )
      DOUBLE PRECISION A( N1, N2 ), DIAGEL
      CHARACTER *(*)   CHBND
C____________________________________________________________________C
C Variable declarations - Working variables .........................C
      INTEGER NR, NH, IRAD, INDFUN, ISN, IEN, IH, IS, IOF,
     1        NRC, IRC, IBC
      DOUBLE PRECISION ZERO
      PARAMETER ( ZERO = 0.0d0 )
C____________________________________________________________________C
C START OF PROGRAM **************************************************C
C____________________________________________________________________C
C
      NR = INARR( 2 )
      NH = INARR( 3 )
C     .
C     . Check on IMF (1 is only acceptable value)
C     .
      IF ( IMF.NE.1 ) THEN
        PRINT *,' Subroutine AMSDEA'
        PRINT *,' IMF = ', IMF
        PRINT *,' Program aborted.'
        STOP
      ENDIF
C     .
C     . Check on size of N2
C     .
      IF ( N2.NE.NH*NR ) THEN
        PRINT *,' Subroutine AMSDEA'
        PRINT *,' N2 = ', N2
        PRINT *,' NR = ', NR
        PRINT *,' NH = ', NH
        PRINT *,' Program aborted.'
        STOP
      ENDIF
C     .
C     . Set inner outer flag
C     .
      IF ( CHBND(1:1).EQ.'I' .OR. CHBND(1:1).EQ.'i' ) THEN
        IOF = 1
        GOTO 50
      ENDIF
C     .
      IF ( CHBND(1:1).EQ.'O' .OR. CHBND(1:1).EQ.'o' ) THEN
        IOF = 2
        GOTO 50
      ENDIF
C     .
      PRINT *,' Subroutine AMSDEA'
      PRINT *,' CHBND = ', CHBND
      PRINT *,' Program aborted.'
      STOP
C     .
 50   CONTINUE
C     .
C     . Begin loop around harmonics
C     .
      DO IH = 1, NH
        IS  = MHP( IH )
        IF ( IS.LT.1 .OR. IS.GT.NDCS ) THEN
          PRINT *,' Subroutine AMSDEA'
          PRINT *,' IS   = ', IS
          PRINT *,' NDCS = ', NDCS
          PRINT *,' Program aborted.'
          STOP
        ENDIF
        IBC = IBCARR( IS )
        IF ( IBC.EQ.1 ) GOTO 51
C ............................ inner boundary case
C
        IF ( IOF.EQ.1 ) THEN
C         .
          IF ( IBC.EQ.2 .OR. IBC.EQ.3 .OR. IBC.EQ.6
     1           .OR. IBC.EQ.7 ) THEN
             ISN = 1
             IEN = 1
             GOTO 61
          ENDIF
C         .
          IF ( IBC.EQ.4 .OR. IBC.EQ.5 ) THEN
             ISN = 1
             IEN = 2
             GOTO 61
          ENDIF
C         .
          PRINT *,' Subroutine AMSDEA'
          PRINT *,' IBCARR(',IS,') = ', IBC
          PRINT *,' Program aborted.'
          STOP
        ENDIF
C ............................ outer boundary case
C
        IF ( IOF.EQ.2 ) THEN
C         .
          IF ( IBC.EQ.2 .OR. IBC.EQ.3 .OR. IBC.EQ.6
     1           .OR. IBC.EQ.7 ) THEN
             ISN = NR
             IEN = NR
             GOTO 61
          ENDIF
C         .
          IF ( IBC.EQ.4 .OR. IBC.EQ.5 ) THEN
             ISN = NR - 1
             IEN = NR
             GOTO 61
          ENDIF
C         .
          PRINT *,' Subroutine AMSDEA'
          PRINT *,' IBCARR(',IS,') = ', IBC
          PRINT *,' Program aborted.'
          STOP
        ENDIF
C ............................ 
 61     CONTINUE
        DO IRAD = ISN, IEN
C           .
C           . Zero the row of the matrix
C           .
            IRC = 1
            NRC = INDFUN( IRAD, IH, INARR )
            CALL BMRCOP( KL, KU, KLE, N2, IRC, NRC, A,
     1                   ZERO, ZERO )
C           .
C           . Enter value of the diagonal element
C           .
            A( KLE + KU + 1, NRC ) = DIAGEL
C           .
        ENDDO
C ............................ 
 51   CONTINUE
      ENDDO
C     .
      RETURN
      END
C*********************************************************************
