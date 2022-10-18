C*********************************************************************
C subroutine Voritcity Equation Solution Routine *********************
C            -         -        -        -       *********************
C Steve Gibbons Tue Nov 16 10:01:06 GMT 1999                         C
C____________________________________________________________________C
C                                                                    C
C Solves the matrix equation A x = b for a solution vector x         C
C containing velocity and temperature to some formulation of the     C
C vorticity equation. The matrix A is stored in LAPACK band storage  C
C (IMF=1) with KU = KL = KLE and N1 = 3*KL + 1.                      C
C                                                                    C
C VESR may be used in two cases.                                     C
C                                                                    C
C The first of these is that the matrix has only just been formed    C
C and so we must perform an LU decomposition. This is specified by   C
C setting ILUDF = 1.                                                 C
C                                                                    C
C In this case, it will be assumed that rows corresponding to IR = 1 C
C and IR = NR (temperature and poloidal rows) and IR = 1, 2, NR-1    C
C and NR (toroidal rows) have been left blank.                       C
C                                                                    C
C Diagonal elements (value DIAGEL) are supplied to prevent           C
C singularity of the matrix.                                         C
C This is achieved by a call to AMSDEA.                              C
C                                                                    C
C If ILUDF = 2, then it will be assumed that the matrix has already  C
C been LU decomposed and the system merely needs to be solved        C
C one more time.                                                     C
C                                                                    C
C VESR also takes in the logical flag OSBR returned by SBRRFC        C
C which specifies whether we need to impose a condition of no solid  C
C body rotation. This is achieved using the vector SBRVEC, also      C
C returned by SBRRFC.                                                C
C                                                                    C
C ITGN is the number of the grid node whose equation must be         C
C replaced with the solid body rotation condition.                   C
C                                                                    C
C                                                                    C
C____________________________________________________________________C
C                                                                    C
C Input variables :-                                                 C
C ===============                                                    C
C  Integer                                                           C
C  -------                                                           C
C                                                                    C
C     N1        : First dimension of matrix, A.                      C
C     N2        : Second dimension of matrix, A.                     C
C     KL        : Number of lower diagonals in matrix                C
C                                                                    C
C     ILUDF     : =1 if matrix has just been formed and LU decomp.   C
C                    is required.                                    C
C                 =2 if LU decomp. has already been done and we      C
C                    simply require a further solution.              C
C                                                                    C
C     INARR     : Int. parameter array corresponding rows and cols.  C
C                 Dim ( * ) but the following must be true.          C
C                                                                    C
C                 INARR( 1 ) = IFORMF                                C
C                 INARR( 2 ) = NRR     See INDFUN for details        C
C                 INARR( 3 ) = NH                                    C
C                                                                    C
C     MHT       : Array length ( * ) - atleast length NH             C
C                 See below for key. (corresponds to cols, not rows) C
C                                                                    C
C     MHT( i ) = 1 for a poloidal velocity harmonic, i.              C
C     MHT( i ) = 2 for a toroidal velocity harmonic, i.              C
C     MHT( i ) = 3 for a temperature harmonic, i.                    C
C     MHT( i ) = 4 for a poloidal magnetic field harmonic, i.        C
C     MHT( i ) = 5 for a toroidal magnetic field harmonic, i.        C
C                                                                    C
C     MHL       : Array length ( * ) - atleast length NH             C
C                  Sph. harm. degree, l.                             C
C     MHM       : Array length ( * ) - atleast length NH             C
C                  Sph. harm. order, m, for cos m phi dep.           C
C                 -Sph. harm. order, m, for sin m phi dep.           C
C                                                                    C
C     MHP       : Array length ( * ) - atleast length NH             C
C                  Pointer array to finite difference coefficients.  C
C                  MHP( ih ) = is, which is the 4th index of         C
C                  array SVFDC - indicates f.d. scheme used.         C
C                                                                    C
C     IPIV      : Work array dim ( N2 ) (LAPACK pivotting info)      C
C                                                                    C
C     ITGN      : Number of grid node whose equation is to be        C
C                 replaced by the zero net ang. mom. condition       C
C                 on the T_1^0 radial function.                      C
C                 Not referenced if OSBR is .FALSE.                  C
C                                                                    C
C     MHIBC     : Inner boundary condition. Dim ( NDCS )             C
C                                                                    C
C  MHIBC( ih ) = 1 --> No condition is to be assumed at the bndry.   C
C  MHIBC( ih ) = 2 --> Function must vanish at the bndry.            C
C  MHIBC( ih ) = 3 --> First derivative must vanish at the bndry.    C
C  MHIBC( ih ) = 4 --> Function must vanish at the bndry AND         C
C                       first derivative must vanish at the bndry.   C
C  MHIBC( ih ) = 5 --> Function must vanish at the bndry AND         C
C                       second derivative must vanish at the bndry.  C
C  MHIBC( ih ) = 6 --> r df/dr - f(r) = 0   at the bndry.            C
C  MHIBC( ih ) = 7 --> r df/dr - l f(r) = 0 at the bndry.            C
C                        where L = LARR( ih )                        C
C                                                                    C
C     MHOBC     : Outer boundary condition. Dim ( NDCS )             C
C                                                                    C
C  MHOBC( ih ) = 1 --> No condition is to be assumed at the bndry.   C
C  MHOBC( ih ) = 2 --> Function must vanish at the bndry.            C
C  MHOBC( ih ) = 3 --> First derivative must vanish at the bndry.    C
C  MHOBC( ih ) = 4 --> Function must vanish at the bndry AND         C
C                       first derivative must vanish at the bndry.   C
C  MHOBC( ih ) = 5 --> Function must vanish at the bndry AND         C
C                       second derivative must vanish at the bndry.  C
C  MHOBC( ih ) = 6 --> r df/dr - f(r) = 0   at the bndry.            C
C  MHOBC( ih ) = 7 --> r df/dr + (l+1) f(r) = 0 at the bndry.        C
C                        where L = LARR( ih )                        C
C                                                                    C
C     NDCS      : Number of distinct finite difference schemes.      C
C                                                                    C
C  Double Precision                                                  C
C  ----------------                                                  C
C                                                                    C
C     A         : Matrix. Dim ( N1, N2 ).                            C
C     SBRVEC    : Additional row of matrix which must be set to zero C
C                 to enforce no solid-body rotation.                 C
C                 Returned by SBRRFC. Dim ( N2 )                     C
C                 Note that the element corresponding to ITGN will   C
C                 be altered after the first call ( ILUDF.EQ.1 ).    C
C                                                                    C
C     RHS       : Right hand side of equation. Dim ( N2 )            C
C                                                                    C
C     WVEC      : Work vector. Dim ( N2 )                            C
C                                                                    C
C     DIAGEL    : Arbitrary constant to be put in the diagonal       C
C                 element of redundant nodes. Must be non-zero to    C
C                 prevent singularity of matrix and must not be      C
C                 close to a part of the eigenspectra which is of    C
C                 interest.                                          C
C                                                                    C
C  Logical                                                           C
C  -------                                                           C
C                                                                    C
C     OSBR      : .TRUE. if and only if we must enforce the no       C
C                 solid body rotation condition.                     C
C                                                                    C
C____________________________________________________________________C
C
C*********************************************************************
      SUBROUTINE VESR( N1, N2, KL, ILUDF, INARR, MHT, MHL, MHM, MHP,
     1                 IPIV, ITGN, MHIBC, MHOBC, NDCS, A, SBRVEC,
     2                 RHS, WVEC, DIAGEL, OSBR )
      IMPLICIT NONE
C____________________________________________________________________C
C Variable declarations - Parameters ................................C
      INTEGER N1, N2, KL, ILUDF, INARR( 3 ), MHT( * ), MHL( * ),
     1        MHM( * ), MHP( * ), IPIV( * ), ITGN, MHIBC( * ),
     2        MHOBC( * ), NDCS
      DOUBLE PRECISION A( N1, N2 ), SBRVEC( * ), RHS( * ), WVEC( * ),
     1                 DIAGEL
      LOGICAL OSBR
C____________________________________________________________________C
C Variable declarations - Working variables .........................C
      INTEGER INDFUN, IMF, NR, NH, ILEN, INDSBR, IH, IRC,
     1        INFO, NRHS, I
      DOUBLE PRECISION ZERO, LOW, VY, VZ, FAC
      CHARACTER *(1) TRANS
      PARAMETER ( IMF = 1, ZERO = 0.0d0, LOW = 1.0d-9, NRHS = 1,
     1            TRANS = 'N' )
C____________________________________________________________________C
C START OF PROGRAM **************************************************C
C____________________________________________________________________C
C
      NR   = INARR( 2 )
      NH   = INARR( 3 )
      ILEN = NR*NH
C     .
      IF ( N2.NE.ILEN ) THEN
        PRINT *,' Subroutine VESR.'
        PRINT *,' N2 = ', N2
        PRINT *,' NR = ', NR
        PRINT *,' NH = ', NH
        PRINT *,' Program aborted.'
        STOP
      ENDIF
C     .
      IF ( N1.NE.(3*KL+1) ) THEN
        PRINT *,' Subroutine VESR.'
        PRINT *,' N1 = ', N1
        PRINT *,' KL = ', KL
        PRINT *,' Program aborted.'
        STOP
      ENDIF
C     .
      IF ( ILUDF.NE.1 .AND. ILUDF.NE.2 ) THEN
        PRINT *,' Subroutine VESR.'
        PRINT *,' ILUDF = ', ILUDF
        PRINT *,' Program aborted.'
        STOP
      ENDIF
C     .
      IF ( ABS( DIAGEL ).LT.LOW ) THEN
        PRINT *,' Subroutine VESR.'
        PRINT *,' DIAGEL = ', DIAGEL
        PRINT *,' Program aborted.'
        STOP
      ENDIF
C     .
      IF ( OSBR .AND. (ITGN.LT.2 .OR. ITGN.GE.NR) ) THEN
        PRINT *,' Subroutine VESR.'
        PRINT *,' ITGN = ', ITGN
        PRINT *,' Program aborted.'
        STOP
      ENDIF
C     .
C     . Find out which row of the matrix needs to
C     . be replaced to solve the solid body rotation.
C     .
      IF ( OSBR ) THEN
        DO IH = 1, NH
          IF ( MHT( IH ).EQ.2 .AND. MHL( IH ).EQ.1 .AND.
     1         MHM( IH ).EQ.0                  ) THEN
             INDSBR = INDFUN( ITGN, IH, INARR )
          ENDIF
        ENDDO
      ENDIF
C     .
C     . See if we need to do an LU decomposition.
C     . If not, then jump to line 50
C     .
      IF ( ILUDF.EQ.2 ) GOTO 50
C     .
C     . OK - we need to do an LU decomposition.
C     . First zero the row of the matrix corresponding
C     . to the solid body rotation - if required.
C     .
      IF ( OSBR ) THEN
        IRC = 1
        CALL BMRCOP( KL, KL, KL, N2, IRC, INDSBR, A, ZERO, ZERO )
C       .
C       . A is now singular so put in 1.0d0 in the
C       . appropriate diagonal element.
C       .
        A( 2*KL + 1, INDSBR ) = 1.0d0
C       .
C       . Subtract the same amount from solid body
C       . rotation vector.
C       .
        SBRVEC( INDSBR ) = SBRVEC( INDSBR ) - 1.0d0
C       .
      ENDIF
C     .
C     . OK - now we must make sure the matrix isn't
C     . singular because of redundant radial grid nodes
C     .
      CALL AMSDEA( A, N1, N2, KL, KL, KL, IMF, INARR, MHP,
     1             MHIBC, 'Inner boundary', DIAGEL, NDCS )
C     .
      CALL AMSDEA( A, N1, N2, KL, KL, KL, IMF, INARR, MHP,
     1             MHOBC, 'Outer boundary', DIAGEL, NDCS )
C     .
C     . Finally, proceed with the LU decomposition
C     .
      CALL DGBTRF( N2, N2, KL, KL, A, N1, IPIV, INFO )
C     .
C     . Check for an error from LU decomp.
C     .
      IF ( INFO.NE.0 ) THEN
        PRINT *,' Subroutine VESR.'
        PRINT *,' The LAPACK subroutine DGBTRF has'
        PRINT *,' returned ',INFO,' as a value of '
        PRINT *,' INFO in LU decomposition of matrix.'
        PRINT *,' Program aborted.'
        STOP
      ENDIF
C     .
 50   CONTINUE
C     .
C     . OK - LU decomposition is done -
C     . can now solve the system
C     .
C     . Zero row INDSBR of right hand side vector if
C     . necessary
C     .
      IF ( OSBR ) RHS( INDSBR ) = 0.0d0
C     .
C     . Now solve the equation A. x = RHS
C     .
      CALL DGBTRS( TRANS, N2, KL, KL, NRHS, A, N1,
     1             IPIV, RHS, N2, INFO )
C     .
C     . Check the solution has occured without error ...
C     .
      IF ( INFO.NE.0 ) THEN
        PRINT *,' Subroutine VESR.'
        PRINT *,' The LAPACK subroutine DGBTRS'
        PRINT *,' has returned ',INFO,' as a value of'
        PRINT *,' INFO solving A.y = RHS .'
        PRINT *,' Program aborted.'
        STOP
      ENDIF
C     .
C     . If solution was successful and we do not need
C     . to use the Sherman-Morrison formula then RHS
C     . already contains our solution and we may exit
C     .
      IF ( .NOT. OSBR ) RETURN
C     .
C     . OK - construct unit vector in the work array
C     . WVEC. This must have a 1 in the row of INDSBR
C     .
      DO I = 1, N2
        WVEC( I ) = ZERO
      ENDDO
      WVEC( INDSBR ) = 1.0d0
C     .
C     . Now solve A.z = WVEC equation (c.f. Press et al.)
C     .
      CALL DGBTRS( TRANS, N2, KL, KL, NRHS, A, N1,
     1             IPIV, WVEC, N2, INFO )
C     .
C     . Check the solution has occured without error ...
C     .
      IF ( INFO.NE.0 ) THEN
        PRINT *,' Subroutine VESR.'
        PRINT *,' The LAPACK subroutine DGBTRS'
        PRINT *,' has returned ',INFO,' as a value of'
        PRINT *,' INFO solving A.z = WVEC .'
        PRINT *,' Program aborted.'
        STOP
      ENDIF
C     .
C     . Now calculate VY ( the scalar product of
C     . vectors RHS and SBRVEC
C
      VY = 0.0d0
      DO I = 1, N2
        VY = VY + SBRVEC( I )*RHS( I )
      ENDDO
C     .
C     . Now calculate VZ ( the scalar product of
C     . vectors WVEC and SBRVEC
C     .
      VZ = 0.0d0
      DO I = 1, N2
        VZ = VZ + SBRVEC( I )*WVEC( I )
      ENDDO
C     .
C     . Now calculate FAC = [ VY / ( 1 + VZ ) ]
C     .
      FAC = VY/(1.0d0+VZ)
C     .
C     . Finally, calculate X from the Sherman Morrison formula
C     .
      DO I = 1, N2
        RHS( I ) = RHS( I ) - FAC*WVEC( I )
      ENDDO
C     .
      RETURN
      END
C*********************************************************************
