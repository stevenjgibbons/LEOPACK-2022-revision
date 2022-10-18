C*********************************************************************
C subroutine Magnetic Field Generalised Eigenvalue Problem Solve *****
C            -        -     -           -          -       -     *****
C Steve Gibbons Mon Jan 31 10:48:13 GMT 2000                         C
C____________________________________________________________________C
C                                                                    C
C Solves the generalised eigenvalue problem for critical stability   C
C of a time-independent magnetic field.                              C
C                                                                    C
C curl ( u x B ) = sigma CMD ( -\nabla^2 ) B                         C
C                                                                    C
C where CMD is passed in as a double precision scalar (normally =1)  C
C                                                                    C
C____________________________________________________________________C
C                                                                    C
C  Integer                                                           C
C  -------                                                           C
C                                                                    C
C     NR        : Number of radial grid nodes.                       C
C     INARR     : Int. parameter array corresponding to VI.          C
C                 Dim ( * ) but the following must be true.          C
C                                                                    C
C                 INARR( 1 ) = IFORMF                                C
C                 INARR( 2 ) = NRR      See INDFUN for details       C
C                 INARR( 3 ) = NH                                    C
C                                                                    C
C     MHT      : Array length ( * ) - atleast length NH              C
C                                                                    C
C         MHT( IH ) = 1 --> harmonic is poloidal velocity            C
C         MHT( IH ) = 2 --> harmonic is toroidal velocity            C
C         MHT( IH ) = 3 --> harmonic is temperature.                 C
C         MHT( IH ) = 4 --> harmonic is poloidal magnetic field      C
C         MHT( IH ) = 5 --> harmonic is toroidal magnetic field      C
C                                                                    C
C  MHT( ih ) must only equal 4 and 5 for this routine ...            C
C                                                                    C
C     MHL      : Array length ( * ) - atleast length NH              C
C                  Sph. harm. degree, l.                             C
C     MHM      : Array length ( * ) - atleast length NH              C
C                  Sph. harm. order, m, for cos m phi dep.           C
C                 -Sph. harm. order, m, for sin m phi dep.           C
C                                                                    C
C     MHP      : Array length ( * ) - atleast length NH              C
C                  Pointer array to finite difference coefficients.  C
C                  MPI( ih ) = is, which is the 4th index of         C
C                  array SVFDC - indicates f.d. scheme used.         C
C                                                                    C
C     NBN       : Number of nodes on each side of point for          C
C                  central differences.                              C
C                                                                    C
C     KL        : Number of lower diagonals in matrix                C
C                                                                    C
C     NCFM      : Leading dim of SVFDC. See SVFDCF.                  C
C                  (Must be atleast 2*NBN + 1 )                      C
C                                                                    C
C     NDRVM     : Limit on NDRVS. Array bound for SVFDC.             C
C      (NDRVM must be atleast 2 and NDRVS must be 2 for atleast      C
C       grid nodes IR = 2, NR - 1 ).                                 C
C                                                                    C
C     N1        : First dimension of matrix, A.                      C
C     N2        : Second dimension of matrix, A.                     C
C                                                                    C
C     NDCS      : Maximum distinct finite difference schemes         C
C                  stored by the array SVFDC.                        C
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
C                                                                    C
C     NEV       : Number of requested eigenvalues.                   C
C                                                                    C
C     NCV       : Length of Arnoldi factorisation.                   C
C     NCVM      : Maximum length of Arnoldi factorisation.           C
C                                                                    C
C ncv must be less than or equal to N2 and greater than or equal     C
C to (nev+2).                                                        C
C                                                                    C
C     MXIT      : Maximum number of Arnoldi iterations allowed.      C
C                                                                    C
C     NCE       : Number of converged eigenvalues obtained.          C
C                                                                    C
C     IPIV      : Dim ( N2 ). Array for pivotting.                   C
C                                                                    C
C  Double Precision                                                  C
C  ----------------                                                  C
C                                                                    C
C     SVFDC     : Finite difference coefficient matrix.              C
C                  Dimension ( NCFM, NR, NDRVM+1, NDCS ).            C
C                   Array is generated by the routine svfdcf         C
C                 See documentation for SVFDCF for details.          C
C                                                                    C
C     A         : Matrix. Dim ( N1, N2 )                             C
C                                                                    C
C     XARR      : Dim ( NR ). i^{th} element gives x_i.              C
C                                                                    C
C     DIAGEL    : Arbitrary constant to be put in the diagonal       C
C                 element of redundant nodes. Must be non-zero to    C
C                 prevent singularity of matrix and must not be      C
C                 close to a part of the eigenspectra which is of    C
C                 interest. (For growth rate problems, DIAGEL will   C
C                 ideally be very negative.)                         C
C                                                                    C
C     DRSV      : Real part of the applied shift.                    C
C                                                                    C
C     DR        : Returned containing real part of eigenvalues.      C
C     DI        : Returned containing imag. part of eigenvalues.     C
C     D3        : Returned containing direct residuals.              C
C                                                                    C
C     WORKEV    : Dimension ( 3*NCVM ). Work array.                  C
C     WORKD     : Dimension ( 3*N2 ). Work array.                    C
C     WORKL     : Dimension ( 3*NCVM*NCVM + 6*NCVM ). Work array.    C
C                                                                    C
C     ARTOL     : Stopping criterion for Arnoldi iteration.          C
C                 (Make this small).                                 C
C                                                                    C
C     RESID     : Dimension ( N2 ). Work array.                      C
C                                                                    C
C     V         : Dimension ( N2, NCVM ). Returned eigenvectors.     C
C                                                                    C
C     W2        : Dimension ( N2 ). Work array.                      C
C                                                                    C
C     CMD       : Coefficient of magnetic diffusion (=1?)            C
C                                                                    C
C  Logical                                                           C
C  -------                                                           C
C                                                                    C
C     SELECT    : Dimension ( NCVM ). Work array.                    C
C                                                                    C
C____________________________________________________________________C
C
C*********************************************************************
      SUBROUTINE MFGEPS( NR, INARR, MHT, MHL, MHM, MHP, NBN, KL,
     1       NCFM, NDRVM, SVFDC, A, N1, N2, NDCS, XARR, 
     2       DIAGEL, MHIBC, MHOBC, NEV, NCV, NCVM, MXIT, NCE,
     3       DRSV, SELECT, DR, DI, D3, WORKEV, WORKD, WORKL, IPIV,
     4       ARTOL, RESID, V, W2, CMD )
      IMPLICIT NONE
C____________________________________________________________________C
C Variable declarations - Parameters ................................C
      INTEGER NR, INARR( 3 ), MHT( * ), MHL( * ), MHM( * ), MHP( * ),
     1        NBN, KL, NCFM, NDRVM, N1, N2, NDCS, MHIBC( NDCS ),
     2        MHOBC( NDCS ), NEV, NCV, NCVM, MXIT, NCE, IPIV( * )
      DOUBLE PRECISION A( N1, N2 ), XARR( NR ), DIAGEL,
     1                 DR( NCVM ), DI( NCVM ), D3( NCVM ), DRSV,
     2                 WORKEV( 3*NCVM ), WORKD( 3*N2 ), CMD,
     3                 WORKL( 3*NCVM*NCVM + 6*NCVM ), ARTOL
      DOUBLE PRECISION RESID( N2 ), V( N2, NCVM ), 
     1                 W2( N2 ), SVFDC( NCFM, NR, NDRVM+1, NDCS )
      LOGICAL SELECT( NCVM )
C____________________________________________________________________C
C Variable declarations - Working variables .........................C
C
      INTEGER NH, NRR, IFORMF, IONE, 
     1        ICMP, ILN, J, NDRVS, IRN, IMF, IOP, ILEN
      DOUBLE PRECISION DZERO, DISV, DLOW, FAC, FACTR, FACTI,
     1                 DMONE, DLAPY2, DNRM2
      PARAMETER ( IONE = 1, DZERO = 0.0d0, DMONE = -1.0d0,
     1            DISV = 0.0d0, DLOW = 1.0d-7, IMF = 1 )
C____________________________________________________________________C
C Variable declarations - Local variables for ARPACK routines  ......C
C
      INTEGER IPARAM( 11 ), IPNTR( 14 ), IDO, INFO, LWORKL,
     1        IPN1, IPN2, NRHS
      CHARACTER *(1) HOWMNY, BMAT, TRANS
      CHARACTER *(2) WHICH
      LOGICAL EVECF, FIRST
      PARAMETER ( WHICH = 'LM', HOWMNY = 'A', BMAT = 'I',
     1            TRANS = 'N', EVECF = .TRUE., NRHS = 1 )
C____________________________________________________________________C
C START OF PROGRAM **************************************************C
C____________________________________________________________________C
C
C Set necessary values of IPARAM
C
      IPARAM( 1 ) = 1
      IPARAM( 3 ) = MXIT
      IPARAM( 4 ) = 1
      IPARAM( 7 ) = 1
C
      LWORKL = 3*NCVM*NCVM + 6*NCVM
C
      IFORMF = INARR( 1 )
      NRR    = INARR( 2 )
      NH     = INARR( 3 )
      ILEN   = NR*NH
      NDRVS  = 2
C     .
      ILN    = 2
      IRN    = NR - 1
C     .
      IF ( IFORMF.NE.3 ) THEN
        PRINT *,' Subroutine MFGEPS.'
        PRINT *,' IFORMF     = ', IFORMF
        PRINT *,' Program aborted.'
        STOP
      ENDIF
C     .
      IF ( NR.NE.NRR ) THEN
        PRINT *,' Subroutine MFGEPS.'
        PRINT *,' NR         = ', NR
        PRINT *,' INARR( 2 ) = ', NRR
        PRINT *,' Program aborted.'
        STOP
      ENDIF
C     .
      IF ( N2.NE.ILEN ) THEN
        PRINT *,' Subroutine MFGEPS.'
        PRINT *,' N2 = ', N2
        PRINT *,' NR = ', NR
        PRINT *,' NH = ', NH
        PRINT *,' Program aborted.'
        STOP
      ENDIF
C     .
      IF ( N1.NE.(3*KL+1) ) THEN
        PRINT *,' Subroutine MFGEPS.'
        PRINT *,' N1 = ', N1
        PRINT *,' KL = ', KL
        PRINT *,' Program aborted.'
        STOP
      ENDIF
C     .
      IF ( ABS( DIAGEL ).LT.DLOW ) THEN
        PRINT *,' Subroutine MFGEPS.'
        PRINT *,' DIAGEL = ', DIAGEL
        PRINT *,' Program aborted.'
        STOP
      ENDIF
C     .
      IF ( CMD.LT.DLOW ) THEN
        PRINT *,' Subroutine MFGEPS.'
        PRINT *,' CMD = ', CMD
        PRINT *,' Program aborted.'
        STOP
      ENDIF
C     .
C     . Now we must evaluate ( A - DRSV * B )
C     .
      FAC = DRSV*CMD
      ICMP = 4
      CALL AMLP( NR, INARR, MHT, MHL, MHM, MHP, ICMP, MHT,
     1           MHL, MHM, FAC, NBN, NDRVS, NDRVM, NCFM,
     2           ILN, IRN, SVFDC, XARR, NDCS, A, N1, N2, IMF,
     3           KL, KL, KL )
C     .
      ICMP = 5
      CALL AMLP( NR, INARR, MHT, MHL, MHM, MHP, ICMP, MHT,
     1           MHL, MHM, FAC, NBN, NDRVS, NDRVM, NCFM,
     2           ILN, IRN, SVFDC, XARR, NDCS, A, N1, N2, IMF,
     3           KL, KL, KL )
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
        PRINT *,' Subroutine MFGEPS.'
        PRINT *,' The LAPACK subroutine DGBTRF has'
        PRINT *,' returned ',INFO,' as a value of '
        PRINT *,' INFO in LU decomposition of matrix.'
        PRINT *,' Program aborted.'
        STOP
      ENDIF
C     .
C     . OK - LU decomposition is done -
C     . Now, we begin main Arnoldi loop.
C     .
      IDO = 0
C     .
 10   CONTINUE
C     .
      CALL DNAUPD ( IDO, BMAT, N2, WHICH, NEV, ARTOL, RESID,
     1              NCV, V, N2, IPARAM, IPNTR, WORKD,
     2              WORKL, LWORKL, INFO )
C     .
C     . Check for error on exit from dnaupd
C     .
      IF ( INFO.NE.0 ) THEN
        PRINT *,' Subroutine MFGEPS.'
        PRINT *,' The ARPACK subroutine DNAUPD has'
        PRINT *,' returned ',INFO,' as a value of INFO.'
        PRINT *,' Program aborted.'
        STOP
      ENDIF
C     .
      IF ( IDO.EQ.-1 .OR. IDO.EQ.1 ) THEN
C       .
C       . Here we must calculate y such that
C       . y = INV[ A - drsv*B ]*B*x
C       .
C       . X = WORKD( ipntr(1) )
C       . Y = WORKD( ipntr(2) )
C       .
        IPN1 = IPNTR( 1 )
        IPN2 = IPNTR( 2 )
C       .
C       . Zero 'y' vector -- i.e. workd( ipn2 )
C       .
        IOP = 0
        CALL VECOP( WORKD( IPN2 ), DZERO, N2, IOP )
C       .
C       . Now add B*x to y -
C       .
        FAC  = (-1.0d0)*CMD
        ICMP = 4
        CALL ASVLP( NR, NDCS, WORKD( IPN1 ), INARR, MHT, MHL, MHM,
     1              ICMP, MHP, WORKD( IPN2 ), INARR, MHT, MHL, MHM,
     2              FAC, NBN, NDRVS, NDRVM, ILN, IRN, SVFDC, XARR,
     3              NCFM )
C       .
        ICMP = 5
        CALL ASVLP( NR, NDCS, WORKD( IPN1 ), INARR, MHT, MHL, MHM,
     1              ICMP, MHP, WORKD( IPN2 ), INARR, MHT, MHL, MHM,
     2              FAC, NBN, NDRVS, NDRVM, ILN, IRN, SVFDC, XARR,
     3              NCFM )
C       .
C       . WORKD( IPN2 ) now contains B*x so call DGBTRS to solve
C       .
        CALL DGBTRS( TRANS, N2, KL, KL, NRHS, A, N1,
     1               IPIV, WORKD( IPN2 ), N2, INFO )
C       .
        IF ( INFO.NE.0 ) THEN
          PRINT *,' Subroutine MFGEPS.'
          PRINT *,' The LAPACK subroutine DGBTRS has'
          PRINT *,' returned ',INFO,' as a value of INFO.'
          PRINT *,' Program aborted.'
          STOP      
        ENDIF
C       .
        GOTO 10
C       .
      ENDIF
C     .
C     . Either there is convergence or an error has occured.
C     . Call DNEUPD.
C     .
      CALL DNEUPD ( EVECF, HOWMNY, SELECT, DR, DI, V, N2, 
     1              DRSV, DISV, WORKEV, BMAT, N2, WHICH, NEV,
     2              ARTOL, RESID, NCV, V, N2, IPARAM, IPNTR,
     3              WORKD, WORKL, LWORKL, INFO )
C     .
C     . Check call to DNEUPD has been successful.
C     .
      IF ( INFO.NE.0 ) THEN
        PRINT *,' Subroutine MFGEPS.'
        PRINT *,' The ARPACK subroutine DNEUPD has'
        PRINT *,' returned ',INFO,' as a value of INFO.'
        PRINT *,' Program aborted.'
        STOP       
      ENDIF
C     .
      NCE   = IPARAM( 5 )
      FIRST = .TRUE.
C     .
C     . Let's transform eigenvalues to those of original problem
C     . Also, 'complete' the vector ready to calculate
C     . direct residuals ...
C     .
      DO J = 1, NCE
        FACTR   = DLAPY2( DR( J ), DI( J ) )
        FACTR   = FACTR*FACTR
        DR( J ) = DRSV + DR( J )/FACTR
        DI( J ) = (-1.0d0)*DI( J )/FACTR
C       .
        CALL ASVCPL( V( 1, J ), NR, NDCS, INARR, MHP, MHIBC,
     1               MHOBC, NCFM, NDRVS, NDRVM, NBN, SVFDC )
C       .
      ENDDO
C     .
C     . Calculate direct residuals
C     .
      DO J = 1, NCE
C       .
C       . Zero the workspace array W2.
C       .
        IOP = 0
        CALL VECOP( W2, DZERO, N2, IOP )
C       .
C       . First the case of real eigenvalues ...
C       .
        IF ( DI( J ).EQ.DZERO ) THEN
          FACTR = DR( J ) - DRSV
C         .
C         . V( 1, j ) contains eigenvector so
C         . calculate, in W2, eval*B*x
C         .
          FAC = FACTR*CMD*(-1.0d0)
          ICMP = 4
          CALL ASVLP( NR, NDCS, V( 1, J ), INARR, MHT, MHL, MHM,
     1              ICMP, MHP, W2, INARR, MHT, MHL, MHM, FAC, NBN,
     2              NDRVS, NDRVM, ILN, IRN, SVFDC, XARR, NCFM )
C         .
          ICMP = 5
          CALL ASVLP( NR, NDCS, V( 1, J ), INARR, MHT, MHL, MHM,
     1              ICMP, MHP, W2, INARR, MHT, MHL, MHM, FAC, NBN,
     2              NDRVS, NDRVM, ILN, IRN, SVFDC, XARR, NCFM )
C         .
C         . W2 now contains eval*B*x: call DGBTRS to solve
C         .
          CALL DGBTRS( TRANS, N2, KL, KL, NRHS, A, N1,
     1                 IPIV, W2, N2, INFO )
C         .
          IF ( INFO.NE.0 ) THEN
            PRINT *,' Subroutine MFGEPS.'
            PRINT *,' The LAPACK subroutine DGBTRS has'
            PRINT *,' returned ',INFO,' as a value of INFO.'
            PRINT *,' Program aborted.'
            STOP
          ENDIF
C         .
C         . Complete newly solved vector
C         .
          CALL ASVCPL( W2, NR, NDCS, INARR, MHP, MHIBC,
     1               MHOBC, NCFM, NDRVS, NDRVM, NBN, SVFDC )
C         .
C         . Subtract V( 1, J ) from W2
C         .
          CALL DAXPY( N2, DMONE, V( 1, J ), IONE, W2, IONE )
C         .
          D3( J ) = DNRM2( N2, W2, IONE )
C         .
        ELSE
C         .
C         . Now the case of complex eigenvalues ...
C         .
          IF ( FIRST ) THEN
            FACTR = DR( J ) - DRSV
            FACTI = DI( J )
C           .
            FAC = FACTR*CMD*(-1.0d0)
            ICMP = 4
            CALL ASVLP( NR, NDCS, V( 1, J ), INARR, MHT, MHL, MHM,
     1              ICMP, MHP, W2, INARR, MHT, MHL, MHM, FAC, NBN,
     2              NDRVS, NDRVM, ILN, IRN, SVFDC, XARR, NCFM )
C           .
            ICMP = 5
            CALL ASVLP( NR, NDCS, V( 1, J ), INARR, MHT, MHL, MHM,
     1              ICMP, MHP, W2, INARR, MHT, MHL, MHM, FAC, NBN,
     2              NDRVS, NDRVM, ILN, IRN, SVFDC, XARR, NCFM )
C           .
            FAC = FACTI*CMD
            ICMP = 4
            CALL ASVLP( NR, NDCS, V( 1, J+1 ), INARR, MHT, MHL, MHM,
     1              ICMP, MHP, W2, INARR, MHT, MHL, MHM, FAC, NBN,
     2              NDRVS, NDRVM, ILN, IRN, SVFDC, XARR, NCFM )
C           .
            ICMP = 5
            CALL ASVLP( NR, NDCS, V( 1, J+1 ), INARR, MHT, MHL, MHM,
     1              ICMP, MHP, W2, INARR, MHT, MHL, MHM, FAC, NBN,
     2              NDRVS, NDRVM, ILN, IRN, SVFDC, XARR, NCFM )
C           .
C           . W2 now contains eval*B*x: call DGBTRS to solve
C           .
            CALL DGBTRS( TRANS, N2, KL, KL, NRHS, A, N1,
     1                   IPIV, W2, N2, INFO )
C           .
            IF ( INFO.NE.0 ) THEN
              PRINT *,' Subroutine MFGEPS.'
              PRINT *,' The LAPACK subroutine DGBTRS has'
              PRINT *,' returned ',INFO,' as a value of INFO.'
              PRINT *,' Program aborted.'
              STOP
            ENDIF
C           .
C           . Complete newly solved vector
C           .
            CALL ASVCPL( W2, NR, NDCS, INARR, MHP, MHIBC,
     1               MHOBC, NCFM, NDRVS, NDRVM, NBN, SVFDC )
C           .
C           . Subtract V( 1, J ) from W2
C           .
            CALL DAXPY( N2, DMONE, V( 1, J ), IONE, W2, IONE )
C           .
            D3( J ) = DNRM2( N2, W2, IONE )
C           .
C           . Zero the vector W2
C           .
            IOP = 0
            CALL VECOP( W2, DZERO, N2, IOP )
C           .
            FAC = FACTR*CMD*(-1.0d0)
            ICMP = 4
            CALL ASVLP( NR, NDCS, V( 1, J+1 ), INARR, MHT, MHL, MHM,
     1              ICMP, MHP, W2, INARR, MHT, MHL, MHM, FAC, NBN,
     2              NDRVS, NDRVM, ILN, IRN, SVFDC, XARR, NCFM )
C           .
            ICMP = 5
            CALL ASVLP( NR, NDCS, V( 1, J+1 ), INARR, MHT, MHL, MHM,
     1              ICMP, MHP, W2, INARR, MHT, MHL, MHM, FAC, NBN,
     2              NDRVS, NDRVM, ILN, IRN, SVFDC, XARR, NCFM )
C           .
            FAC = FACTI*CMD*(-1.0d0)
            ICMP = 4
            CALL ASVLP( NR, NDCS, V( 1, J ), INARR, MHT, MHL, MHM,
     1              ICMP, MHP, W2, INARR, MHT, MHL, MHM, FAC, NBN,
     2              NDRVS, NDRVM, ILN, IRN, SVFDC, XARR, NCFM )
C           .
            ICMP = 5
            CALL ASVLP( NR, NDCS, V( 1, J ), INARR, MHT, MHL, MHM,
     1              ICMP, MHP, W2, INARR, MHT, MHL, MHM, FAC, NBN,
     2              NDRVS, NDRVM, ILN, IRN, SVFDC, XARR, NCFM )
C           .
C           . W2 now contains eval*B*x: call DGBTRS to solve
C           .
            CALL DGBTRS( TRANS, N2, KL, KL, NRHS, A, N1,
     1                   IPIV, W2, N2, INFO )
C           .
            IF ( INFO.NE.0 ) THEN
              PRINT *,' Subroutine MFGEPS.'
              PRINT *,' The LAPACK subroutine DGBTRS has'
              PRINT *,' returned ',INFO,' as a value of INFO.'
              PRINT *,' Program aborted.'
              STOP
            ENDIF
C           .
C           . Complete newly solved vector
C           .
            CALL ASVCPL( W2, NR, NDCS, INARR, MHP, MHIBC,
     1               MHOBC, NCFM, NDRVS, NDRVM, NBN, SVFDC )
C           .
C           . Subtract V( 1, J+1 ) from W2
C           .
            CALL DAXPY( N2, DMONE, V( 1, J+1 ), IONE, W2, IONE )
C           .
            FAC       = DNRM2( N2, W2, IONE )
            D3( J )   = DLAPY2( D3( J ), FAC )
            D3( J+1 ) = D3( J )
            FIRST     = .FALSE.
C           .
          ELSE
            FIRST     = .TRUE.
          ENDIF
C         .
        ENDIF
C       .
      ENDDO
C     .
      RETURN
      END
C*********************************************************************

