C*********************************************************************
C subroutine EigenVALue Array Sort ***********************************
C            -    ---   -     -    ***********************************
C Steve Gibbons Sun Nov 21 12:39:31 GMT 1999                         C
C____________________________________________________________________C
C                                                                    C
C Subroutines such as VMEPS which utilise the Arnoldi eigensolvers   C
C return three double precision arrays containing the eigenvalue     C
C information.                                                       C
C                                                                    C
C  They are; DR which contains the real parts of the eigenvalues,    C
C            DI which contains the imag. parts of the eigenvalues,   C
C       and  D3  which contains their direct residuals.              C
C                                                                    C
C  The purpose of EVALAS is to sort through these arrays and         C
C  return as single variables:                                       C
C                                                                    C
C   GRR - the maximum value of the real part.                        C
C   GRI - the corresponding value of the imag. part.                 C
C                                                                    C
C   IEV is returned with the number of eigenvalue with the           C
C   largest real part.                                               C
C                                                                    C
C  If the direct residual for any given eigenvalue is greater        C
C  than a given tolerance value, RESTOL, then an error message is    C
C  returned, implying that there has been some problem with the      C
C  solution problem.                                                 C
C                                                                    C
C  D3 is overwritten on exit with work information.                  C
C                                                                    C
C  DRSV is the real shift applied in the solution routine.           C
C                                                                    C
C  EVALAS seeks to find if any of the other eigenvalues (whose       C
C  real part was not the greatest) are closer in the complex plane   C
C  to the point ( DRSV, 0.0 ).                                       C
C                                                                    C
C  If so, then it is possible that a further solution with modified  C
C  parameters will miss this mode and find a lower mode which,       C
C  although of smaller real part, is easier to find because of a     C
C  smaller imaginary part.                                           C
C                                                                    C
C  In this case, DSRSV (the suggested shift value) is returned       C
C  conaining DABS( GRI ).                                            C
C                                                                    C
C                                                                    C
C____________________________________________________________________C
C                                                                    C
C  Integer                                                           C
C  -------                                                           C
C                                                                    C
C     NCE       : Number of converged eigenvalues.                   C
C     IEV       : [Output] number of eigenvalue with largest Re part C
C                 Returned with -1 if NCE is less than 1.            C
C                                                                    C
C  Double Precision                                                  C
C  ----------------                                                  C
C                                                                    C
C     DRSV      : Real shift applied when solving eigensystem.       C
C     RESTOL    : Tolerance level for direct residuals.              C
C     DR        : Dim ( NCE ). Real parts of eigenvalues.            C
C     DI        : Dim ( NCE ). Imag. parts of eigenvalues.           C
C     D3        : Dim ( NCE ). Direct residuals.                     C
C                                                                    C
C     GRR       : Largest real part of eigenvalues.                  C
C     GRI       : Corresponding imaginary part of eigenvalue.        C
C                                                                    C
C     DSRSV     : Suggested real shift to be applied when solving    C
C                   eigensystem.                                     C
C                                                                    C
C____________________________________________________________________C
C
C*********************************************************************
      SUBROUTINE EVALAS( NCE, IEV, DRSV, RESTOL, DR, DI, D3,
     1                   GRR, GRI, DSRSV )
      IMPLICIT NONE
C____________________________________________________________________C
C Variable declarations - Parameters ................................C
      INTEGER NCE, IEV
      DOUBLE PRECISION DR( * ), DI( * ), D3( * ), GRR, GRI, DSRSV,
     1                 DRSV, RESTOL
C____________________________________________________________________C
C Variable declarations - Working variables .........................C
C
      INTEGER IE
      DOUBLE PRECISION ZERO, DX, DY, DIST, DLAPY2
      LOGICAL FIRST
      PARAMETER ( ZERO = 0.0d0 )
C____________________________________________________________________C
C START OF PROGRAM **************************************************C
C____________________________________________________________________C
C
      GRR = -1.0d8
      IEV = -1
      IF ( NCE.LT.1 ) RETURN
C     .
      IF ( RESTOL.LT.ZERO ) THEN
        PRINT *,' Subroutine EVALAS.'
        PRINT *,' RESTOL = ', RESTOL
        PRINT *,' Program aborted.'
        STOP
      ENDIF
C     .
      FIRST = .TRUE.
      DO IE = 1, NCE
        IF ( DI( IE ).EQ.ZERO ) THEN
C         .
C         . Case of real eigenvalues
C         . First: check direct residual
C         .
          IF ( D3( IE ).GT.RESTOL ) THEN
            PRINT *,' Subroutine EVALAS.'
            PRINT *,' D3(',IE,') = ', D3( IE )
            PRINT *,' Probable fault in solution of eigenproblem.'
            PRINT *,' Program aborted.'
            STOP
          ENDIF
C         .
          IF ( DR( IE ).GT.GRR ) THEN
            GRR      = DR( IE )
            GRI      = 0.0d0
            IEV      = IE
          ENDIF
          DY         = DRSV - DR( IE )
          DIST       = DABS( DY )
          D3( IE )   = DIST
C         .
        ELSE
C         .
C         . Case of complex eigenvalues
C         . First: check direct residual
C         .
          IF ( FIRST ) THEN
C          .
           IF ( D3( IE ).GT.RESTOL .OR. D3(IE+1).GT.RESTOL ) THEN
             PRINT *,' Subroutine EVALAS.'
             PRINT *,' D3(',IE,') = ', D3( IE )
             PRINT *,' D3(',IE+1,') = ', D3( IE+1 )
             PRINT *,' Probable fault in solution of eigenproblem.'
             PRINT *,' Warning only ... but be aware of this.'
c            PRINT *,' Program aborted.'
c            STOP
           ENDIF
C          .
           IF ( DR( IE ).GT.GRR ) THEN
             GRR      = DR( IE )
             GRI      = DI( IE )
             IEV      = IE
           ENDIF
           DY         = DRSV - DR( IE )
           DX         = DI( IE )
           DIST       = DLAPY2( DX, DY )
           D3( IE )   = DIST
           D3( IE+1 ) = DIST + 1.0d0
C          .
           FIRST = .FALSE.
          ELSE
           FIRST = .TRUE.
          ENDIF
C         .
        ENDIF
      ENDDO
C     .
      DIST  = D3( IEV )
      DSRSV = DRSV
      DO IE = 1, NCE
        IF ( IE.NE.IEV .AND. D3( IE ).LT.DIST ) DSRSV = DABS( GRI )
      ENDDO
C     .
      RETURN
      END
C*********************************************************************

