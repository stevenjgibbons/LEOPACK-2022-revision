C*********************************************************************
C subroutine Chebyshev Polynomial Change Of Variable Routine *********
C            -         -          -      -  -        -       *********
C Steve Gibbons Thu Apr 27 16:30:28 BST 2000                         C
C____________________________________________________________________C
C                                                                    C
C If your variable x lies in the interval [a,b], then to be able     C
C to use Chebyshev polynomials, you need to transform to the         C
C variable y in the interval [-1,1].                                 C
C                                                                    C
C      x - 0.5(b+a)                                                  C
C y =  ------------                                                  C
C       0.5(b - a)                                                   C
C                                                                    C
C and so ofcourse                                                    C
C                                                                    C
C x = 0.5( b - a )y + 0.5(b+a)                                       C
C                                                                    C
C If IDIR = 1 then x is input and y is output, and if IDIR = -1      C
C then y is input and x is output.                                   C
C                                                                    C
C____________________________________________________________________C
C                                                                    C
C Input variables :-                                                 C
C ===============                                                    C
C  Integer                                                           C
C  -------                                                           C
C                                                                    C
C     IDIR      : Set to 1 for y as a function of x.                 C
C                 Set to -1 for x as a function of y.                C
C                                                                    C
C  Double Precision                                                  C
C  ----------------                                                  C
C                                                                    C
C     X         : Variable in interval [a,b]                         C
C     Y         : Variable in interval [-1,1]                        C
C     A         : Lower limit for x.                                 C
C     B         : Upper limit for x.                                 C
C                                                                    C
C                                                                    C
C____________________________________________________________________C
C
C*********************************************************************
      SUBROUTINE CPCOVR( IDIR, X, Y, A, B )
      IMPLICIT NONE
C____________________________________________________________________C
C Variable declarations - Parameters ................................C
      INTEGER IDIR
      DOUBLE PRECISION X, Y, A, B
C____________________________________________________________________C
C Variable declarations - Working variables .........................C
      DOUBLE PRECISION DLOW
      PARAMETER ( DLOW = 1.0d-7 )
C____________________________________________________________________C
C START OF PROGRAM **************************************************C
C____________________________________________________________________C
C
C Check input parameters ...
C
      IF ( IDIR.EQ.1 ) THEN
        IF ( X.LT.A .OR. X.GT.B .OR. DABS( B - A ).LT.DLOW ) THEN
          PRINT *,' Subroutine CPCOVR.'
          PRINT *,' A = ', A,' B = ', B,' X = ', X
          PRINT *,' Program aborted.'
          STOP
        ENDIF
        Y = 2.0d0*(X-0.5d0*(B+A))/(B-A)
        RETURN
      ENDIF
C
      IF ( IDIR.EQ.-1 ) THEN
        IF ( Y.LT.(-1.0d0) .OR. Y.GT.1.0d0 ) THEN
          PRINT *,' Subroutine CPCOVR.'
          PRINT *,' Y = ', Y
          PRINT *,' Program aborted.'
          STOP
        ENDIF
        X = 0.5d0*(B-A)*Y + 0.5d0*(B+A)
        RETURN
      ENDIF
C
      PRINT *,' Subroutine CPCOVR.'
      PRINT *,' IDIR = ', IDIR
      PRINT *,' Program aborted.'
      STOP
      END
C*********************************************************************
