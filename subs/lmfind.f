C********************************************************************
C subroutine L and M FIND *******************************************
C            -     - ---- *******************************************
C Steve Gibbons 12.4.97                                             C
C Modified 12.6.97 to have IT = 1 for cos and 2 for sine            C
C___________________________________________________________________C
C Given a number of harmonic, N; LMFIND will return the level, L,   C
C the order, M, and IT - which is equal to 2 for sin                C
C spherical harmonics and 1 for cosine ones.                        C
C All the above are integers - no point in a variable list .....    C
C     N = L*L for M = 0, IT = 1                                     C
C     N = L*L + 2*M - 1 for non-zero M and IT = 1                   C
C     N = L*L + 2*M for non-zero M and IT = 2                       C
C If N = 0, we have a monopole: L = 0, M = 0 and IT = 1.            C
C___________________________________________________________________C
C
C********************************************************************
      SUBROUTINE LMFIND ( N, L, M, IT)
      IMPLICIT NONE
C___________________________________________________________________C
C Variable Declarations - Parameters ...............................C
      INTEGER N,L,M,IT
C___________________________________________________________________C
C Variable declarations - Working Variables ........................C
      INTEGER LL,IDIFF,N2,ITWO
      PARAMETER (ITWO=2)
C___________________________________________________________________C
C First put N into N2 so that N is not altered
      N2=N
      IF ( N2.EQ.0 ) THEN
         L=0
         M=0
         IT=1
         RETURN
      ENDIF
      IF ( N2.LT.1 ) THEN 
         WRITE (6,989)
         STOP
 989  FORMAT (' Subroutine LMFIND. N<1 - Program stopped.')
      ENDIF
      IF ( N2.EQ.1 ) THEN
         L=1
         M=0
         IT=1
         RETURN
      ENDIF
      L=1
 500  CONTINUE
      LL=L*L
      IDIFF = N2 - LL
      IF ( IDIFF.GT.0 ) THEN
         L=L+1
         GOTO 500
      ENDIF
      IF ( IDIFF.EQ.0 ) THEN
         M=0
         IT=1
         RETURN
      ENDIF
      L=L-1
      LL=L*L
      N2=N2-LL
C so we know now that N2 is equal to either 2*M or 2*M-1
C corresponding to IT=1 and IT=2 respectively
      IDIFF = MOD ( N2, ITWO )
      IF ( IDIFF.EQ.1) THEN
         IT = 1
         M = (N2+1)/ITWO
      ELSE
         IT = 2
         M = N2/ITWO
      ENDIF
      RETURN
      END
C********************************************************************
