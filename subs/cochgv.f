C*********************************************************************
C subroutine COefficients for CHebyshev basis (Gibbons) version ******
C            --               --               -                ******
C Steve Gibbons Fri Apr 28 12:38:41 BST 2000                         C
C                                                                    C
C This routine is almost entirely derived from COCHGA written        C
C by D. Funaro of Consiglio Nazionale Delle Ricerche,                C
C Via Abbiategrasso, 209 - 27100 Pavia, Italy.                       C
C                                                                    C
C I have merely changed the IMPLICIT DOUBLE PRECISION statement to   C
C an implicit none statement and have declared all subsequent        C
C variables.                                                         C
C                                                                    C
C If QZ( i ) contains the value f( x_i ) where x_i is the i^{th}     C
C zero of the Chebyshev polynomial T_N( x ) with x in the interval   C
C ( -1, 1 ), then CO( k ) contains the Chebshev coefficient c_k      C
C where                                                              C
C                                                                    C
C f( x ) = \sum_{k=0}^{N-1} c_k T_k( x )                             C
C                                                                    C
C The array CO is indexed (0:N-1) and not (1:N)                      C
C                                                                    C
C____________________________________________________________________C
C                                                                    C
C Input variables :-                                                 C
C ===============                                                    C
C                                                                    C
C  Integer                                                           C
C  -------                                                           C
C                                                                    C
C     N         : Degree of Chebshev polynomial whose zeros define   C
C                 the x values.                                      C
C                                                                    C
C  Double Precision                                                  C
C  ----------------                                                  C
C                                                                    C
C     QZ        : QZ( i ) is value of function at T_N zero x_i.      C
C     CO        : CO( k ) is output as coeff. c_k (k=0,n-1).         C
C                                                                    C
C____________________________________________________________________C
C
C*********************************************************************
      SUBROUTINE COCHGV( N, QZ, CO )
C*********************************************************************
C Following five lines are Funaro notes.
C   COMPUTES THE CHEBYSHEV FOURIER COEFFICIENTS OF A POLYNOMIAL
C   INDIVIDUATED BY THE VALUES ATTAINED AT THE CHEBYSHEV ZEROES
C   N  = THE NUMBER OF ZEROES
C   QZ = VALUES OF THE POLYNOMIAL AT THE ZEROES, QZ(I), I=1,N
C   CO = FOURIER COEFFICIENTS OF THE POLYNOMIAL, CO(I), I=0,N-1
C*********************************************************************
      IMPLICIT NONE
      INTEGER N
      DOUBLE PRECISION QZ(N), CO(0:N-1)
C____________________________________________________________________C
C Variable declarations - Working variables .........................C
      INTEGER K, J
      DOUBLE PRECISION PH, DN, SU, SK, DJ, DK
C____________________________________________________________________C
C START OF PROGRAM **************************************************C
C____________________________________________________________________C
C
      IF (N .EQ. 0) RETURN
      PH = 1.57079632679489661923D0
      DN = DFLOAT(N)
      SU = 0.0D0
      DO 1 J = 1, N
          SU = SU+QZ(J)
 1    CONTINUE
          CO(0) = SU/DN
      IF (N .EQ. 1) RETURN
          SK = -2.0D0
      DO 2 K=1,N-1
          DK = DFLOAT(K)
          SU = 0.0D0
      DO 3 J = 1, N
          DJ = 2.0D0*DFLOAT(J)-1.0D0
          SU = SU+QZ(J)*DCOS(DK*DJ*PH/DN)
 3    CONTINUE
          CO(K) = SK*SU/DN
          SK = -SK
 2    CONTINUE
C
      RETURN
      END
C*********************************************************************
