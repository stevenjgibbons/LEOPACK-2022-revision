C*********************************************************************
C subroutine Zeroes of Chebyshev Polynomial Abscissa Allocation Sub. *
C            -         -         -          -        -          -    *
C Steve Gibbons Tue Sep 21 17:46:55 BST 1999                         C
C This routine is almost entirely derived from ZECHGA written        C
C by D. Funaro of Consiglio Nazionale Delle Ricerche,                C
C Via Abbiategrasso, 209 - 27100 Pavia, Italy.                       C
C However, rather than giving the N zeros of the Chebyshev polynom.  C
C of degree N in the interval (-1,1), it returns XARR( 1 ) = RI,     C
C XARR( NR ) = RO and XARR( i + 1 ) as the i^{th} zero of the        C
C Chebyshev polynomial of degree (NR-2) scaled into the interval     C
C (RI,RO) ie. x := (ro-ri)(x+1.0d0)/2.0d0 + ri                       C
C                                                                    C
C____________________________________________________________________C
C                                                                    C
C Input variables :-                                                 C
C ===============                                                    C
C  Double Precision                                                  C
C  ----------------                                                  C
C                                                                    C
C     XARR      : Array of dimension (  NR  ).                       C
C                  XARR( i ) contains the value of x or r at the     C
C                   i^{th} radial grid node.                         C
C                                                                    C
C     RI        : Radius of inner boundary                           C
C     RO        : Radius of outer boundary                           C
C                                                                    C
C  Integer                                                           C
C  -------                                                           C
C                                                                    C
C     NR        : Number of grid nodes.                              C
C                                                                    C
C____________________________________________________________________C
C
C*********************************************************************
      SUBROUTINE ZCPAAS( NR, XARR, RI, RO)
      IMPLICIT NONE
C____________________________________________________________________C
C Variable declarations - Parameters ................................C
      INTEGER NR
      DOUBLE PRECISION RI, RO, XARR( NR )
C____________________________________________________________________C
C Variable declarations - Working variables .........................C
      INTEGER N, N2, I, IN
      DOUBLE PRECISION PH, DN, C, SI, DI, CSX, RIPRO2,
     1                 RIMRO2
C____________________________________________________________________C
C START OF PROGRAM **************************************************C
C____________________________________________________________________C
C
      IF ( NR.LT.4 ) THEN
         PRINT *,' Subroutine ZCPAAS.'
         PRINT *,' NR = ', NR
         STOP
      ENDIF
C
      IF ( RI.GE.RO ) THEN
         PRINT *,' Subroutine ZCPAAS.'
         PRINT *,' RI = ', RI
         PRINT *,' RO = ', RO
         STOP
      ENDIF
C
      RIPRO2 = 0.5d0*(RO + RI)
      RIMRO2 = 0.5d0*(RO - RI)
C
      XARR(  1 ) = RI
      XARR( NR ) = RO
      N = NR - 2
      XARR(2) = 0.D0
      N2 = N/2
      IN = 1+4*N2-2*N
      PH = 1.57079632679489661923D0
      DN = DFLOAT(N)
      C  = PH/DN
      SI = -1.D0
      DO 10 I = 1, N2
         DI = DFLOAT(I)
         CSX = DCOS(C*(2.D0*DI-1.D0))
         XARR(I+1) = (1.0d0-CSX)*RIMRO2 + RI
         XARR(N-I+2) = (1.0d0+CSX)*RIMRO2 + RI
         SI = -SI
 10   CONTINUE
C
      IF (IN .EQ. 1) RETURN
      XARR(N2+2) = RIPRO2
C
      RETURN
      END
C*********************************************************************
