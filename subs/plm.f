C*********************************************************************
C function PLM *******************************************************
C Steve Gibbons 16.4.97                                              C
C____________________________________________________________________C
C Calculates the Schmidt Normalised Legendre Function P_l^m (x)      C
C given P_(l-1)^m and P_(l-2)^m according to equation 183 in my notesC
C i.e.                                                               C
C   P_l^m( X ) = { A * P_(l-1)^m - B * P_(l-2)^m }/C                 C
C                                                                    C
C where A = (2*l - 1)*X ,                                            C
C                                                                    C
C B = SQRT( (L+M-1)*(L-M-1) ) and C = SQRT( (L+M)*(L-M) )            C
C                                                                    C
C____________________________________________________________________C
C Input Variables :-                                                 C
C ===============                                                    C
C  Integer                                                           C
C  -------                                                           C
C     M         : Well it's M isn't it!                              C
C     L         : Well it's L isn't it!                              C
C  Double Precision                                                  C
C  ----------------                                                  C
C     X         : Cos(theta)                                         C
C     PLMIN1    : P_(l-1)^m ( X )                                    C
C     PLMIN2    : P_(l-2)^m ( X )                                    C
C____________________________________________________________________C
C Local Variables :-                                                 C
C ===============                                                    C
C  Integer                                                           C
C  -------                                                           C
C     RL        : DBLE  ( L )                                        C
C     RM        : DBLE  ( M )                                        C
C____________________________________________________________________C
C
C*********************************************************************
      FUNCTION PLM ( L, M, X, PLMIN1, PLMIN2 )
      IMPLICIT NONE
      DOUBLE PRECISION PLM
C____________________________________________________________________C
C Variable declarations - Parameters ................................C
      INTEGER M,L
      DOUBLE PRECISION X, PLMIN1, PLMIN2
C____________________________________________________________________C
C Variable declarations - Working variables .........................C
      DOUBLE PRECISION RM,RL
C____________________________________________________________________C
C START OF FUNCTION *************************************************C
C____________________________________________________________________C
C
      IF ( L.LT.2 ) THEN
         PRINT *,' You are trying to run function PLM with L < 2.'
         PRINT *,' Program aborted.'
         STOP
      ENDIF
      IF ( M.LT.0 .OR. M.GT.L ) THEN
         PRINT *,' M is out of range in function PLM.'
         PRINT *,' Program aborted.'
         STOP
      ENDIF
      IF ( L.EQ.M ) THEN
         PRINT *,' PLM function called with L = M .'
         PRINT *,' Division by zero would follow.'
         PRINT *,' Program aborted.'
         STOP
      ENDIF
      RM = DBLE( M )
      RL = DBLE( L )
      PLM = ( 2.0d0*RL - 1.0d0 )*X*PLMIN1
      PLM = PLM - PLMIN2*DSQRT( ( RL+RM-1.0d0 )*( RL-RM-1.0d0 ) )
      PLM = PLM/DSQRT( (RL+RM)*(RL-RM) )
      RETURN
      END
C*********************************************************************
