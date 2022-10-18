C*********************************************************************
C subroutine Toroidal Solid Body Rotation Vector Form ****************
C            -        -     -    -        -      -    ****************
C Steve Gibbons Wed Oct 31 13:40:16 WET 2001                         C
C____________________________________________________________________C
C                                                                    C
C Use the routine TSBRVA to add the vector TSBRV to VT2.             C
C                                                                    C
C____________________________________________________________________C
C                                                                    C
C Input variables :-                                                 C
C ===============                                                    C
C                                                                    C
C  Integer                                                           C
C  -------                                                           C
C                                                                    C
C     NR        : Number of radial grid nodes in main soln.          C
C     NH2       : Number of toroidal harmonic functions.             C
C     ML2       : Array dim ( NH2 ). Sph. harm degree, L.            C
C     MM2       : Array dim ( NH2 ). Sph. harm order, M, or -M.      C
C            ( ml2 and mm2 describe toroidal velocity )              C
C                                                                    C
C  Double Precision                                                  C
C  ----------------                                                  C
C                                                                    C
C     TSBRV     : Dim (NR). Toroidal Solid Body Rotation Vector.     C
C     VT2       : Dim (NR*NH2). Toroidal velocity function.          C
C                                                                    C
C____________________________________________________________________C
C
C*********************************************************************
      SUBROUTINE TSBRVA( NR, NH2, ML2, MM2, TSBRV, VT2 )
      IMPLICIT NONE
C____________________________________________________________________C
C Variable declarations - Parameters ................................C
      INTEGER          NR, NH2, ML2( NH2 ), MM2( NH2 )
      DOUBLE PRECISION TSBRV( NR ), VT2( NR*NH2 )
C____________________________________________________________________C
C Variable declarations - Working variables .........................C
      INTEGER          IR, IND, IHARM, IND0, IH
C____________________________________________________________________C
C START OF PROGRAM **************************************************C
C____________________________________________________________________C
C
      IHARM = 0
      DO IH = 1, NH2
        IF ( ML2( IH ).EQ.1 .AND. MM2( IH ).EQ.0 ) THEN
          IHARM = IH
          GOTO 500
        ENDIF
      ENDDO
      IF ( IHARM.EQ.0 ) THEN
        PRINT *,' Subroutine TSBRVA.'
        PRINT *,' T_1^0 harmonic not found.'
        PRINT *,' Program aborted.'
        STOP
      ENDIF
C
 500  CONTINUE
      IND0   = ( IHARM - 1 )*NR
      DO IR = 1, NR
        IND    = IND0 + IR
        VT2( IND ) = VT2( IND ) + TSBRV( IR )
      ENDDO
C
      RETURN
      END
C*********************************************************************
