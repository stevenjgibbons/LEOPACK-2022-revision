C*********************************************************************
C subroutine Euler Angles of Rotation Cosine Matrix Calculate ********
C            -     -         -        -      -      -         ********
C Steve Gibbons Sat Mar 17 14:47:55 MET 2001                         C
C____________________________________________________________________C
C                                                                    C
C Here are the definitions of our Euler angles for rotating a        C
C coordinate axis system.                                            C
C                                                                    C
C (1) x_1', x_2' and x_3' axes are rotated about the x_3 axis        C
C through an angle ALPHA, anti-clockwise relative to x_1, x_2 and    C
C x_3 ( x_3 and x_3' coincide ).                                     C
C                                                                    C
C (2) x_1'', x_2'' and x_3'' axes are rotated about the x_2' axis    C
C through an angle BETA, anti-clockwise relative to x_1', x_2'       C
C and x_3' ( x_2' and x_2'' coincide ).                              C
C                                                                    C
C (3) Rotate about the x_3'' axis through an angle GAMMA -           C
C this yields x_1''', x_2''' and x_3''' (x_3'' and x_3''' coincide). C
C                                                                    C
C                                                                    C
C____________________________________________________________________C
C                                                                    C
C Input variables :-                                                 C
C ===============                                                    C
C                                                                    C
C  Double Precision                                                  C
C  ----------------                                                  C
C                                                                    C
C     ALPHAD    : Angle alpha in degrees.                            C
C     BETAD     : Angle beta in degrees.                             C
C     GAMMAD    : Angle gamma in degrees.                            C
C                                                                    C
C Output variables :-                                                C
C ================                                                   C
C                                                                    C
C  Double Precision                                                  C
C  ----------------                                                  C
C                                                                    C
C     EARCM     : Euler angle matrix Dim( 3, 3 ).                    C
C                                                                    C
C____________________________________________________________________C
C
C*********************************************************************
      SUBROUTINE EARCMC( ALPHAD, BETAD, GAMMAD, EARCM )
      IMPLICIT NONE
C____________________________________________________________________C
C Variable declarations - Parameters ................................C
      DOUBLE PRECISION ALPHAD, BETAD, GAMMAD, EARCM( 3, 3 )
C____________________________________________________________________C
C Variable declarations - Working variables .........................C
      DOUBLE PRECISION ALPHAR, BETAR, GAMMAR, PI,
     1                 COSA, SINA, COSB, SINB, COSG, SING
      PARAMETER ( PI=3.14159265358979312D0 )
C____________________________________________________________________C
C START OF PROGRAM **************************************************C
C____________________________________________________________________C
C
      ALPHAR = ALPHAD*PI/180.0d0
      BETAR  = BETAD*PI/180.0d0
      GAMMAR = GAMMAD*PI/180.0d0
C
      COSA   = DCOS( ALPHAR )
      SINA   = DSIN( ALPHAR )
      COSB   = DCOS( BETAR  )
      SINB   = DSIN( BETAR  )
      COSG   = DCOS( GAMMAR )
      SING   = DSIN( GAMMAR )
C
      EARCM( 1, 1 ) = COSG*COSB*COSA - SING*SINA
      EARCM( 1, 2 ) = COSG*COSB*SINA + SING*COSA
      EARCM( 1, 3 ) = COSG*SINB*(-1.0d0)
C
      EARCM( 2, 1 ) = (-1.0d0)*SING*COSB*COSA - COSG*SINA
      EARCM( 2, 2 ) = COSG*COSA - SING*COSB*SINA
      EARCM( 2, 3 ) = SING*SINB
C
      EARCM( 3, 1 ) = SINB*COSA
      EARCM( 3, 2 ) = SINB*SINA
      EARCM( 3, 3 ) = COSB
C
      RETURN
      END
C*********************************************************************
