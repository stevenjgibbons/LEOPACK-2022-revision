C*********************************************************************
C function CHandrasekhar's ORthogonal CHaracteristic function ********
C          --              --         --                      ********
C Steve Gibbons Thu Oct  7 15:00:53 BST 1999                         C
C____________________________________________________________________C
C                                                                    C
C If IFLAG = INTARR( 1 ) = 1 then                                    C
C              CHORCH = TANH( 0.5*X ) + TAN( 0.5*X )                 C
C                                                                    C
C If IFLAG = INTARR( 1 ) = 2 then                                    C
C              CHORCH = 1.0d0/TANH( 0.5*X ) - 1.0d0/TAN( 0.5*X )     C
C                                                                    C
C DPRARR is not referred to by CHORCH. It is merely there so that    C
C roots to this equation may be found by SMITSL.                     C
C                                                                    C
C____________________________________________________________________C
C
C*********************************************************************
      FUNCTION CHORCH( X, INTARR, DPRARR )
      IMPLICIT NONE
C____________________________________________________________________C
C Variable declarations - Parameters ................................C
      DOUBLE PRECISION CHORCH, DPRARR( * ), X
      INTEGER INTARR( * )
C____________________________________________________________________C
C Variable declarations - Working variables .........................C
      INTEGER IFLAG
C____________________________________________________________________C
C START OF PROGRAM **************************************************C
C____________________________________________________________________C
C
      IFLAG = INTARR( 1 )
C     .
C     . Check IFLAG is valid
C     .
      IF ( IFLAG.NE.1 .AND. IFLAG.NE.2 ) THEN
        PRINT *,' Function CHORCH.'
        PRINT *,' INTARR( 1 ) = ', IFLAG
        PRINT *,' Program aborted.'
        STOP
      ENDIF
C     .
      IF ( IFLAG.EQ.1 ) CHORCH = TANH( 0.5*X ) + TAN( 0.5*X )
      IF ( IFLAG.EQ.2 ) CHORCH = 1.0d0/TANH( 0.5*X ) -
     1                 1.0d0/TAN( 0.5*X )
C     .
      RETURN
      END
C*********************************************************************
