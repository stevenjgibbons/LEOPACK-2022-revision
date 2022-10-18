C*********************************************************************
C THE NUMERICAL RECIPES NUMERICAL INTEGRATION ROUTINES ADAPTED FOR ***
C DOUBLE PRECISION AND A GENERALISED FUNCTION OF THE FORM ....       *
C								     *
C FUNC ( X ,INTPAR, DPPAR )					     *
C where INTPAR and DPPAR are respectively INTEGER and DOUBLE         *
C PRECISION arrays with dimension declaration (*).                   *
C								     *
C These codes are (C) Copr. 1986-92 Numerical Recipes Software       *
C Adapted by Steve Gibbons 8.5.97.                                   *
C*********************************************************************
C calling sequence:-						     *
C  CALL NUMINT ( FUNC, X1, X2, DPINT, LOW, INTPAR, DPPAR )           *
C____________________________________________________________________C
C Input variables :-                                                 C
C ===============                                                    C
C  Double Precision                                                  C
C  ----------------                                                  C
C     FUNC	: Declared EXTERNALly. Has the form                  C
C		  FUNC( X, INTPAR, DPPAR )			     C
C     X1	: Lower limit of integration.                        C
C     X2	: Upper limit of integration.                        C
C     DPINT	: The integral.					     C
C     LOW	: The accepted difference between successive         C
C                  integrations. Suggested value about 1.0d-7        C
C     DPPAR	: Double Precision array of parameters for FUNC.     C
C                  with dimension ( * ).                             C
C  Integer                                                           C
C  -------                                                           C
C     INTPAR	: Integer array of parameters for FUNC. Dim (*)      C
C____________________________________________________________________C
C Calls the numerical recipes routine TRAPZD.                        C
C____________________________________________________________________C
C
C*********************************************************************
      SUBROUTINE NUMINT ( FUNC, X1, X2, DPINT, LOW, INTPAR, DPPAR )
      IMPLICIT NONE
C____________________________________________________________________C
C Variable declarations - Parameters ................................C
      DOUBLE PRECISION FUNC, X1, X2, DPINT, LOW, DPPAR( * )
      EXTERNAL FUNC
      INTEGER INTPAR( * )
C____________________________________________________________________C
C Variable declarations - Working variables .........................C
      INTEGER J, JMAX
      DOUBLE PRECISION ST,OS,OST
      PARAMETER ( JMAX = 20 )
C____________________________________________________________________C
C START OF PROGRAM **************************************************C
C____________________________________________________________________C
C Check validity of arguments ....
      IF ( LOW.LT.1.0d-20 .OR. LOW.GT.1.0d-4 ) THEN
         PRINT *,' Subroutine NUMINT.'
         PRINT *,' LOW = ',LOW,'. Should be approx. 1.0d-7 '
         PRINT *,' Program aborted.'
         STOP
      ENDIF
C
      OST = -1.0d30
      OS = -1.0d30
C
      DO J = 1, JMAX
        CALL TRAPZD( FUNC, X1, X2, ST, J , INTPAR, DPPAR )
        DPINT = (4.0d0*ST-OST)/3.0d0
        IF ( DABS( DPINT - OS ).LT.LOW*DABS(OS) ) RETURN
        OS = DPINT
        OST = ST
      ENDDO
C ........ integration hasn't worked;
      PRINT *,' Subroutine NUMINT.'
      PRINT *,' Integration failed.'
      PRINT *,' Program aborted.'
      STOP
C
      END
C*********************************************************************
