C*********************************************************************
C function Index of Spherical Harm. Coeff. from Index Array **********
C          -        -         -     -           -     -     **********
C Steve Gibbons Wed Feb  2 07:39:05 GMT 2000                         C
C____________________________________________________________________C
C                                                                    C
C Performs the action of the function INDSHC but without having      C
C to repeat the lines of code to convert MHL and MHM to L, M and ICS C
C                                                                    C
C____________________________________________________________________C
C                                                                    C
C Input variables :-                                                 C
C ===============                                                    C
C  Integer                                                           C
C  -------                                                           C
C     IH        : Number of current spherical harmonic in set.       C
C     MHL       : Dim ( * ). Array of sph. hrm. degree l values.     C
C     MHM       : Dim ( * ). Array of sph. hrm. orders, m, for cos   C
C                  phi dependency and (-m) for sin m phi depend.     C
C                                                                    C
C____________________________________________________________________C
C
C*********************************************************************
      FUNCTION ISHCIA( IH, MHL, MHM )
      IMPLICIT NONE
C____________________________________________________________________C
C Variable declarations - Parameters ................................C
C
      INTEGER ISHCIA, IH, MHL( * ), MHM( * )
C____________________________________________________________________C
C Variable declarations - Working variables .........................C
C
      INTEGER L, M, ICS, INDSHC
C____________________________________________________________________C
C START OF PROGRAM **************************************************C
C____________________________________________________________________C
C
      L  = MHL( IH )
      IF ( MHM( IH ).LT.0 ) THEN
        M    = -MHM( IH )
        ICS  = 2
      ELSE
        M    = MHM( IH )
        ICS  = 1
      ENDIF
      ISHCIA = INDSHC( L, M, ICS )
      RETURN
      END
C*********************************************************************
