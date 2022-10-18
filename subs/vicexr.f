C*********************************************************************
C subroutine Vector Interaction Coefficient EXtraction Routine *******
C            -      -           -           --         -       *******
C Steve Gibbons Mon Jan 17 10:19:39 GMT 2000                         C
C____________________________________________________________________C
C                                                                    C
C For a given 4 character code string REQVIT, index numbers          C
C IRAHN, IRBHN and IRGHN; VICEXR will search through the arrays      C
C IHNALP, IHNBET, IHNGAM and TVHI to look for a non-zero interaction C
C                                                                    C
C If for a given 'in' (interaction number),                          C
C                                                                    C
C   TVHI( in ) .EQ. REQVIT,                                          C
C   IHNALP( in ) .EQ. IRAHN,                                         C
C   IHNBET( in ) .EQ. IRBHN,                                         C
C   IHNGAM( in ) .EQ. IRGHN  then                                    C
C                                                                    C
C COEF is returned with CVI( in ) and ONZIC is returned .TRUE.       C
C                                                                    C
C Otherwise, it is assumed that there is no such non-zero interact.  C
C and ONZIC is returned .FALSE.                                      C
C                                                                    C
C____________________________________________________________________C
C                                                                    C
C Input variables :-                                                 C
C ===============                                                    C
C  Integer                                                           C
C  -------                                                           C
C                                                                    C
C     NVI      : Number of vector interactions.                      C
C     IHNALP   : Number of alpha harmonics. Dim ( * )                C
C     IHNBET   : Number of beta harmonics. Dim ( * )                 C
C     IHNGAM   : Number of gamma harmonics. Dim ( * )                C
C                                                                    C
C     IRAHN    : Number of requested alpha harmonic                  C
C     IRBHN    : Number of requested beta harmonic                   C
C     IRGHN    : Number of requested gamma harmonic                  C
C                                                                    C
C  Character                                                         C
C  ---------                                                         C
C                                                                    C
C     TVHI      : *(4) Type of vector interaction. Dim. ( * )        C
C     REQVIT    : Requested type of vector interaction.              C
C                                                                    C
C  Double Precision                                                  C
C  ----------------                                                  C
C                                                                    C
C     CVI       : Coefficients of vector interaction. Dim ( * )      C 
C     COEF      : Returned coeffcient.                               C
C                                                                    C
C  Logical                                                           C
C  -------                                                           C
C                                                                    C
C     ONZIC     : Returned .TRUE. if COEF contains the coefficient   C
C                 Returned .FALSE. if interaction was not found.     C
C____________________________________________________________________C
C
C*********************************************************************
      SUBROUTINE VICEXR( NVI, IHNALP, IHNBET, IHNGAM, IRAHN, IRBHN, 
     1                   IRGHN, TVHI, REQVIT, CVI, COEF, ONZIC )
      IMPLICIT NONE
C____________________________________________________________________C
C Variable declarations - Parameters ................................C
      INTEGER NVI, IHNALP( * ), IHNBET( * ), IHNGAM( * ), IRAHN,
     1        IRBHN, IRGHN
      CHARACTER *(4) TVHI( * ), REQVIT
      DOUBLE PRECISION CVI( * ), COEF
      LOGICAL ONZIC
C____________________________________________________________________C
C Variable declarations - Working variables .........................C
C
      INTEGER IN
C____________________________________________________________________C
C START OF PROGRAM **************************************************C
C____________________________________________________________________C
C     .
      ONZIC = .FALSE.
      IN    = 0
      COEF  = 0.0d0
C
 50   CONTINUE
      IN = IN + 1
      IF (    IN.GT.NVI    ) RETURN
      IF (    IRAHN.NE.IHNALP( IN )    ) GOTO 50
      IF (    IRBHN.NE.IHNBET( IN )    ) GOTO 50
      IF (    IRGHN.NE.IHNGAM( IN )    ) GOTO 50
      IF (    REQVIT.NE.TVHI( IN )     ) GOTO 50
C
C OK it seems we have our desired interaction
C
      COEF  = CVI( IN )
      ONZIC = .TRUE.
      RETURN
      END
C*********************************************************************

