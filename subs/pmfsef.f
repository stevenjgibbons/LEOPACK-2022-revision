C*********************************************************************
C function Poloidal Magnetic Field Spectral Expansion Function *******
C          -        -        -     -        -         -        *******
C Steve Gibbons Sun Oct 10 15:54:27 BST 1999                         C
C____________________________________________________________________C
C                                                                    C
C Auxiliary function to MNEWTR to solve for alpha and beta           C
C when F1( alpha, beta) = 0  and F2( alpha, beta) = 0.               C
C                                                                    C
C Here F1( alpha, beta) = -ri alpha sin ( alpha ri - beta )          C
C                         - l cos ( alpha ri - beta )                C
C                                                                    C
C Here F2( alpha, beta) = -ro alpha sin ( alpha ro - beta )          C
C                         + (l + 1 ) cos ( alpha ri - beta )         C
C                                                                    C
C Also gives partial derivatives of F1 and F2 with respect to        C
C alpha and beta.                                                    C
C                                                                    C
C____________________________________________________________________C
C Input variables :-                                                 C
C ===============                                                    C
C  Double Precision                                                  C
C  ----------------                                                  C
C     XVEC      : XVEC( 1 ) contains alpha.                          C
C                 XVEC( 2 ) contains beta.                           C
C                                                                    C
C     DPRARR    : DPRARR( 1 ) contains RI.                           C
C                 DPRARR( 2 ) contains RO. ( Array of dim ( * ) )    C
C                                                                    C
C  Integer                                                           C
C  -------                                                           C
C                                                                    C
C     NEQN      : Number of equations and unknowns. Must be 2.       C
C                                                                    C
C     IEQN      : Number of current equation. Must be either 1 or 2. C
C                 (Corresponds to F1 and F2 above resp.)             C
C                                                                    C
C     ID        : ID = 0 --> return value of F( IEQN )               C
C                 ID = 1 --> return value of dF( IEQN )/d alpha      C
C                 ID = 2 --> return value of dF( IEQN )/d beta       C
C                                                                    C
C     INTARR    : INTARR( 1 ) contains L. Array of dim ( * )         C
C                                                                    C
C____________________________________________________________________C
C Output :-                                                          C
C ======                                                             C
C  Double Precision                                                  C
C  ----------------                                                  C
C     PMFSEF    : Value of function or partial derivative.           C
C____________________________________________________________________C
C
C*********************************************************************
      FUNCTION PMFSEF( NEQN, IEQN, ID, XVEC, INTARR, DPRARR )
      IMPLICIT NONE
C____________________________________________________________________C
C Variable declarations - Parameters ................................C
      INTEGER NEQN, IEQN, ID, INTARR( * )
      DOUBLE PRECISION PMFSEF, XVEC( * ), DPRARR( * )
C____________________________________________________________________C
C Variable declarations - Working variables .........................C
      INTEGER L
      DOUBLE PRECISION FACINN, FACOUT, RI, RO, ALPHA, BETA
C____________________________________________________________________C
C START OF PROGRAM **************************************************C
C____________________________________________________________________C
C
C First, trap any bad values of NEQN, IEQN and ID.
C
      IF ( NEQN.NE.2 ) THEN
         PRINT *,' Function PMFSEF.'
         PRINT *,' NEQN = ', NEQN,' and must be 2.'
         PRINT *,' Program aborted.'
         STOP
      ENDIF
C
      IF ( IEQN.NE.1 .AND. IEQN.NE.2 ) THEN
         PRINT *,' Function PMFSEF.'
         PRINT *,' IEQN = ', IEQN,' and must be 1 or 2.'
         PRINT *,' Program aborted.'
         STOP
      ENDIF
C
      IF ( ID.NE.0 .AND. ID.NE.1 .AND. ID.NE.2 ) THEN
         PRINT *,' Function PMFSEF.'
         PRINT *,' ID = ', ID,' and must be 0, 1 or 2.'
         PRINT *,' Program aborted.'
         STOP
      ENDIF
C
C ok - our command is valid
C
      L      = INTARR( 1 )
      RI     = DPRARR( 1 )
      RO     = DPRARR( 2 )
      ALPHA  = XVEC( 1 )
      BETA   = XVEC( 2 )
C
      FACINN = ALPHA*RI - BETA
      FACOUT = ALPHA*RO - BETA
C
C First do the functions
C
      IF ( ID.EQ.0 ) THEN
        IF ( IEQN.EQ.1 ) THEN
          PMFSEF = (-1.0d0)*RI*ALPHA*SIN( FACINN ) -
     1             DBLE( L )*COS( FACINN )
        ENDIF
        IF ( IEQN.EQ.2 ) THEN
          PMFSEF = (-1.0d0)*RO*ALPHA*SIN( FACOUT ) +
     1             DBLE( L + 1 )*COS( FACOUT )
        ENDIF
      ENDIF
C
C Now do the derivatives w.r.t. alpha
C
      IF ( ID.EQ.1 ) THEN
        IF ( IEQN.EQ.1 ) THEN
          PMFSEF = (-1.0d0)*RI*RI*ALPHA*COS( FACINN ) -
     1             RI*SIN( FACINN ) +
     2             DBLE( L )*RI*SIN( FACINN )
        ENDIF
        IF ( IEQN.EQ.2 ) THEN
          PMFSEF = (-1.0d0)*RO*RO*ALPHA*COS( FACOUT ) -
     1             RO*SIN( FACOUT ) -
     2             DBLE( L + 1 )*RO*SIN( FACOUT )
        ENDIF
      ENDIF
C
C Now do the derivatives w.r.t. beta
C
      IF ( ID.EQ.2 ) THEN
        IF ( IEQN.EQ.1 ) THEN
          PMFSEF = RI*ALPHA*COS( FACINN ) -
     1             DBLE( L )*SIN( FACINN )
        ENDIF
        IF ( IEQN.EQ.2 ) THEN
          PMFSEF = RO*ALPHA*COS( FACOUT ) +
     1             DBLE( L + 1 )*SIN( FACOUT )
        ENDIF
      ENDIF
C
      RETURN
      END
C*********************************************************************
