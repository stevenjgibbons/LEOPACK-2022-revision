C*********************************************************************
C subroutine Single Harmonic 2 QST coefficients **********************
C            -      -        - ---              **********************
C Steve Gibbons Fri Sep 24 17:08:43 BST 1999                         C
C____________________________________________________________________C
C                                                                    C
C  Returns a zero QST array except for a (1.0d0) in the element      C
C  specified by IL, IM and IT as indicated below.                    C
C____________________________________________________________________C
C Input variables :-                                                 C
C ===============                                                    C
C  Integer                                                           C
C  -------                                                           C
C     LH        : Maximum spherical harmonic degree, l.              C
C     IL        : Contains l of the single vector spherical harm.    C
C     IM        : Contains the order, m, of this vector spherical    C
C                 harmonic if it has a cos( m phi) dependence and    C
C                 contains (-m) if it has a sin( m phi) dependence.  C
C     IT        : Type of harmonic. IT = 1 --> scaloidal.            C
C                                   IT = 2 --> spheroidal.           C
C                                   IT = 3 --> toroidal.             C
C                                                                    C
C____________________________________________________________________C
C Output variables :-                                                C
C ================                                                   C
C  Double Precision                                                  C
C  ----------------                                                  C
C     QST       : Output array containing scaloidal/spheroidal and   C
C                 toroidal decomposition of vector. Is entirely zero C
C                 except for a 1.0d0 in the element indicated by     C
C                 IL, IM and IT.                                     C
C                                                                    C
C                 IF IM.GE.0 then ICS = 1 and                        C
C                 if IM.LT.0 then ICS = 2                            C
C                                                                    C
C                 M = IABS( IM )                                     C
C                                                                    C
C                 IND = INDSHC( IL, M, ICS )                         C
C                                                                    C
C                 and QST( IND, IT ) is then 1.0d0 and all other     C
C                 elements are zero.                                 C
C     ZCOEF     : Coefficient of Q_0^0                               C
C                                                                    C
C____________________________________________________________________C
C
C*********************************************************************
      SUBROUTINE SH2QST ( LH, IL, IM, IT, QST, ZCOEF )
      IMPLICIT NONE
C____________________________________________________________________C
C Variable declarations - Parameters ................................C
      INTEGER LH, IL, IM, IT
      DOUBLE PRECISION QST(  LH*(LH+2) , 3), ZCOEF
C____________________________________________________________________C
C Variable declarations - Working variables .........................C
      INTEGER M, ICS, INDSHC, IND
      DOUBLE PRECISION DZERO
      PARAMETER ( DZERO = 0.0d0 )
C____________________________________________________________________C
C START OF PROGRAM **************************************************C
C____________________________________________________________________C
C
      IF ( IL.GT.LH ) THEN
        PRINT *,' Subroutine SH2QST. IL = ',IL
        PRINT *,' LH = ', LH,' Program aborted.'
        STOP
      ENDIF
C
      IF ( IT.NE.1 .AND. IT.NE.2 .AND. IT.NE.3 ) THEN
        PRINT *,' Subroutine SH2QST. IT = ',IT
        PRINT *,' Program aborted.'
        STOP
      ENDIF
C
C Zero QST array
C
      ICS = 0
      CALL MATOP( QST, DZERO, LH*(LH+2), 3, ICS )
      ZCOEF = 0.0d0
C
C Find non-zero element
C
      M = IABS( IM )
      IF ( IM.GE.0 ) ICS = 1
      IF ( IM.LT.0 ) ICS = 2
C
      IF ( IL.EQ.0 .AND. IM.EQ.0 .AND. IT.EQ.1 ) THEN
        ZCOEF = 1.0d0
      ELSE
        IND = INDSHC( IL, M, ICS )
        QST( IND, IT ) = 1.0d0
      ENDIF
C
      RETURN
      END
C*********************************************************************
