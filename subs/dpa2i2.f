C*********************************************************************
C subroutine Double Precision Array 2 Idl format 2 *******************
C            -      -         -     - -          - *******************
C Steve Gibbons Thu May 11 13:58:39 BST 2000                         C
C____________________________________________________________________C
C                                                                    C
C Given a double precision array DARR( N1, N2 ), DPA2I2 will write   C
C to a file with logical unit LU a crude input for an IDL script     C
C containing the data.                                               C
C                                                                    C
C dpa2if appears to put the indices in the wrong order, so the       C
C purpose of DPA2IF is to correct this.                              C
C____________________________________________________________________C
C                                                                    C
C Input variables :-                                                 C
C ===============                                                    C
C  Integer                                                           C
C  -------                                                           C
C                                                                    C
C     N1        : First dimension of array DARR.                     C
C     N2        : Second dimension of array DARR.                    C
C     LU        : Logical unit number of output file.                C
C     NCH       : Number of characters in variable name.             C
C                                                                    C
C  Character                                                         C
C  ---------                                                         C
C                                                                    C
C     VARNAM    : Variable name (*)                                  C
C                                                                    C
C  Double Precision                                                  C
C  ----------------                                                  C
C                                                                    C
C     DARR      : Data array with bounds ( N1, N2 )                  C
C____________________________________________________________________C
C
C*********************************************************************
      SUBROUTINE DPA2I2( N1, N2, LU, NCH, VARNAM, DARR )
      IMPLICIT NONE
C____________________________________________________________________C
C Variable declarations - Parameters ................................C
      INTEGER N1, N2, LU, NCH
      CHARACTER *(*) VARNAM
      DOUBLE PRECISION DARR( N1, N2 )
C____________________________________________________________________C
C Variable declarations - Working variables .........................C
      INTEGER NER, I1, I2, NROWS, IREM, IROW, ISR, I
      CHARACTER *(80) LINE
      CHARACTER *(3)  CHEL
      DOUBLE PRECISION EL( 4 )
C____________________________________________________________________C
C START OF PROGRAM **************************************************C
C____________________________________________________________________C
C
      LINE(1:NCH)       = VARNAM(1:NCH)
      LINE(NCH+1:NCH+7) = ' = [ $ '
      WRITE( LU, 80 ) LINE(1:NCH+7)
C
      NROWS = N1/4
      IREM   = MOD( N1, 4 )
      IF ( IREM.NE.0 ) NROWS = NROWS + 1
      IF ( IREM.EQ.0 ) IREM = 4
C
      DO I2 = 1, N2
        ISR = 0
        IF ( I2.EQ.N2 ) THEN
          CHEL = ' ] '
        ELSE
          CHEL = ', $'
        ENDIF
        DO IROW = 1, NROWS
          IF ( IROW.LT.NROWS ) NER = 4
          IF ( IROW.EQ.NROWS ) NER = IREM
          DO I = 1, NER
            I1 = ISR + I
            EL( I ) = DARR( I1, I2 )
          ENDDO
          IF ( IROW.GT.1 .AND. IROW.LT.NROWS ) THEN
            WRITE ( LU, 82 ) EL( 1 ), EL( 2 ), EL( 3 ), EL( 4 )
          ENDIF
          IF ( IROW.EQ.1 .AND. NROWS.GT.1 ) THEN
            WRITE ( LU, 81 ) EL( 1 ), EL( 2 ), EL( 3 ), EL( 4 )
          ENDIF
          IF ( NROWS.EQ.1 ) THEN
            IF ( NER.EQ.4 ) WRITE ( LU, 87 ) EL( 1 ), EL( 2 ),
     1                     EL( 3 ), EL( 4 ), CHEL
            IF ( NER.EQ.3 ) WRITE ( LU, 88 ) EL( 1 ), EL( 2 ),
     1                     EL( 3 ), CHEL
            IF ( NER.EQ.2 ) WRITE ( LU, 89 ) EL( 1 ), EL( 2 ),
     1                     CHEL
            IF ( NER.EQ.1 ) WRITE ( LU, 90 ) EL( 1 ), CHEL
          ENDIF
          IF ( IROW.EQ.NROWS .AND. NROWS.GT.1 ) THEN
            IF ( NER.EQ.4 ) WRITE ( LU, 83 ) EL( 1 ), EL( 2 ),
     1                     EL( 3 ), EL( 4 ), CHEL
            IF ( NER.EQ.3 ) WRITE ( LU, 84 ) EL( 1 ), EL( 2 ),
     1                     EL( 3 ), CHEL
            IF ( NER.EQ.2 ) WRITE ( LU, 85 ) EL( 1 ), EL( 2 ),
     1                     CHEL
            IF ( NER.EQ.1 ) WRITE ( LU, 86 ) EL( 1 ), CHEL
          ENDIF
          ISR = ISR + 4
        ENDDO
 50   CONTINUE
      ENDDO
C
 80   FORMAT(A)
 81   FORMAT('[ ',1PD14.5,',',1PD14.5,',',1PD14.5,',',1PD14.5,', $ ')
 82   FORMAT('  ',1PD14.5,',',1PD14.5,',',1PD14.5,',',1PD14.5,', $ ')
 83   FORMAT('  ',1PD14.5,',',1PD14.5,',',1PD14.5,',',1PD14.5,' ] ',A3)
 84   FORMAT('  ',1PD14.5,',',1PD14.5,',',1PD14.5,' ] ',A3)
 85   FORMAT('  ',1PD14.5,',',1PD14.5,' ] ',A3)
 86   FORMAT('  ',1PD14.5,' ] ',A3)
 87   FORMAT(' [',1PD14.5,',',1PD14.5,',',1PD14.5,',',1PD14.5,' ]',A3)
 88   FORMAT(' [',1PD14.5,',',1PD14.5,',',1PD14.5,' ] ',A3)
 89   FORMAT(' [',1PD14.5,',',1PD14.5,' ] ',A3)
 90   FORMAT(' [',1PD14.5,' ] ',A3)
C
      RETURN
      END
C*********************************************************************
