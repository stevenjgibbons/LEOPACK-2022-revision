C*********************************************************************
C subroutine Radial QST array DisPlay ********************************
C            -      ---       -  -    ********************************
C Steve Gibbons Wed Dec  1 13:15:23 GMT 1999                         C
C____________________________________________________________________C
C                                                                    C
C Opens a file with name FNAME and logical unit number LU and writes C
C out the contents of the array RQST. This is a VERY verbose action  C
C and should only be used at the testing stage of programs.          C
C____________________________________________________________________C
C                                                                    C
C Input variables :-                                                 C
C ===============                                                    C
C  Integer                                                           C
C  -------                                                           C
C     LH        : Maximum spherical harmonic degree, l.              C
C     NR        : Number of radial grid nodes.                       C
C     LU        : Logical unit number of file.                       C
C     IF        : Format flag. Current options:-                     C
C       if = 1 --> 1PD11.4                                           C
C                                                                    C
C  Double Precision                                                  C
C  ----------------                                                  C
C                                                                    C
C     RQST      : Output array containing scaloidal/spheroidal       C
C                  decomposition of vector.                          C
C                  Has dimensions (  LH*(LH+2) ,3, NR ).             C
C              RQST (l*l+2m,1,I) = q_l^ms(r_i)                       C
C   RQST (l*l+2m-1+delta_m0,1,I) = q_l^mc(r_i)                       C
C              RQST (l*l+2m,2,I) = s_l^ms(r_i)                       C
C   RQST (l*l+2m-1+delta_m0,2,I) = s_l^mc(r_i)                       C
C              RQST (l*l+2m,3,I) = t_l^ms(r_i)                       C
C   RQST (l*l+2m-1+delta_m0,3,I) = t_l^mc(r_i)                       C
C                                                                    C
C  Character                                                         C
C  ---------                                                         C
C                                                                    C
C     FNAME     : Filename for output.                               C
C                                                                    C
C____________________________________________________________________C
C
C*********************************************************************
      SUBROUTINE RQSTDP( RQST, LH, NR, IF, LU, FNAME )
      IMPLICIT NONE
C____________________________________________________________________C
C Variable declarations - Parameters ................................C
      INTEGER LH, NR, IF, LU
      DOUBLE PRECISION RQST( LH*(LH+2), 3, NR )
      CHARACTER *(*) FNAME
C____________________________________________________________________C
C Variable declarations - Working variables .........................C
C
      INTEGER IR, IH, L, M, ICS, IWR
      DOUBLE PRECISION Q, S, T, TOTAL, LOW, QTOT, STOT, TTOT
      PARAMETER ( LOW = 1.0d-7 )
      LOGICAL OQ, OS, OT
C____________________________________________________________________C
C START OF PROGRAM **************************************************C
C____________________________________________________________________C
C
C Check input parameters ...
C
      IF ( IF.EQ.1 ) GOTO 50
      PRINT *,' Subroutine RQSTDP.'
      PRINT *,' IF = ', IF,'. Illegal value. Program aborted.'
      STOP
C
 50   CONTINUE
C
C Open file
C
      IWR = 3
      CALL FOPEN ( LU, FNAME, IWR )
C
      DO IH = 1, LH*(LH+2)
        CALL LMFIND( IH, L, M, ICS )
        TOTAL = 0.0d0
        QTOT  = 0.0d0
        STOT  = 0.0d0
        TTOT  = 0.0d0
        OQ    = .TRUE.
        OS    = .TRUE.
        OT    = .TRUE.
        DO IR = 1, NR
          Q = RQST( IH, 1, IR )
          S = RQST( IH, 2, IR )
          T = RQST( IH, 3, IR )
          TOTAL = TOTAL + Q*Q + S*S + T*T
          QTOT  = QTOT  + Q*Q              
          STOT  = STOT        + S*S      
          TTOT  = TTOT              + T*T
        ENDDO
        IF ( TOTAL.LT.LOW ) THEN
          WRITE( LU, 100 ) L, M, ICS
        ELSE
          IF ( QTOT.LT.LOW ) OQ = .FALSE.
          IF ( STOT.LT.LOW ) OS = .FALSE.
          IF ( TTOT.LT.LOW ) OT = .FALSE.
          DO IR = 1, NR
            Q = RQST( IH, 1, IR )
            S = RQST( IH, 2, IR )
            T = RQST( IH, 3, IR )
            IF ( IF.EQ.1 .AND. OQ .AND. OS .AND. OT )
     1              WRITE ( LU, 101 ) L, M, ICS, IR, Q, S, T
            IF ( IF.EQ.1 .AND. OQ .AND. (.NOT. OS) .AND. (.NOT. OT ) )
     1              WRITE ( LU, 102 ) L, M, ICS, IR, Q
            IF ( IF.EQ.1 .AND. OS .AND. (.NOT. OQ) .AND. (.NOT. OT ) )
     1              WRITE ( LU, 103 ) L, M, ICS, IR, S
            IF ( IF.EQ.1 .AND. OT .AND. (.NOT. OQ) .AND. (.NOT. OS ) )
     1              WRITE ( LU, 104 ) L, M, ICS, IR, T
            IF ( IF.EQ.1 .AND. OQ .AND. OS .AND. (.NOT. OT) )
     1              WRITE ( LU, 105 ) L, M, ICS, IR, Q, S
            IF ( IF.EQ.1 .AND. OT .AND. OS .AND. (.NOT. OQ) )
     1              WRITE ( LU, 106 ) L, M, ICS, IR, S, T
            IF ( IF.EQ.1 .AND. OT .AND. OQ .AND. (.NOT. OS) )
     1              WRITE ( LU, 107 ) L, M, ICS, IR, Q, T
          ENDDO
        ENDIF
      ENDDO
C
      CALL FCLOSE ( LU, FNAME, 'Error closing file.' )
 100  FORMAT('l= ',I4,' m= ',I4,' c= ',I1,' Q, S and T all zero.')
 101  FORMAT('l= ',I4,' m= ',I4,' c= ',I1,' r= ',I4,' Q: ',
     1        1PD11.4,' S: ',1PD11.4,' T: ',1PD11.4 )
 102  FORMAT('l= ',I4,' m= ',I4,' c= ',I1,' r= ',I4,' Q: ',
     1        1PD11.4,' S:             T:            ' )
 103  FORMAT('l= ',I4,' m= ',I4,' c= ',I1,' r= ',I4,' Q:            ',
     1        ' S: ',1PD11.4,' T:            ' )
 104  FORMAT('l= ',I4,' m= ',I4,' c= ',I1,' r= ',I4,' Q:            ',
     1        ' S:             T: ',1PD11.4 )
 105  FORMAT('l= ',I4,' m= ',I4,' c= ',I1,' r= ',I4,' Q: ',
     1        1PD11.4,' S: ',1PD11.4,' T:            ' )
 106  FORMAT('l= ',I4,' m= ',I4,' c= ',I1,' r= ',I4,' Q:            ',
     1       ' S: ',1PD11.4,' T: ',1PD11.4 )
 107  FORMAT('l= ',I4,' m= ',I4,' c= ',I1,' r= ',I4,' Q: ',
     1        1PD11.4,' S:             T: ',1PD11.4 )
      RETURN
      END
C*********************************************************************
