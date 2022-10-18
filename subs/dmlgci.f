C*********************************************************************
C subroutine Distributed Matrix Local and Global Coordinates Inter.  *
C            -           -      -         -      -           -       *
C Steve Gibbons Fri Aug 17 10:49:06 MET DST 2001                     C
C____________________________________________________________________C
C                                                                    C
C  If ( IDIR .ge. 0 ) then DMLGCI will take the global coordinate    C
C IGL - which may be row or column (i or j). Given NBLOCK - the      C
C length of the block in this direction - and NPROC - the number of  C
C processes in this direction - DMLGCI will calculate IPCRD (the row C
C - or column if appropriate - coordinate of the process in the      C
C process grid which is responsible for this element, and ILO - the  C
C local coordinate.                                                  C
C                                                                    C
C  If ( IDIR .lt. 0 ) then DMLGCI will take IPCRD, NBLOCK, NPROC and C
C ILO and return IGL.                                                C
C                                                                    C
C      Example:                                                      C
C      --------                                                      C
C                                                                    C
C     MB     = row block-size                                        C
C     NB     = column block-size                                     C
C     NPROW  = number of processes in row of process grid            C
C     NPCOL  = number of processes in col of process grid            C
C                                                                    C
C     MYROW  = process row coordinate for element with global        C
C               coordinate  IGLOBAL                                  C
C                                                                    C
C     MYCOL  = process col coordinate for element with global        C
C               coordinate  JGLOBAL                                  C
C                                                                    C
C     IDIR   = 1     
C     CALL DMLGCI( IDIR, MB, NPROW, MYROW, IGLOBAL, ILO )
C     CALL DMLGCI( IDIR, NB, NPCOL, MYCOL, JGLOBAL, JLO )
C     PRINT *,' Element is dealt with by process with coordinates'
C     PRINT *,' (',MYROW,',',MYCOL,')'
C     PRINT *,' Local coordinates are (',ILO,',',JLO,')'
C     IDIR   = -1     
C     CALL DMLGCI( IDIR, MB, NPROW, MYROW, IGLOBAL, ILO )
C     CALL DMLGCI( IDIR, NB, NPCOL, MYCOL, JGLOBAL, JLO )
C     PRINT *,' Global coordinates are (',IGLOBAL,',',JGLOBAL,')'
C____________________________________________________________________C
C                                                                    C
C Input variables :-                                                 C
C ===============                                                    C
C                                                                    C
C  Integer                                                           C
C  -------                                                           C
C                                                                    C
C     IDIR      : Set .ge. 0 to do global --> local                  C
C                 Set .lt. 0 to do local --> global                  C
C                                                                    C
C     NBLOCK    : Block length in the given direction.               C
C                                                                    C
C     NPROC     : Number of processes in the given direction.        C
C                                                                    C
C     IPCRD     : Process coordinate for the current element.        C
C                                                                    C
C     IGL       : Global coordinate of element.                      C
C                                                                    C
C     ILO       : Local coordinate of element.                       C
C                                                                    C
C____________________________________________________________________C
C
C*********************************************************************
      SUBROUTINE DMLGCI( IDIR, NBLOCK, NPROC, IPCRD, IGL, ILO )
      IMPLICIT NONE
C____________________________________________________________________C
C Variable declarations - Parameters ................................C
      INTEGER IDIR, NBLOCK, NPROC, IPCRD, IGL, ILO
C____________________________________________________________________C
C Variable declarations - Working variables .........................C
      INTEGER IM1, IM1BNB
C____________________________________________________________________C
C START OF PROGRAM **************************************************C
C____________________________________________________________________C
C
      IF ( IDIR.GE.0 ) THEN
C
C First do the global --> local
C
        IM1    = IGL - 1
        IM1BNB = IM1/NBLOCK
        IPCRD  = MOD( IM1BNB, NPROC )
        ILO    = IM1BNB/NPROC*NBLOCK + MOD( IM1, NBLOCK ) + 1
C
      ELSE
C
C Now do the local --> global
C
        IM1    = ILO - 1
        IM1BNB = IM1/NBLOCK
        IGL    = IPCRD*NBLOCK + IM1BNB*NBLOCK*NPROC + 
     1            MOD( IM1, NBLOCK ) + 1
C
      ENDIF
C
      RETURN
      END
C*********************************************************************
