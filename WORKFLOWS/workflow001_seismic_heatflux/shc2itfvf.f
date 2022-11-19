C**************************************************************
C Steve Gibbons   Spherical Harmonic Coef. 2 ITFVF input file *
C                 -         -        -     - -----            *
C Mon Jun 18 12:23:20 BST 2001                                C
C As the name might suggest; this converts a file  of         C
C spherical harmonic coefficients to an input file to be read C
C by itfvf.exe                                                C
C**************************************************************
	program shc2itfvf
        IMPLICIT NONE
C-------------------------------------------------------------
C Variable list
C-------------------------------------------------------------
      INTEGER LU, LH, LHMAX, L, M, ICS, IHARM, IWR
      PARAMETER ( LHMAX=20 )
      DOUBLE PRECISION SHC( LHMAX*( LHMAX + 2 ) )
      CHARACTER *(80) FNAME1, FNAME2
C-------------------------------------------------------------
C Beginning of program. ********************* START **********
C-------------------------------------------------------------
c
C First let's read in the spherical harmonic coefficients
      FNAME2 = ' File containing spherical harmonic coeff.s '
      CALL FNAMER ( FNAME1, FNAME2 )
c
      LU = 87
      CALL SHSFRD ( LH, LHMAX, LU, FNAME1, SHC, 1 )
c
C  Now let's open the output file
      FNAME1 = ' Name for itfvf input file.'
      CALL FNAMER ( FNAME2, FNAME1 )
      IWR = 3
      CALL FOPEN ( LU, FNAME2, IWR )
      DO IHARM = 1, LH*(LH+2)
        CALL LMFIND( IHARM, L, M, ICS )
        WRITE ( LU, 80 ) 'OB', L, M, ICS, SHC( IHARM )
      ENDDO
      CALL FCLOSE ( LU, FNAME2, FNAME1 )
 80   FORMAT(A2,' ',I5,' ',I5,' ',I2,' ',1pd16.7)
      STOP
      END

C*********************************************************************
C subroutine Spherical Harmonic Scalar File ReaD *********************
C            -         -        -      -    -  - *********************
C Steve Gibbons 13.1.98                                              C
C____________________________________________________________________C
C Reads in an 'nview' format file containing spherical harmonic      C
C coefficients.                                                      C
C____________________________________________________________________C
C Input Variables :-                                                 C
C ===============                                                    C
C  Integer                                                           C
C  -------                                                           C
C     LHMAX     : Maximum degree, l, of spherical harmonics.         C
C     LU        : Logical unit number.                               C
C     IFORM     : Format flag for data in file.                      C
C                                                                    C
C     Currently, only option is IFORM = 1 for a format (5d16.7)      C
C                                                                    C
C  Character                                                         C
C  ---------                                                         C
C     FNAME     : Filename - length not specified                    C
C____________________________________________________________________C
C Output Variables :-                                                C
C ================                                                   C
C  Integer                                                           C
C  -------                                                           C
C     LH        : Maximum degree, l, of spherical harmonics          C
C                  in file. Error is returned if this exceeds LHMAX. C
C                                                                    C
C  Double Precision                                                  C
C  ----------------                                                  C
C     SHC       : Spherical harmonic coefficients.                   C
C____________________________________________________________________C
C
C*********************************************************************
      SUBROUTINE SHSFRD ( LH, LHMAX, LU, FNAME, SHC, IFORM )
      IMPLICIT NONE
C____________________________________________________________________C
C Variable declarations - Parameters ................................C
      INTEGER LU, LH, LHMAX, IFORM
      DOUBLE PRECISION SHC( * )
      CHARACTER *( * ) FNAME
C____________________________________________________________________C
C Variable Declarations - Working Variables .........................C
      CHARACTER *( 80) LINE
      INTEGER IWR, NH, I
C____________________________________________________________________C
C START OF PROGRAM **************************************************C
C____________________________________________________________________C
C
      IF ( IFORM.NE.1 ) THEN
         PRINT *,' Subroutine SHSFRD.'
         PRINT *,' IFORM = ', IFORM
         PRINT *,' Not a recognised option.'
         STOP
      ENDIF
C
      IWR = 1
C open file ...
      CALL FOPEN ( LU, FNAME, IWR )
C ignore first line ...
      READ ( LU, 88 ) LINE
C read second line ind read LH ...
      READ ( LU, 88 ) LINE
      READ ( LINE, * ) LH
      IF ( LH.GT.LHMAX ) THEN
         PRINT *,' Subroutine SHSFRD. LHMAX = ', LHMAX
         PRINT *,' Program aborted.'
         STOP
      ENDIF
      NH = LH * ( LH + 2 )
      IF ( IFORM.EQ.1 ) READ ( LU, 89 ) ( SHC( I ), I = 1, NH )
      CALL FCLOSE ( LU, FNAME, FNAME )

 88   FORMAT (A80)
 89   FORMAT (5d16.7)
      RETURN
      END
C*********************************************************************
C*********************************************************************
C subroutine File NAME giveR *****************************************
C            -    ----     - *****************************************
C Steve Gibbons 14.4.97 (Adapted from Dan Gordon's course.)          C
C____________________________________________________________________C
C Asks the user for a file name FNAME. Pretty simple really ...      C
C____________________________________________________________________C
C Input Variable :-						     C
C ==============   						     C
C  Character							     C
C  ---------							     C
C     LABEL	: Message arbitrary length  			     C
C____________________________________________________________________C
C Output Variable :-						     C
C ===============   						     C
C  Character							     C
C  ---------							     C
C     FNAME	: Filename arbitrary length			     C
C____________________________________________________________________C
C
C*********************************************************************
      SUBROUTINE FNAMER ( FNAME, LABEL )
      IMPLICIT NONE
      CHARACTER *(*) FNAME
      CHARACTER *(*) LABEL

      PRINT *,' Please enter a filename for FNAME.'
      PRINT *, LABEL
      READ (5, 200) FNAME
 200  FORMAT (A)
      RETURN
      END
C*********************************************************************
C********************************************************************
C SUBROUTINE File OPEN **********************************************
C            -    ---- **********************************************
C Steve Gibbons 14.4.97 (Adapted from Dan Gordon's Code)            C
C Routine modified 15.3.99
C___________________________________________________________________C
C Opens a file with number LU, name FNAME, access OACCES  and       C
C a flag IRW to indicate whether the file is to be read or written  C
C to. ( IRW=1 ==> read only, IRW=2 ==> write but only if the file   C
C doesn't already exist, IRW=3 ==> write regardless of whether file C
C exists or not.)						    C
C___________________________________________________________________C
C Input Variables :-						    C
C ===============   						    C
C  Integer							    C
C  -------							    C
C     LU	: File number					    C
C     IRW	: Read / Write Flag 				    C
C                  = 1 for read only		                    C
C                  = 2 for write (provided that file doesn't exist. C
C                  = 3 for write (regardless of existence of file.  C
C                  = 4 for append status.                           C
C  Character							    C
C  ---------							    C
C     FNAME	: File name					    C
C___________________________________________________________________C
C Working Variables :-						    C
C =================   						    C
C  Character							    C
C  ---------							    C
C     OACCES 	: Access flag - should be set to 'OLD' for read     C
C                          and 'UNKNOWN' for write                  C
C     CONTYN    : For a yes/no to IWR = 2 option.		    C
C     LABEL     : Null string to pass into FNAMER option.	    C
C  Logical							    C
C  -------							    C
C     LEXIST	: Existence of file. File present <==> LEXIST=.TRUE.C
C___________________________________________________________________C
C Subroutines Used :-                                               C
C ================                                                  C
C     FNAMER	: For the case of IRW = 2; trying to write to an    C
C		   existing file. Used if alternative filename is   C
C                   asked for.					    C
C___________________________________________________________________C
C
C********************************************************************
      SUBROUTINE FOPEN ( LU, FNAME, IRW)
      IMPLICIT NONE
C___________________________________________________________________C
C Variable declarations - Parameters ...............................C
      INTEGER LU, IRW
      CHARACTER *(*) FNAME
C___________________________________________________________________C
C Variable declarations - Working Variables ........................C
      LOGICAL LEXIST
      CHARACTER *(7) OACCES
      CHARACTER *(1) CONTYN
      CHARACTER *(1) LABEL
C___________________________________________________________________C
      IF ( LU.EQ.0 ) THEN
         PRINT *,' Subroutine FOPEN'
         PRINT *,' I bet you ve forgotten to set LU ...??'
         PRINT *,' Think again and come back when you have'
         PRINT *,' remembered that LU must be a non zero integer!'
         PRINT *,' See you later. Bye for now!!'
         STOP
      ENDIF
C************************
C temporary code : SJG Thu Jun  1 08:07:06 BST 2000
C The Linux compiler will not allow file opening with
C lu.ge.100, so the following line should prevent it.
C 
      IF ( LU.GT.99 ) THEN
         PRINT *,' Subroutine FOPEN'
         PRINT *,' LU = ', LU,' too large.'
         PRINT *,' Program aborted.'
         STOP
      ENDIF
C************************
C--------------
 600  CONTINUE
      INQUIRE (FILE=FNAME, EXIST=LEXIST)
C Case of read only - LEXIST must be = .TRUE.
      IF ( IRW.EQ.1 ) THEN
         OACCES = 'OLD'
         IF ( .NOT. LEXIST ) THEN
            PRINT *,' Subroutine FOPEN. You are trying to open an'
            PRINT *,' old file which does not exist.'
            PRINT *,' Filename = ',FNAME,' LU= ',LU
            PRINT *,' Program aborted.'
            STOP
         ELSE
            GOTO 500
         ENDIF
      ENDIF
C Case of write to file provided that it doesn't exist
      IF ( IRW.EQ.2 ) THEN
         OACCES = 'UNKNOWN'
         IF ( LEXIST ) THEN
            PRINT *, ' Subroutine FOPEN. You are trying to write'
            PRINT *, ' to an existing file with IRW set to 2.'
            PRINT *,' Filename = ',FNAME,' LU= ',LU
            PRINT *,' Do you wish to give an alternative FNAME?'
            PRINT *,' Type y or n.'
            READ ( 5, 267) CONTYN
 267         FORMAT (A)
            IF (CONTYN.NE.'y'.AND.CONTYN.NE.'Y') THEN
               PRINT *, ' Program Aborted.'
               STOP
            ELSE
               LABEL=' '
               CALL FNAMER ( FNAME, LABEL )
               GOTO 600
            ENDIF
         ELSE
            GOTO 500
         ENDIF
      ENDIF
C Case of write to file regardless of the existence of file.
      IF ( IRW.EQ.3 ) THEN
         OACCES = 'UNKNOWN'
         GOTO 500
      ENDIF
C Treat appendment case
      IF ( IRW.EQ.4 ) THEN
         OACCES = 'UNKNOWN'
         OPEN ( UNIT=LU , FILE=FNAME , STATUS=OACCES,
     1          ACCESS='APPEND', ERR=999 )
         RETURN
      ENDIF
C___________________________________________________________________C
C All the IRW cases as of 14.4.97 have now been covered.
      PRINT *,' Subroutine FOPEN. IRW must be set to 1, 2, 3 or 4.'
      PRINT *,' Program aborted.'
      STOP

 500  CONTINUE
      OPEN ( UNIT=LU , FILE=FNAME , STATUS=OACCES, ERR=999 )
      RETURN

 999  PRINT *,' Subroutine FOPEN. Error in opening file ',FNAME
      STOP

      END
C********************************************************************
C___________________________________________________________________C
C********************************************************************
C subroutine File CLOSE *********************************************
C            -    ----- *********************************************
C Steve Gibbons 14.4.97                                             C
C  ( note that this is essentially the routine of Dan Gordon        C
C___________________________________________________________________C
C Closes file with integer logical unit LU, filename FNAME.         C
C LABEL contains any other information regarding the nature of the  C
C the file.							    C
C___________________________________________________________________C
C Input Variables :-						    C
C ===============						    C
C  Integer							    C
C  -------							    C
C     LU	: Number of file.				    C
C  Character							    C
C  ---------							    C
C     FNAME	: Name of file. Undefined length		    C
C     LABEL	: Any further information. Undefined length         C
C___________________________________________________________________C
C
C********************************************************************
      SUBROUTINE FCLOSE ( LU, FNAME, LABEL )
      IMPLICIT NONE
C___________________________________________________________________C
C Variable Declarations - Parameters ...............................C
      INTEGER LU
      CHARACTER *(*) FNAME
      CHARACTER *(*) LABEL
C___________________________________________________________________C
C START OF PROGRAM *************************************************C
C___________________________________________________________________C
      IF ( LU.EQ.0 ) THEN
         PRINT *,' Subroutine FCLOSE '
         PRINT *,' I bet you ve forgotten to set LU ...??'
         PRINT *,' Think again and come back when you have'
         PRINT *,' remembered that LU must be a non zero integer!'
         PRINT *,' See you later. Bye for now!!'
         STOP
      ENDIF
C----------------------------
      CLOSE (UNIT=LU, STATUS='KEEP', ERR=989 )
      RETURN
C
 989  PRINT *,' Error.  Failed to close ', LABEL, ' file ', FNAME
      STOP
      END
C********************************************************************
C********************************************************************
C subroutine L and M FIND *******************************************
C            -     - ---- *******************************************
C Steve Gibbons 12.4.97                                             C
C Modified 12.6.97 to have IT = 1 for cos and 2 for sine            C
C___________________________________________________________________C
C Given a number of harmonic, N; LMFIND will return the level, L,   C
C the order, M, and IT - which is equal to 2 for sin                C
C spherical harmonics and 1 for cosine ones.                        C
C All the above are integers - no point in a variable list .....    C
C     N = L*L for M = 0, IT = 1                                     C
C     N = L*L + 2*M - 1 for non-zero M and IT = 1                   C
C     N = L*L + 2*M for non-zero M and IT = 2                       C
C If N = 0, we have a monopole: L = 0, M = 0 and IT = 1.            C
C___________________________________________________________________C
C
C********************************************************************
      SUBROUTINE LMFIND ( N, L, M, IT)
      IMPLICIT NONE
C___________________________________________________________________C
C Variable Declarations - Parameters ...............................C
      INTEGER N,L,M,IT
C___________________________________________________________________C
C Variable declarations - Working Variables ........................C
      INTEGER LL,IDIFF,N2,ITWO
      PARAMETER (ITWO=2)
C___________________________________________________________________C
C First put N into N2 so that N is not altered
      N2=N
      IF ( N2.EQ.0 ) THEN
         L=0
         M=0
         IT=1
         RETURN
      ENDIF
      IF ( N2.LT.1 ) THEN 
         WRITE (6,989)
         STOP
 989  FORMAT (' Subroutine LMFIND. N<1 - Program stopped.')
      ENDIF
      IF ( N2.EQ.1 ) THEN
         L=1
         M=0
         IT=1
         RETURN
      ENDIF
      L=1
 500  CONTINUE
      LL=L*L
      IDIFF = N2 - LL
      IF ( IDIFF.GT.0 ) THEN
         L=L+1
         GOTO 500
      ENDIF
      IF ( IDIFF.EQ.0 ) THEN
         M=0
         IT=1
         RETURN
      ENDIF
      L=L-1
      LL=L*L
      N2=N2-LL
C so we know now that N2 is equal to either 2*M or 2*M-1
C corresponding to IT=1 and IT=2 respectively
      IDIFF = MOD ( N2, ITWO )
      IF ( IDIFF.EQ.1) THEN
         IT = 1
         M = (N2+1)/ITWO
      ELSE
         IT = 2
         M = N2/ITWO
      ENDIF
      RETURN
      END
C********************************************************************
