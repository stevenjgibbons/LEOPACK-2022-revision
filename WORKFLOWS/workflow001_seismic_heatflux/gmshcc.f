C*********************************************************************
C Steve Gibbons 13.1.98 
C
C gmshcc - Guy Master's Spherical Harmonic Coefficient Convert
C
C Reads a file of complex spherical harmonic coefficients
C as formatted by T. G. Masters 
C and converts them into Schmidt Quasi Normalised harmonics
C ( \int_0^{\pi} ( P_l^m )^2 \sin \theta d \theta =
C \frac{ 4 \pi }{ 2l + 1 } as used by SG convection code.
C
C The file containing Guy's coefficients is (lmax*(lmax+1))/2 - 1
C lines long. lmax here is 10 and so our coefficient file is
C (10*11)/2 - 1 = 65 lines long.
C
C*********************************************************************
      PROGRAM gmshcc
      IMPLICIT NONE
C-----------------------------------------------------
      INTEGER LH, LHMAX, NPOLMX, LU, IWR, NPOLH,J1, J2,
     1        L, M, ICS, J, LU2, KLM, IND
      PARAMETER ( LHMAX = 10, NPOLMX = LHMAX * ( LHMAX + 2 )  )
      INTEGER LMIPOL( NPOLMX, 3 )
      DOUBLE PRECISION BHCOEF( NPOLMX ), A, B, 
     1                 COEFR( NPOLMX ), COEFI( NPOLMX )
      CHARACTER *(80) FNAME, TITLE, INFILE
C
      LU = 87
      LU2 = 89
C-------------------------------------------------------------
C This is the file name of the output file.
C
      FNAME = 'output_coeffs.txt'
C
C Change it as is appropriate ...
C-------------------------------------------------------------
      IWR = 3
      TITLE = 'Scripps shear wave velocity anomaly map'
C
c     PRINT *,' Enter LH =< 10 '
c     READ ( 5, * ) LH
c     IF ( LH.GT.10 ) THEN
c       PRINT *,' You have asked to truncated at degree ',LH
c       PRINT *,' Program aborted.'
c       STOP
c     ENDIF
      LH = 10
C-------------------------------------------------------------
C This is the file name of the input file.
C
      INFILE = 'lower_mantle_coefs.txt'
C
C Change it as is appropriate ...
C-------------------------------------------------------------
      CALL FOPEN ( LU2, INFILE, 1 )
C      
C Next section simply loops around harmonics and stores
C the indices (l, m, cos/sin ) into an array. This is
C merely a convention for my convection codes.
C
      NPOLH = LH*(LH+2)
      DO J = 1, NPOLH
         CALL LMFIND ( J, L, M, ICS )
         LMIPOL( J, 1 ) = L
         LMIPOL( J, 2 ) = M
         LMIPOL( J, 3 ) = ICS
      ENDDO
C
C We now read in Guy's coefficients and
C store the real parts in the array COEFR
C and the imaginary parts in the array COEFI.
C They have been rearranged into an order
C governed by the function KLM.
C This is another convention of mine
C to be compatible with my Schmidt normalised
C Associated Legendre Function routines.
C     
 500  CONTINUE
      READ ( LU2, *, END=501 ) J1, J2, A, B
      IND = KLM( J1, J2 )
      COEFR( IND ) = A
      COEFI( IND ) = B
      GOTO 500
 501  CONTINUE
c     .
C
C We now pass the arrays COEFR and COEFI into
C the routine MAS2SN where the coefficients
C will be converted into harmonics with
C Schmidt Quasi Normalisation and stored
C in an array BHCOEF.
C The coefficient in BHCOEF( i )
C has l given by LMIPOL( i, 1 )
C has m given by LMIPOL( i, 2 )
C and LMIPOL( i, 3 ) = 1 for a cos ( m phi ) harmonic
C and LMIPOL( i, 3 ) = 2 for a sin ( m phi ) harmonic
C
      CALL MAS2SN ( LH, COEFR, COEFI, BHCOEF)
c     .
      CALL FCLOSE ( LU2, INFILE , INFILE )
c     .
      CALL SHFWR ( LH, NPOLH, NPOLMX, LMIPOL, LU, IWR,
     1                   BHCOEF, FNAME, TITLE )
c     .
      DO J = 1, NPOLH
         CALL LMFIND ( J, L, M, ICS )
         LMIPOL( J, 1 ) = L
         LMIPOL( J, 2 ) = M
         LMIPOL( J, 3 ) = ICS
        WRITE ( 6, * ) L, M, ICS, BHCOEF( J )
      ENDDO
c     .
      STOP
      END
c     .
C*********************************************************************
C*********************************************************************
C subroutine MASters 2 Schmidt Normalisation *************************
C            ---     - -       -             *************************
C Steve Gibbons 31.7.98                                              C
C____________________________________________________________________C
C                                                                    C
C Converts the real and imaginary complex spherical harmonic         C
C coefficients COEFR and COEFI (from Guy Masters shear wave velocity C
C anomaly model); using the routine GXFCN and the function SHMPLG    C
C to calculate the differently normalised spherical harmonic coeff.s C
C                                                                    C
C____________________________________________________________________C
C                                                                    C
C Input variables :-                                                 C
C ===============                                                    C
C  Integer                                                           C
C  -------                                                           C
C     LH        : Maximum degree of spherical harmonic.              C
C                                                                    C
C  Double Precision                                                  C
C  ----------------                                                  C
C                                                                    C
C     COEFR     : Real part of complex spherical harmonic coeff.     C
C     COEFI     : Imag part of complex spherical harmonic coeff.     C
C                                                                    C
C  ( coefr and coefi have dimension ( LH + 1 )*( LH + 2 )/2 )        C
C  ( their index is controlled by the function KLM( L, M )           C
C                                                                    C
C Output variables :-                                                C
C ================                                                   C
C  Double Precision                                                  C
C  ----------------                                                  C
C                                                                    C
C     SHC       : Real Schmidt-Normalised Spherical Harmonic coeff.s C
C                                                                    C
C  ( shc has dimension ( LH*(LH + 2 ) )                              C
C  ( their position is defined by the function INDSHC )              C
C                                                                    C
C____________________________________________________________________C
C
C*********************************************************************
      SUBROUTINE MAS2SN ( LH, COEFR, COEFI, SHC )
      IMPLICIT NONE
C____________________________________________________________________C
C Variable declarations - Parameters ................................C
      INTEGER LH
      DOUBLE PRECISION COEFR(*), COEFI(*), SHC(*)
C____________________________________________________________________C
C Variable declarations - Working variables .........................C
      INTEGER INDSHC, KLM, L, M, ICS, LHMAX, INDGM, INDSG
      DOUBLE PRECISION SHMPLG, C, S, TOL
      PARAMETER ( LHMAX = 14, TOL = 1.0d-6 )
      DOUBLE PRECISION X( 0:LHMAX ), X1( 0:LHMAX ), X2( 0:LHMAX )
C____________________________________________________________________C
C START OF PROGRAM **************************************************C
C____________________________________________________________________C
C     .
C     . Check inputs
C     .
      IF ( LH.GT.LHMAX ) THEN
         PRINT *,' Subroutine MAS2SM. LH = ', LH
         PRINT *,' Too large. Program aborted.'
         STOP
      ENDIF
C     .
      C = 0.4123d0
      S = DSQRT( 1.0d0 - C*C )
C     .
      DO L = 1, LH
 500     CONTINUE
         CALL GXFCN ( L, C, S, X(0), X1(0), X2(0) )
         DO M = 0, L
           IF ( ABS(X( M )).LT.TOL ) THEN
             C = C + 0.005
             S = DSQRT( 1.0d0 - C*C )
             GOTO 500
           ENDIF
C
C Note - I am now overwriting X1(M) with
C a Schmidt Quasi Normalised Associated Legendre Function
C which is given by the function SHMPLG as a function
C of C = cos(theta).
C If you want to convert this program to directly
C deal with a S.H. normalisation you use, you
C must supply the appropriate function instead of
C SHMPLG in the next line. Steve Gibbons 1st Sept. 1999
C
           X1( M ) = SHMPLG( L, M, C )
         ENDDO
c        .
c        . ok let's translate the coefficients
c        .
         M = 0
         ICS = 1
         INDGM = KLM( L, M )
         INDSG = INDSHC( L, M, ICS )
         SHC( INDSG ) = COEFR( INDGM )*X( M )/X1( M )
c        .
         DO M = 1, L
           ICS = 1
           INDGM = KLM( L, M )
           INDSG = INDSHC( L, M, ICS )
           SHC( INDSG ) = COEFR( INDGM )*X( M )/X1( M )
c        .
           ICS = 2
           INDGM = KLM( L, M )
           INDSG = INDSHC( L, M, ICS )
           SHC( INDSG ) = (-1.0d0)*COEFI( INDGM )*X( M )/X1( M )
         ENDDO
c        .
      ENDDO
C     .     end do l = 1, lh
C     .
      RETURN
      END
C*********************************************************************
C*********************************************************************
C subroutine GXFCN ( don't know the acronym. It is a TGM one).       C
C            -----                                                   C
C Subroutine written by T.G.Masters to calculate associated          C
C Legendre Functions.   I'm not sure what normalisation we're        C
C dealing with here; but by accurately copying the code and          C
C converting to double precision etc. I will be able to use this     C
C routine to convert them to Schmidt Quasi-Normalised.               C
C steve gibbons 31.7.98                                              C
C____________________________________________________________________C
C                                                                    C
C Here are Guy Masters' notes ...                                    C
C                                                                    C
C    A TGM special.                                                  C
C    c=cos(theta), s=sin(theta),theta=colatitude,l=angular order,    C
C    x(1) contains m=0, x(2) contains m=1, x(k+1) contains m=k       C
C    m=azimuthal order 0.le.m.le.l .  x1=dx/dtheta                   C
C    x(l,m,theta)*exp(imphi) are a set of orthonormal                C
C    spherical harmonics.                                            C
C                                                                    C
C    modified to also calculate the 2nd derivatives of               C
C    associated Legendre fns.                                        C
C                                                                    C
C____________________________________________________________________C
C                                                                    C
C Input variables :-                                                 C
C ===============                                                    C
C  Integer                                                           C
C  -------                                                           C
C                                                                    C
C     L         : Level of spherical harmonics being caluclated.     C
C                                                                    C
C  Double Precision                                                  C
C  ----------------                                                  C
C                                                                    C
C     C         : Cosine of the angle.  ( in radians )               C
C     S         : Sine of the angle.    ( in radians )               C
C                                                                    C
C Output variables :-                                                C
C ================                                                   C
C                                                                    C
C     X         : Associated Legendre Functions.                     C
C     X1        : Derivatives of the above.                          C
C     X2        : Second derivatives                                 C
C                                                                    C
C  ( x, x1 and x2 must have dimension atleast L + 1 )                C
C                                                                    C
C (Note - send arrays declared X(0:L) etc. and in                    C
C calling sequence do CALL GXFCN ( L, C, S, X(0), X1(0), X2(0) )     C
C This makes things easier with the indexing ...)                    C
C____________________________________________________________________C
C
C*********************************************************************
      SUBROUTINE GXFCN ( L, C, S, X, X1, X2 )
      IMPLICIT NONE
C____________________________________________________________________C
C Variable declarations - Parameters ................................C
      INTEGER L
      DOUBLE PRECISION C, S,
     1                 X( * ), X1( * ), X2( * )
C____________________________________________________________________C
C Variable declarations - Working variables .........................C
      INTEGER I, M, MM, LP1, LL, MI, MP1
      DOUBLE PRECISION PI, PI4, PMM, PLL, F, CON, COT, FL2P1
C____________________________________________________________________C
C START OF PROGRAM **************************************************C
C____________________________________________________________________C
C
      PI=DATAN(1.d0)*4.d0
      PI4 = PI*4.d0
      LP1=L+1
      FL2P1=L+LP1
      CON=SQRT((FL2P1)/PI4)
      IF(ABS(S).LT.1.D-20) GOTO 300
      COT=C/S
c*** try m recursions ; first compute XLL
      F=1.d0
      DO 100 i=1,L
  100 F=F*(I+I-1.d0)/(I+I)
      X(LP1)=SQRT(FL2P1*F/PI4)*(-S)**L
C*** use L recursions IF XLL underflows
      IF(ABS(X(LP1)).LT.1.D-295) GOTO 200
      X1(LP1)=X(LP1)*L*COT
c*** use m recurrence starting from m=L
      DO 110 I=1,L
      M=LP1-I
      MP1=M+1
      F=SQRT(I*(FL2P1-I))
      X(M)=-(X1(MP1)+M*COT*X(MP1))/F
  110 X1(M)=(M-1)*X(M)*COT+X(MP1)*F
      GOTO 1000
c
c*** resort to l recursions ; first compute Pmm's
  200 X(1)=1.d0
      DO 210 I=1,L
  210 X(I+1)=X(I)*S*(I+I-1.d0)
c*** perform l recurrence for each m value
      DO 220 MM=1,L
      M=MM-1
      PMM=0.d0
      DO 220 LL=MM,L
      PLL=(C*(LL+LL-1.d0)*X(MM)-(LL+M-1.d0)*PMM)/(LL-M)
      PMM=X(MM)
  220 X(MM)=PLL
c*** normalise and use m recurrence to compute derivative
      X(1)=X(1)*CON
      X1(1)=0.d0
      IF(L.GT.0) X1(1)=-X(2)*CON
      DO 230 MM=2,LP1
      M=MM-1
      F=DSQRT(DBLE((LP1-M)*(L+M)))
      CON=-CON/F
      X(MM)=X(MM)*CON
  230 X1(MM)=-M*COT*X(MM)-F*X(M)
      GOTO 1000
c
c*** handle very small arguments
  300 DO 310 I=1,LP1
      X(I)=0.d0
  310 X1(I)=0.d0
      X(1)=CON
      IF(L.GT.0) X1(2)=-0.5d0*CON*DSQRT(DBLE(L*LP1))
      GOTO 1000
c
c now use Legendre's equation for 2nd derivative
c
 1000 DO 1001 MI=1,LP1+1
      M=MI-1
c     X2(MI)=( 2.*C*X1(MI) - ( DBLE(L)*DBLE(L+1)
c    .       - DBLE(M)*DBLE(M)/(1.-C*C))*X(MI)) /(1.-C*C)
      X2(MI)=(M*M/(S*S) -  DBLE(L)*DBLE(L+1))*X(MI) - COT*X1(MI)
c     WRITE(1,*) M,X2(MI)
 1001 CONTINUE
c
      RETURN
      END
C*********************************************************************
C********************************************************************
C SUBROUTINE File OPEN **********************************************
C            -    ---- **********************************************
C Steve Gibbons 14.4.97 (Adapted from Dan Gordon's Code)            C
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
C     LABLE     : Null string to pass into FNAMER option.	    C
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
      CHARACTER *(1) LABLE
C___________________________________________________________________C
      IF ( LU.EQ.0 ) THEN
         PRINT *,' Subroutine FOPEN'
         PRINT *,' I bet you ve forgotten to set LU ...??'
         PRINT *,' Think again and come back when you have'
         PRINT *,' remembered that LU must be a non zero integer!'
         PRINT *,' See you later. Bye for now!!'
         STOP
      ENDIF
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
               LABLE=' '
               CALL FNAMER ( FNAME, LABLE )
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
C___________________________________________________________________C
C All the IRW cases as of 14.4.97 have now been covered.
      PRINT *,' Subroutine FOPEN. IRW must be set to 1, 2 or 3.'
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
C     LABLE	: Message arbitrary length  			     C
C____________________________________________________________________C
C Output Variable :-						     C
C ===============   						     C
C  Character							     C
C  ---------							     C
C     FNAME	: Filename arbitrary length			     C
C____________________________________________________________________C
C
C*********************************************************************
      SUBROUTINE FNAMER ( FNAME, LABLE )
      IMPLICIT NONE
      CHARACTER *(*) FNAME
      CHARACTER *(*) LABLE

      PRINT *,' Please enter a filename for FNAME.'
      PRINT *, LABLE
      READ (5, 200) FNAME
 200  FORMAT (A)
      RETURN
      END
C*********************************************************************
C********************************************************************
C subroutine File CLOSE *********************************************
C            -    ----- *********************************************
C Steve Gibbons 14.4.97                                             C
C  ( note that this is essentially the routine of Dan Gordon        C
C___________________________________________________________________C
C Closes file with integer logical unit LU, filename FNAME.         C
C LABLE contains any other information regarding the nature of the  C
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
C     LABLE	: Any further information. Undefined length         C
C___________________________________________________________________C
C
C********************************************************************
      SUBROUTINE FCLOSE ( LU, FNAME, LABLE )
      IMPLICIT NONE
C___________________________________________________________________C
C Variable Declarations - Parameters ...............................C
      INTEGER LU
      CHARACTER *(*) FNAME
      CHARACTER *(*) LABLE
C___________________________________________________________________C
C START OF PROGRAM *************************************************C
C___________________________________________________________________C
      IF ( LU.EQ.0 ) THEN
         PRINT *,' Subroutine FCLOSE '
         PRINT *,' I bet you ve forgotten to set LU ...??'
         PRINT *,' Think again and come back when you have'
         PRINT *,' remembered that LU must be a non zero integer!'
         STOP
      ENDIF
C----------------------------
      CLOSE (UNIT=LU, STATUS='KEEP', ERR=989 )
      RETURN

 989  PRINT *,' Error.  Failed to close ', LABLE, ' file ', FNAME
      STOP
      END
C********************************************************************
C*********************************************************************
C function ScHMidt normalised (Polynomial) LeGendre function *********
C          - --                -           - -               *********
C Steve Gibbons 8.5.97						     C
C____________________________________________________________________C
C Adapted from SCHNLF to give a single Schmidt Normalised Assoc.     C
C Legendre Function as a func. of integers L and M.                  C
C____________________________________________________________________C
C Input variables :-                                                 C
C ===============                                                    C
C  Integer                                                           C
C  -------                                                           C
C     L		: Level of Assoc. Legendre Polynomial.               C
C     M		: Order of Assoc. Legendre Polynomial.               C
C                                                                    C
C  Double Precision                                                  C
C  ----------------                                                  C
C     COSTH	: Cosine of theta.				     C
C____________________________________________________________________C
C Functions called:-						     C
C ----------------                                                   C
C PMM ( M, S)                                                        C
C PMM1 ( M, X, PMM0 )                                                C
C PLM ( L, M, X, PLMIN1, PLMIN2 )                                    C
C____________________________________________________________________C
C
C*********************************************************************
      FUNCTION SHMPLG ( L, M, COSTH )
      IMPLICIT NONE
C____________________________________________________________________C
C Variable declarations - Parameters ................................C
      INTEGER L,M
      DOUBLE PRECISION SHMPLG, COSTH
C____________________________________________________________________C
C Variable declarations - Working variables .........................C
      DOUBLE PRECISION SINE,POLD,POLD1
      INTEGER I
C____________________________________________________________________C
C Variable declarations - Functions called ..........................C
      DOUBLE PRECISION PMM, PMM1, PLM
C
C____________________________________________________________________C
C START OF PROGRAM **************************************************C
C____________________________________________________________________C
C Check validity of arguments ....
      IF (L.LT.1) THEN
         PRINT *,' Function SHMPLG. L.LT.1 Program Stopped. '
         STOP
      ENDIF
      IF (M.LT.0 .OR. M.GT.L ) THEN
         PRINT *,' Function SHMPLG.'
         PRINT *,' M must be between 0 and L.'
         PRINT *,' Program aborted.'
         STOP
      ENDIF
      IF ( COSTH.LT.-1.0d0 .OR. COSTH.GT.1.0d0 ) THEN
         PRINT *,' Function SHMPLG.'
         PRINT *,' Cos(theta) out of range.'
         PRINT *,' Program aborted.'
         STOP
      ENDIF
C ....................................................................
      SINE = (1.0d0 - COSTH)*(1.0d0 + COSTH )
      SINE = SQRT( SINE )
      SHMPLG = PMM ( M, SINE )
      IF ( L.EQ.M ) RETURN
      POLD = SHMPLG
      SHMPLG = PMM1 ( M, COSTH, POLD )
      IF ( L.EQ.M+1 ) RETURN
C ............................ so L is greater than M + 1.
      DO I = M + 2, L
         POLD1 = POLD
         POLD = SHMPLG
         SHMPLG = PLM( I, M, COSTH, POLD, POLD1 )
      ENDDO

      RETURN
      END
C*********************************************************************
C*********************************************************************
C function PMM *******************************************************
C Steve Gibbons 16.4.97                                              C
C____________________________________________________________________C
C Gives the Schmidt Normalised Legendre Function P_m^m ( X )         C
C from eqn. 175 in my notes ie                                       C
C                                                                    C
C                   ( 2m - 1)!! (1- XX)^(m/2) * SQRT (2.0d0 )        C
C   P_m^m( X ) =  ---------------------------------------------      C
C                        SQRT ( (2m)! )                              C
C                                                                    C
C       for m non-zero and                                           C
C                                                                    C
C   P_0^0( X ) = 1.0d0                                               C
C                                                                    C
C N.B. The double factorial sign means the product of all ODD        C
C integers between 1 and ( 2m - 1 ).                                 C
C Best to use the form in eq 184 to calculate PMM.                   C
C____________________________________________________________________C
C Input Variables :-                                                 C
C ===============                                                    C
C  Integer                                                           C
C  -------                                                           C
C     M         : Well it's M isn't it!                              C
C  Double Precision                                                  C
C  ----------------                                                  C
C     SINE      : Sin (theta)                                        C
C____________________________________________________________________C
C
C*********************************************************************
      FUNCTION PMM ( M, S )
      IMPLICIT NONE
      DOUBLE PRECISION PMM
C____________________________________________________________________C
C Variable declarations - Parameters ................................C
      INTEGER M
      DOUBLE PRECISION S
C____________________________________________________________________C
C Variable declarations - Working variables .........................C
      INTEGER I
      DOUBLE PRECISION TWOI
C____________________________________________________________________C
C START OF FUNCTION *************************************************C
C____________________________________________________________________C
C
      IF ( M.LT.0 ) THEN
         PRINT *,' Function PMM. M < 0 error. Stopped.'
         STOP
      ENDIF
      IF ( M.EQ.0 ) THEN
         PMM = 1.0d0
         RETURN
      ENDIF
C ................................. so M is greater than 0
      PMM = SQRT ( 2.0d0 )
      DO I = 1, M
         TWOI = 2.0d0*FLOAT(I)
         PMM = PMM * S * SQRT ( (TWOI - 1.0d0 )/TWOI )
      ENDDO
      RETURN
      END
C____________________________________________________________________C
C*********************************************************************
C function PMM1 ******************************************************
C Steve Gibbons 16.4.97                                              C
C____________________________________________________________________C
C Evaluates the Schmidt Normalised Legendre Function P_(m+1)^m (X)   C
C according to equation 179 in my notes ; i.e.                       C
C                                                                    C
C    P_(m+1)^m (X) = SQRT( 2m+1 ).X.P^m_m(X)                         C
C                                                                    C
C____________________________________________________________________C
C Input Variables :-                                                 C
C ===============                                                    C
C  Integer                                                           C
C  -------                                                           C
C     M         : Well it's M isn't it!                              C
C  Double Precision                                                  C
C  ----------------                                                  C
C     X         : Cos(theta)                                         C
C     PMM0      : P_m^m(X) as evaluated by Function PMM.             C
C____________________________________________________________________C
C
C*********************************************************************
      FUNCTION PMM1 ( M, X, PMM0 )
      IMPLICIT NONE
      DOUBLE PRECISION PMM1
C____________________________________________________________________C
C Variable declarations - Parameters ................................C
      INTEGER M
      DOUBLE PRECISION X, PMM0
C____________________________________________________________________C
C Variable declarations - Working variables .........................C
      DOUBLE PRECISION RM
C____________________________________________________________________C
C START OF FUNCTION *************************************************C
C____________________________________________________________________C
C
      IF ( M.LT.0 ) THEN
         PRINT *,' Function PMM1. M < 0 error. Stopped.'
         STOP
      ENDIF
      RM = FLOAT ( M )
      PMM1 = X*PMM0*SQRT( 2.0d0*RM+1.0d0 )
      RETURN
      END
C____________________________________________________________________C
C*********************************************************************
C function PLM *******************************************************
C Steve Gibbons 16.4.97                                              C
C____________________________________________________________________C
C Calculates the Schmidt Normalised Legendre Function P_l^m (x)      C
C given P_(l-1)^m and P_(l-2)^m according to equation 183 in my notesC
C i.e.                                                               C
C   P_l^m( X ) = { A * P_(l-1)^m - B * P_(l-2)^m }/C                 C
C                                                                    C
C where A = (2*l - 1)*X ,                                            C
C                                                                    C
C B = SQRT( (L+M-1)*(L-M-1) ) and C = SQRT( (L+M)*(L-M) )            C
C                                                                    C
C____________________________________________________________________C
C Input Variables :-                                                 C
C ===============                                                    C
C  Integer                                                           C
C  -------                                                           C
C     M         : Well it's M isn't it!                              C
C     L         : Well it's L isn't it!                              C
C  Double Precision                                                  C
C  ----------------                                                  C
C     X         : Cos(theta)                                         C
C     PLMIN1    : P_(l-1)^m ( X )                                    C
C     PLMIN2    : P_(l-2)^m ( X )                                    C
C____________________________________________________________________C
C Local Variables :-                                                 C
C ===============                                                    C
C  Integer                                                           C
C  -------                                                           C
C     RL        : FLOAT ( L )                                        C
C     RM        : FLOAT ( M )                                        C
C____________________________________________________________________C
C
C*********************************************************************
      FUNCTION PLM ( L, M, X, PLMIN1, PLMIN2 )
      IMPLICIT NONE
      DOUBLE PRECISION PLM
C____________________________________________________________________C
C Variable declarations - Parameters ................................C
      INTEGER M,L
      DOUBLE PRECISION X, PLMIN1, PLMIN2
C____________________________________________________________________C
C Variable declarations - Working variables .........................C
      DOUBLE PRECISION RM,RL
C____________________________________________________________________C
C START OF FUNCTION *************************************************C
C____________________________________________________________________C
C
      IF ( L.LT.2 ) THEN
         PRINT *,' You are trying to run function PLM with L < 2.'
         PRINT *,' Program aborted.'
         STOP
      ENDIF
      IF ( M.LT.0 .OR. M.GT.L ) THEN
         PRINT *,' M is out of range in function PLM.'
         PRINT *,' Program aborted.'
         STOP
      ENDIF
      IF ( L.EQ.M ) THEN
         PRINT *,' PLM function called with L = M .'
         PRINT *,' Division by zero would follow.'
         PRINT *,' Program aborted.'
         STOP
      ENDIF
      RM = FLOAT ( M )
      RL = FLOAT ( L )
      PLM = ( 2.0d0*RL - 1.0d0 )*X*PLMIN1
      PLM = PLM - PLMIN2*SQRT( ( RL+RM-1.0d0 )*( RL-RM-1.0d0 ) )
      PLM = PLM/SQRT( (RL+RM)*(RL-RM) )
      RETURN
      END
C*********************************************************************
C*********************************************************************
C function KLM(l,m)=l*(l+1)/2+1+m ... I'm sorry but not even I will  C
C write any more comments for this function!! Steve Gibbons 16.4.97  C
C____________________________________________________________________C
      FUNCTION KLM ( L, M )
      IMPLICIT NONE
      INTEGER KLM,L,M
      KLM = L*(L+1)/2+M+1
      RETURN
      END
C*********************************************************************
C*********************************************************************
C function INDex for Spherical Harmonic Coefficient ******************
C          ---       -         -        -           ******************
C Steve Gibbons 25.4.97                                              C
C____________________________________________________________________C
C Inputs are integers, L and M are obvious ICS = 1 for a cosine harm C
C and ICS = 2 for a sine harmonic.                                   C
C____________________________________________________________________C
      FUNCTION INDSHC ( L, M, ICS )
      IMPLICIT NONE
      INTEGER INDSHC, L, M, ICS
C 
      IF ( M.GT.L ) THEN
         PRINT *,' Function INDSHC. M > L. Program Aborted.'
         PRINT *,' L = ', L
         PRINT *,' M = ', M
         STOP
      ENDIF
      IF ( ICS.NE.1 .AND. ICS.NE.2 ) THEN
         PRINT *,' Function INDSHC. ICS must be 1 or 2.'
         PRINT *,' Program Aborted.'
         STOP
      ENDIF
      IF ( M.EQ.0 ) THEN
         INDSHC = L*L
         RETURN
      ENDIF
      INDSHC = L*L + 2*M - 2 + ICS
      RETURN
      END
C*********************************************************************
C*********************************************************************
C subroutine Spherical Harmonic Scalar File Write ********************
C            -         -        -      -    -     ********************
C Steve Gibbons 13.1.98                                              C
C____________________________________________________________________C
C Writes out the BHCOEF elements into an array compatible with the   C
C nview program ( ie all L * ( L + 2 ) coefficients are written out) C
C____________________________________________________________________C
C Input Variables :-                                                 C
C ===============                                                    C
C  Integer                                                           C
C  -------                                                           C
C     LH	: Level of harmonics.                                C
C     NPOLH 	: Number of scalar harmonics.                        C
C     NPOLMX	: Maximum number of scalar harmonics.                C
C     LMIPOL	: Dimension ( NPOLMX, 3 ) contains, L, M and ICS.    C
C     LU        : Number of file.                                    C
C     IWR       : Mode of writing file. Set to 2 if you do not wish  C
C                  to overwrite a given file. Set to 3 if it doesn't C
C                  matter.                                           C
C                                                                    C
C  Double Precision                                                  C
C  ----------------                                                  C
C     BHCOEF    : Coefficients for the boundary harmonics.           C
C                  Dimension ( NPOLMX )                              C
C  Character                                                         C
C  ---------                                                         C
C     FNAME     : Filename - length not specified                    C
C     TITLE     : Dimension (80).                                    C
C____________________________________________________________________C
C
C*********************************************************************
      SUBROUTINE SHFWR ( LH, NPOLH, NPOLMX, LMIPOL, LU, IWR,
     1                   BHCOEF, FNAME, TITLE )
      IMPLICIT NONE
C____________________________________________________________________C
C Variable declarations - Parameters ................................C
      INTEGER LU, NPOLH, NPOLMX, LMIPOL( NPOLMX, 3 ), IWR, LH
      DOUBLE PRECISION BHCOEF( NPOLMX )
      CHARACTER *( * ) FNAME
      CHARACTER *( 80 ) TITLE
C____________________________________________________________________C
C Variable Declarations - Working Variables .........................C
      INTEGER L, M, ICS, IHARM, ICOL, J
      DOUBLE PRECISION ROW( 5 )
C____________________________________________________________________C
C START OF PROGRAM **************************************************C
C____________________________________________________________________C
C
      IF ( IWR.EQ.2 ) GOTO 500
      IF ( IWR.EQ.3 ) GOTO 500
      WRITE (6, 989)
 989  FORMAT (' Subroutine BHCWT. IWR not 2 or 3. Stop Program. ')
      STOP
c
 500  CONTINUE
c
      CALL FOPEN ( LU, FNAME, IWR )
      WRITE ( LU, 88 ) TITLE
      WRITE ( LU, 89 ) LH

      ICOL = 0
      DO IHARM = 1, LH * ( LH + 2 )
         ICOL = ICOL + 1
         ROW( ICOL ) = 0.0d0
         CALL LMFIND ( IHARM, L, M, ICS )
         DO J = 1, NPOLH
            IF (  LMIPOL( J, 1 ).EQ.L      .AND.
     1            LMIPOL( J, 2 ).EQ.M      .AND.
     2            LMIPOL( J, 3 ).EQ.ICS             ) THEN
               ROW( ICOL ) = BHCOEF( J )
            ENDIF
         ENDDO
         IF ( ICOL.EQ.5 ) THEN
            WRITE ( LU, 90 ) ( ROW( J ), J = 1, 5 )
            ICOL = 0
         ENDIF
      ENDDO
      IF ( ICOL.NE.0 ) WRITE ( LU, 90 ) ( ROW( J ), J = 1, ICOL )
      
      
c
      CALL FCLOSE ( LU, FNAME, TITLE )
c
 88   FORMAT (a80)
 89   FORMAT (i4)
 90   FORMAT (5d16.7)
      RETURN
      END
C*********************************************************************
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
