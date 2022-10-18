C*********************************************************************
C subroutine Kinematic Dynamo (Triangular) Harmonic Selection Routine*
C            -         -       -           -        -         -      *
C Steve Gibbons Fri Jan 28 14:44:05 GMT 2000                         *
C____________________________________________________________________C
C This code is adapted from Graeme Sarson's code to be compatible    C
C with the variables in my code. I've also rewritten it to have no   C
C implicit variables.                                                C
C____________________________________________________________________C
C Input variables :-                                                 C
C ===============                                                    C
C  Integer                                                           C
C  -------                                                           C
C     LHMAX     : Maximum level of harmonics.                        C
C     LH        : Actual level of harmonics.                         C
C     NHMAX     : Maximum number of magnetic field harmonics.        C
C     NVMAX     : Maximum number of velocity field harmonics.        C
C     NV        : Actual number of velocity field harmonics.         C
C     ITRI      : Set to 1 to invoke triangular truncation.          C
C....................................................................C
C  Description of the Seed Magnetic Field Harmonic.                  C
C....................................................................C
C                                                                    C
C     SFIT      : Type = 4 for poloidal and 5 for toroidal.          C
C     SFIL      : Spherical harmonic degree, l.                      C
C     SFIM      : Wavenumber, m, for cos phi dependence.             C
C                 Or (-m) for sin phi dependence.                    C
C                                                                    C
C....................................................................C
C  Description of the Input velocity field.                          C
C....................................................................C
C                                                                    C
C     MVT( IH ) : Type = 1 for poloidal and 2 for toroidal.          C
C     MVL( IH ) : Spherical harmonic degree, l.                      C
C     MVM( IH ) : Wavenumber, m, for cos phi dependence.             C
C                 Or (-m) for sin phi dependence.                    C
C                                                                    C
C____________________________________________________________________C
C Output variables :-                                                C
C ================                                                   C
C  Integer                                                           C
C  -------                                                           C
C     NH        : Actual number of magnetic field harmonics.         C
C                                                                    C
C....................................................................C
C  Description of the induced magnetic field.                        C
C....................................................................C
C                                                                    C
C     MHT( IH ) : Type = 4 for poloidal and 5 for toroidal.          C
C     MHL( IH ) : Spherical harmonic degree, l.                      C
C     MHM( IH ) : Wavenumber, m, for cos phi dependence.             C
C                 Or (-m) for sin phi dependence.                    C
C                                                                    C
C____________________________________________________________________C
C Working Arrays   :-                                                C
C ================                                                   C
C  Integer                                                           C
C  -------                                                           C
C     IQ	: Dimensions( 2, LHMAX, 2 * ( LHMAX + 1 ) )          C
C____________________________________________________________________C
C
C*********************************************************************
      SUBROUTINE KDTHSR( SFIT, SFIL, SFIM, MVT, MVL, MVM,
     1                   MHT, MHL, MHM, IQ, LHMAX, NHMAX,
     2                   ITRI, NVMAX, LH, NH, NV )
      IMPLICIT NONE
C____________________________________________________________________C
C Variable declarations - Parameters ................................C
      INTEGER    SFIT, SFIL, SFIM, ITRI, NVMAX, LH, NH, NV,
     1           LHMAX, NHMAX
      INTEGER    MVT( NVMAX ), MVL( NVMAX ), MVM( NVMAX ),
     1           IQ( 2,LHMAX,2*(LHMAX+1) ), MHT( NHMAX ),
     2           MHL( NHMAX ), MHM( NHMAX )
C____________________________________________________________________C
C Variable declarations - Working variables .........................C
      INTEGER    LMTRI, T, L, M, C, MM, IMM, I, IMAXL, MMM, 
     1           HCOUNT, IQ1, IQ2, IT, IL, IM, IC,
     2           IFEQ1, IFEQ2, IFEQ3, VTI, VLI, VMI, VCI
      LOGICAL ILTRI, IESUML, ISUMZR,
     1        IESUMT, ITTS, IMZSIN, IESUMC, ITWOEQ, OK
C____________________________________________________________________C
C START OF PROGRAM **************************************************C
C____________________________________________________________________C
C
      OK = .TRUE.
C
C********************************************************************
C Set the variables IT, IL, IM and IC to the values they would
C have taken under Graeme's formulation.
C
      IT    = SFIT - 3
      IL    = SFIL
C
      IF ( SFIM.LT.0 ) THEN
        IM = -SFIM
        IC = 1
      ELSE
        IM = SFIM
        IC = 2
      ENDIF
C
      LMTRI = 2*LH
      IF( ITRI.EQ.1 ) LMTRI = LH
      DO 4, T = 1,2
      DO 4, L = 1,LH
      DO 4, M = 1,LH+1
      DO 4, C = 1,2
       MM = M+(C-1)*(LHMAX+1)
    4 IQ(T,L,MM) = 0
C----
C    CHOOSE THE LHS OF THIS EQUATION.
C----
    5 IMM           = IM+1+(LHMAX+1)*(IC-1)
      IQ(IT,IL,IMM) = 2
C-----
C     CHOOSE A VELOCITY 
C     HARMONIC AND CHECK WHICH RHS TERMS ARE RELEVANT.
C-----
      DO 13, I = 1, NV
C       .
C       . Set integer variables VTI, VLI, VMI and VCI
C       . To the values VT( I ), VL( I ), VM( I ) and
C       . VC( I ) would have taken under Graeme Sarson's
C       . formulation.
C       .
        VTI = MVT( I )
        VLI = MVL( I )
        IF ( MVM( I ).LT.0 ) THEN
          VMI = -MVM( I )
          VCI = 1
C
C Remember this is for sin with c.eq.1
C
        ELSE
          VMI = MVM( I )
          VCI = 2
C
C Remember this is for cos with c.eq.2
C
        ENDIF
C       .
       DO 12, L = 1, LH
C
        ILTRI = .FALSE.
        IMAXL = MAX0( VLI, L, IL)
        IF( IMAXL.LE.(VLI+L+IL-IMAXL) ) ILTRI = .TRUE.
C
        IESUML = .FALSE.
        IQ1    = VLI + L + IL
        IQ2    = IQ1/2
        IQ2    = IQ2*2
        IF( IQ1 .EQ. IQ2 ) IESUML = .TRUE.

        DO 11, MM = 1,L+1
         M = MM-1
         IF( L+M.GT.LMTRI ) GOTO 11
C
         ISUMZR = .FALSE.
         IF( VMI+M.EQ.IM.OR.IABS(VMI-M).EQ.IM) ISUMZR = .TRUE.
C
         DO 10, T = 1,2
          IESUMT = .FALSE.
          IQ1 = VTI + T + IT
          IQ2 = IQ1/2
          IQ2 = IQ2*2
          IF( IQ1.EQ.IQ2 ) IESUMT = .TRUE.
C
          ITTS = .FALSE.
          IF( VTI+IT-T.EQ.3 ) ITTS = .TRUE.
C
          DO 9, C = 1,2
           MMM = MM+(LHMAX+1)*(C-1)
C
           IMZSIN = .FALSE.
           IF( M+C.EQ.1 ) IMZSIN = .TRUE.
C
           IESUMC = .FALSE.
           IF( IABS(VCI+C+IC-4).EQ.1 ) IESUMC = .TRUE.
C
           ITWOEQ=.FALSE.
          IFEQ1=IABS( VTI-T )+IABS( VMI-M )+IABS( VLI-L )+IABS( VCI-C )
          IFEQ2=IABS( IT-T )+IABS( L-IL )+IABS( IM-M )+IABS( IC-C )
          IFEQ3=IABS( VTI-IT )+IABS( VMI-IM )+IABS( VLI-IL )+
     +                  IABS( VCI-IC )
          IF( IFEQ1.EQ.0.OR.IFEQ2.EQ.0.OR.IFEQ3.EQ.0 ) ITWOEQ = .TRUE.
C
          IF( ILTRI.AND.ISUMZR.AND..NOT.IMZSIN ) THEN
C
            IF(((.NOT.IESUMT.AND.IESUML.AND..NOT.IESUMC.AND.
     *          .NOT.ITTS).OR.(IESUMT.AND..NOT.IESUML.AND.
     *          IESUMC.AND..NOT.ITWOEQ)).AND.IQ(T,L,MMM).NE.2) 
     *                       IQ(T,L,MMM)=1
C
           ENDIF
C
    9     CONTINUE
   10    CONTINUE
   11   CONTINUE
   12  CONTINUE
   13 CONTINUE
C-----
C     COUNT THE TOTAL NUMBER OF HARMONIC, AND ORDER THEM
C-----
      HCOUNT=0
      DO 14, T  = 1, 2
      DO 14, L  = 1, LH
      DO 14, MM = 1, L+1
        M = MM-1
      DO 14, C = 1,2
        MMM = MM+(LHMAX+1)*(C-1)
        IF( IQ(T,L,MMM).EQ.1 ) THEN
          IT = T
          IL = L
          IM = M
          IC = C
          GOTO 5
        ENDIF
      IF( IQ(T,L,MMM).EQ.2 ) HCOUNT = HCOUNT+1
   14 CONTINUE
C--------------------------------------------------------
C   Fill the arrays MHT, MHL, MHM with bounds checking. |
C--------------------------------------------------------
      NH = 0
      DO 15, T  = 1, 2
      DO 15, L  = 1, LH
      DO 15, MM = 1, L+1
        M = MM-1
      DO 15, C = 1,2
       MMM = MM+(LHMAX+1)*(C-1)
       IF( IQ(T,L,MMM).NE.0 ) THEN
         CALL CNTRIC( NH, NHMAX, OK )
         IQ(T,L,MMM)=NH
         IF ( OK )              MHT( NH ) = T + 3
         IF ( OK )              MHL( NH ) = L
         IF ( OK .AND. C.EQ.1 ) MHM( NH ) = -M
         IF ( OK .AND. C.EQ.2 ) MHM( NH ) = M
C        .
C        . Remember that c.eq.1 means sin, and
C        .               c.eq.2 means cos.
C        .
       ENDIF
 15   CONTINUE
C
C   IG = GAMMA FROM B.& G., IA = ALPHA, ETC
C
C-----
      IF ( OK ) RETURN
C
C Too many harmonics were required.
C
      PRINT *,' Subroutine KDTHSR.'
      PRINT *,' Number of harmonics required = ',NH
      PRINT *,' Maximum possible number = ', NHMAX
      PRINT *,' Program aborted.'
      STOP
      END
C********************************************************************
