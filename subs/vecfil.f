C*********************************************************************
C subroutine VECtor FILl *********************************************
C            ---    ---  *********************************************
C Steve Gibbons 16.3.99                                              C
C____________________________________________________________________C
C                                                                    C
C____________________________________________________________________C
C Input variables :-                                                 C
C ===============                                                    C
C  Integer                                                           C
C  -------                                                           C
C     NDIM      : Length of the solution vector.                     C
C     NR        : Number of radial grid nodes.                       C
C     NH        : Number of harmonics (all types)                    C
C                                                                    C
C     MHT       : Dimension ( NH )                                   C
C     MHL       : Dimension ( NH )                                   C
C     MHM       : Dimension ( NH )                                   C
C     MHC       : Dimension ( NH )                                   C
C                                                                    C
C  mht, mhl, mhm and mhc define the spherical harmonics present      C
C  in the solution vector. For spherical harmonic number I;          C
C                                                                    C
C   MHT( I ) = 1 for a poloidal velocity vector                      C
C   MHT( I ) = 2 for a toroidal velocity vector                      C
C   MHT( I ) = 3 for a temperature / codensity term                  C
C                                                                    C
C   MHL( I ) = spherical harmonic degree, l                          C
C                                                                    C
C   MHM( I ) = spherical harmonic order, m                           C
C                                                                    C
C   MHC( I ) = 1 for a cosine dependence in phi and                  C
C            = 2  "  "  sine     "        "  "                       C
C                                                                    C
C     IFN       : Function to be entered. Currently                  C
C                 ifn = 1 ---> f(r) = cos( r )                       C
C                 ifn < 0 ---> a trigonometric function which        C
C                automatically satisfies all BCs is chosen with      C
C                multiple of IABS( IFN ) ...                         C
C                                                                    C
C     IMODE     : Mode of choice for output                          C
C             IMODE = 0 fills every harmonic with function IFN.      C
C             IMODE = -1 fills every poloidal harmonic with IFN.     C
C             IMODE = -2 fills every toroidal harmonic with IFN.     C
C             IMODE = -3 fills every theta    harmonic with IFN.     C
C             IMODE .gt. 0 fills only harmonic number IMODE          C
C                                                                    C
C  Double Precision                                                  C
C  ----------------                                                  C
C     RI        : Radius of inner boundary.                          C
C     RO        : Radius of outer boundary.                          C
C     RVEC      : DP vector of dimension ( NDIM )                    C
C                                                                    C
C____________________________________________________________________C
C
C*********************************************************************
      SUBROUTINE VECFIL ( NDIM, NR, NH, MHT, MHL, MHM, MHC, RI, RO,
     1                    RVEC, IFN, IMODE, BCIVEL, BCOVEL, BCITHE,
     2                    BCOTHE )
      IMPLICIT NONE
C____________________________________________________________________C
C Variable declarations - Parameters ................................C
      INTEGER NDIM, NR, NH, IFN, IMODE,
     1        MHT( NH ), MHM( NH ), MHC( NH ), MHL( NH )
      DOUBLE PRECISION RI, RO, RVEC( NDIM )
      CHARACTER *(2) BCIVEL, BCOVEL, BCITHE, BCOTHE
C____________________________________________________________________C
C Variable declarations - Working variables .........................C
      DOUBLE PRECISION H, RAD, PI, A, FFUN, GFUN, AFUN, BFUN, LAMB,
     1                 CONST, X
      PARAMETER (PI=3.14159265358979312D0, LAMB=4.73004074d0)
      INTEGER IRN, IH, IND, NN
      LOGICAL OK
C____________________________________________________________________C
C START OF PROGRAM **************************************************C
C____________________________________________________________________C
C
      CONST = 2.0d0*LAMB*LAMB
C
C Check the length of the vector
C
      IF ( NDIM.NE.NR*NH ) THEN
         PRINT *,' Subroutine VECFIL. Bad dimensions '
         PRINT *,' Array length = ',NDIM
         PRINT *,' NR = ', NR,' NH = ',NH
         STOP
      ENDIF
C
C Check the character variables
C
      IF ( BCITHE.NE.'HF' .AND. BCITHE.NE.'TM' ) THEN
         PRINT *,' Sub. VECFIL. BCITHE = ', BCITHE
      ENDIF
      IF ( BCOTHE.NE.'HF' .AND. BCOTHE.NE.'TM' ) THEN
         PRINT *,' Sub. VECFIL. BCOTHE = ', BCOTHE
      ENDIF
      IF ( BCIVEL.NE.'SF' .AND. BCIVEL.NE.'NS' ) THEN
         PRINT *,' Sub. VECFIL. BCIVEL = ', BCIVEL
      ENDIF
      IF ( BCOVEL.NE.'SF' .AND. BCOVEL.NE.'NS' ) THEN
         PRINT *,' Sub. VECFIL. BCOVEL = ', BCOVEL
      ENDIF
C
      H = ( RO - RI ) / DBLE( NR - 1 )
      NN = IABS( IFN )
      IF ( BCIVEL.NE.BCOVEL .OR.
     1     BCITHE.NE.BCOTHE ) NN = 1
C
      DO IH = 1, NH
        OK = .FALSE.
        IF ( IMODE.EQ.0 ) OK = .TRUE.
        IF ( IMODE.EQ.IH ) OK = .TRUE.
        IF ( (IMODE + MHT( IH )).EQ.0 ) OK = .TRUE.
        IF ( .NOT. OK ) GOTO 500
        DO IRN = 1, NR
          RAD = RI + H*( IRN - 1 )
          IND = ( IRN - 1 )*NH + IH
C-----------------------------------------------------
C  DOING HARMONIC IH at grid node IRN
C
C  First simple option: IFN = 1
C  f(r) = cos( r )
C
          IF ( IFN.EQ.1 ) RVEC( IND ) = DCOS( RAD )
C
C Second option, f(r) = rad
C
          IF ( IFN.EQ.2 ) RVEC( IND ) = RAD
C
C Third option, f(r) = rad^2/2
C
          IF ( IFN.EQ.3 ) RVEC( IND ) = 0.5d0*RAD*RAD
C
C Fourth option, f(r) = rad^3/3
C
          IF ( IFN.EQ.4 ) RVEC( IND ) = RAD*RAD*RAD/3.0d0
C
C Fifth option, f(r) = rad^4/4
C
          IF ( IFN.EQ.4 ) RVEC( IND ) = RAD*RAD*RAD*RAD/4.0d0
C
C Option IFN < 0
C
          IF ( IFN.LT.0 ) THEN
c            .
             A = PI*DBLE( NN )/(RO-RI)
c            .
c            . First, poloidal harmonics
c            .
             IF ( MHT( IH ).EQ.1 ) THEN
               IF ( BCIVEL.EQ.'SF' .AND. BCOVEL.EQ.'SF' ) THEN
                  RVEC( IND ) = DSIN( A*(RAD - RI) )
               ENDIF
               IF ( BCIVEL.EQ.'NS' .AND. BCOVEL.EQ.'NS' ) THEN
                  X = (RAD-RI)/(RO-RI) - 0.5d0
                  RVEC( IND ) = COSH( LAMB*X )/COSH( 0.5d0*LAMB) -
     1                          DCOS( LAMB*X )/DCOS( 0.5d0*LAMB)
               ENDIF
               IF ( BCIVEL.EQ.'SF' .AND. BCOVEL.EQ.'NS' ) THEN
                  X = (RAD-RI)/(RO-RI) - 0.5d0
                  AFUN = X - 0.5d0
                  BFUN = (-2.0d0)*DBLE(NN)*PI/CONST
                  FFUN = DSIN( NN*PI*( X + 0.5d0 )  )
                  GFUN = COSH( LAMB*X )/COSH( 0.5d0*LAMB) -
     1                    DCOS( LAMB*X )/DCOS( 0.5d0*LAMB)
                  RVEC( IND ) = AFUN*FFUN + BFUN*GFUN
               ENDIF
               IF ( BCIVEL.EQ.'NS' .AND. BCOVEL.EQ.'SF' ) THEN
                  X = (RAD-RI)/(RO-RI) - 0.5d0
                  AFUN = X + 0.5d0
                  BFUN = 2.0d0*DBLE(NN)*PI/CONST
                  FFUN = DSIN( NN*PI*( X + 0.5d0 )  )
                  GFUN = COSH( LAMB*X )/COSH( 0.5d0*LAMB) -
     1                    DCOS( LAMB*X )/DCOS( 0.5d0*LAMB)
                  RVEC( IND ) = AFUN*FFUN + BFUN*GFUN
               ENDIF
             ENDIF
c            .
c            . Secondly, toroidal harmonics
c            .
             IF ( MHT( IH ).EQ.2 ) THEN
                FFUN = RAD*DCOS( A*(RAD - RI) )
                GFUN = DSIN( A*(RAD - RI) )
                IF ( BCIVEL.EQ.'SF' .AND. BCOVEL.EQ.'SF' ) THEN
                   AFUN = 1.0d0
                   BFUN = 0.0d0
                ENDIF
                IF ( BCIVEL.EQ.'NS' .AND. BCOVEL.EQ.'NS' ) THEN
                   AFUN = 0.0d0
                   BFUN = 1.0d0
                ENDIF
                IF ( BCIVEL.EQ.'NS' .AND. BCOVEL.EQ.'SF' ) THEN
                   AFUN = (RAD - RI)
                   BFUN = -1.0d0*RO/A
                ENDIF
                IF ( BCIVEL.EQ.'SF' .AND. BCOVEL.EQ.'NS' ) THEN
                   AFUN = (RAD - RO )
                   BFUN = -1.0d0*RI/A
                ENDIF
                RVEC( IND ) = AFUN*FFUN + BFUN*GFUN
                IF ( BCIVEL.EQ.'SF' .AND. BCOVEL.EQ.'SF' .AND.
     1               MHL( IH ).EQ.1 .AND. MHM( IH ).EQ.0 .AND.
     2               MHC( IH ).EQ.1                    ) THEN
                  RVEC( IND ) = 0.0d0
                ENDIF
             ENDIF
c            .
c            . Lastly, temperature harmonics
c            .
             IF ( MHT( IH ).EQ.3 ) THEN
                IF ( BCITHE.EQ.'TM' .AND. BCOTHE.EQ.'TM' ) THEN
                   RVEC( IND ) = DSIN( A*(RAD - RI) )
                ENDIF
                IF ( BCITHE.EQ.'HF' .AND. BCOTHE.EQ.'TM' ) THEN
                   RVEC( IND ) = DCOS( 0.5d0*A*(RAD - RI) )
                ENDIF
                IF ( BCITHE.EQ.'TM' .AND. BCOTHE.EQ.'HF' ) THEN
                   RVEC( IND ) = DSIN( 0.5d0*A*(RAD - RI) )
                ENDIF
                IF ( BCITHE.EQ.'HF' .AND. BCOTHE.EQ.'HF' ) THEN
                   RVEC( IND ) = DCOS( A*(RAD - RI) )
                ENDIF
             ENDIF
c            .
          ENDIF
C
C-----------------------------------------------------
        ENDDO
 500    CONTINUE
      ENDDO
C
      RETURN
      END
C*********************************************************************
