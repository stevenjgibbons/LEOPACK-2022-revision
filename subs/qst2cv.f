C********************************************************************
C subroutine QST to Solution Vector *********************************
C            ---    -        -      *********************************
C Steve Gibbons 19.7.99                                             C
C___________________________________________________________________C
C Reads in QST array and converts values into CV vector.            C
C___________________________________________________________________C
C Input Variables :-                                                C
C ===============                                                   C
C  Integers                                                         C
C  --------                                                         C
C     LH        : Maximum degree of spherical harmonics             C
C     NDIM      : Length of the solution vector.                    C
C     NR        : Number of radial grid nodes.                      C
C     NH        : Number of harmonics (all types)                   C
C     IVLMF     : Velocity / Magnetic field switch.                 C
C                                                                   C
C                  ivlmf = 1 --> put in velocity entries            C
C                  ivlmf = 2 --> put in magnetic field entries      C
C                                                                   C
C     MHT       : Dimension ( NH )                                  C
C     MHL       : Dimension ( NH )                                  C
C     MHM       : Dimension ( NH )                                  C
C     MHC       : Dimension ( NH )                                  C
C                                                                   C
C  mht, mhl, mhm and mhc define the spherical harmonics present     C
C  in the solution vector. For spherical harmonic number I;         C
C                                                                   C
C   MHT( I ) = 1 for a poloidal velocity vector                     C
C   MHT( I ) = 2 for a toroidal velocity vector                     C
C   MHT( I ) = 3 for a temperature / codensity term                 C
C                                                                   C
C   MHL( I ) = spherical harmonic degree, l                         C
C                                                                   C
C   MHM( I ) = spherical harmonic order, m                          C
C                                                                   C
C   MHC( I ) = 1 for a cosine dependence in phi and                 C
C            = 2  "  "  sine     "        "  "                      C
C                                                                   C
C  Double Precision                                                 C
C  ----------------                                                 C
C     QST       : Array containing the QST coefficients             C
C     RAD       : Value of radius.                                  C
C     CV        : DP vector of dimension ( NDIM )                   C
C___________________________________________________________________C
C
C********************************************************************
      SUBROUTINE QST2CV ( LH, QST, CV, IRN, NR, NH, NDIM, IVLMF,
     1                    MHT, MHL, MHM, MHC, RAD )
      IMPLICIT NONE
C___________________________________________________________________C
C Variable declarations - Parameters ...............................C
      INTEGER LH, IRN, NR, NH, NDIM, IVLMF,
     1        MHT( NH ), MHL( NH ), MHM( NH ), MHC( NH )
      DOUBLE PRECISION RAD, QST ( LH*(LH + 2) , 3), CV( NDIM )
C___________________________________________________________________C
C Variable declarations - Working Variables ........................C
      INTEGER L, M, ICS, IH, IT, NOH, INDSHC, IND, IPOL, ITOR
      DOUBLE PRECISION RLL1, SQRLL1, FAC
C___________________________________________________________________C
C START OF PROGRAM *************************************************C
C___________________________________________________________________C
C
      IF ( IVLMF.NE.1 .AND. IVLMF.NE.2 ) THEN
         PRINT *,' Subroutine QST2CV. IVLMF = ', IVLMF
         PRINT *,' Must be either one or two.'
         STOP
      ENDIF
C
      IF ( IVLMF.EQ.1 ) THEN
        IPOL = 1
        ITOR = 2
      ENDIF
C
      IF ( IVLMF.EQ.2 ) THEN
        IPOL = 4
        ITOR = 5
      ENDIF
C
C Check the length of the vector
C
      IF ( NDIM.NE.NR*NH ) THEN
         PRINT *,' Subroutine QST2CV. Bad dimensions '
         PRINT *,' Array length = ',NDIM
         PRINT *,' NR = ', NR,' NH = ',NH
         STOP
      ENDIF
C
C Loop around all harmonics ....
C
      DO IH = 1, NH
        IT   = MHT( IH )
        L    = MHL( IH )
        M    = MHM( IH )
        ICS  = MHC( IH )
        NOH  = INDSHC( L, M, ICS )
        RLL1 = DBLE( L*L + L )
        SQRLL1 = DSQRT( RLL1 )
        IND = ( IRN - 1 )*NH + IH
C
C Do poloidal harmonic
C
        IF ( IT.EQ.IPOL ) THEN
           FAC = QST( NOH, 1) * RAD / RLL1
           CV( IND ) = FAC
        ENDIF
C
C Do toroidal harmonic
C
        IF ( IT.EQ.ITOR ) THEN
           FAC = ( -1.0d0 )*QST( NOH, 3) / SQRLL1
           CV( IND ) = FAC
        ENDIF
C
      ENDDO
C
      RETURN
      END
C********************************************************************
