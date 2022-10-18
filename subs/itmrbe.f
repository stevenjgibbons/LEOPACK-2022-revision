C*********************************************************************
C subroutine Inhomogeneous TeMperature Rhs Boundary cond. Enforce ****
C            -             -           -   -              -       ****
C Steve Gibbons Mon Oct 18 16:04:01 BST 1999                         C
C____________________________________________________________________C
C                                                                    C
C Imposes the boundary conditions upon the right hand side of        C
C the equation A x = rhs, once x has been treated by TMPMBE.         C
C The coefficients for the inhomogeneous temperature of heat flux    C
C boundary conditions are given by the array SHC and the scalar      C
C ZCOEF (for the P^0_0 term).                                        C
C                                                                    C
C____________________________________________________________________C
C Input variables :-                                                 C
C ===============                                                    C
C  Integer                                                           C
C  -------                                                           C
C     N2        : Length of r.h.s. vector.                           C
C     NR        : Number of radial grid nodes.                       C
C                                                                    C
C     INARR     : Array of integers. Dimension ( * )                 C
C                 At present the elements of INARR are as follows.   C
C                                                                    C
C                 INARR( 1 ) = IFORMF. Format flag.                  C
C                  This flag decides how the elements are arranged   C
C                  in the solution vector.                           C
C                                                                    C
C                  The current options are:-                         C
C                                                                    C
C                   IFORMF = 1, 3. INDFUN = ( IR - 1 )*NH + IH       C
C                   IFORMF = 2, 4. INDFUN = ( IH - 1 )*NR + IR       C
C                                                                    C
C  where IR and IH are the current grid node and harmonic resp.      C
C  and NR and NH are the total numbers of nodes and harmonics        C
C  in the solution vector.                                           C
C                                                                    C
C                 INARR( 2 ) = NRR. Number of radial grid nodes.     C
C                    (NRR must be consistent with NR above).         C
C                 INARR( 3 ) = NH. Number of harmonics in sol. vect. C
C                                                                    C
C MHTO, MHLO and MHMO all describe the harmonics of the equations    C
C (i.e. describe the rows of the matrix)                             C
C                                                                    C
C     MHTO      : Array length ( * ) - atleast length NH             C
C                                                                    C
C         MHTI( IH ) = 1 --> harmonic is poloidal velocity           C
C         MHTI( IH ) = 2 --> harmonic is toroidal velocity           C
C         MHTI( IH ) = 3 --> harmonic is temperature.                C
C         MHTI( IH ) = 4 --> harmonic is poloidal magnetic field     C
C         MHTI( IH ) = 5 --> harmonic is toroidal magnetic field     C
C                                                                    C
C     MHLO      : Array length ( * ) - atleast length NH             C
C                  Sph. harm. degree, l.                             C
C     MHMO      : Array length ( * ) - atleast length NH             C
C                  Sph. harm. order, m, for cos m phi dep.           C
C                 -Sph. harm. order, m, for sin m phi dep.           C
C                                                                    C
C     LH        : Maximum spherical harmonic degree, l.              C
C                                                                    C
C  Double Precision                                                  C
C  ----------------                                                  C
C                                                                    C
C     RHS       : Right hand vector. Dim ( N2 )                      C
C     ZCOEF     : Coefficient of P_0^0 term.                         C
C     SHC       : Coefficients for non zero l, m                     C
C                 Dimension ( LH*( LH + 2 ) )                        C
C                                                                    C
C                  Coeff for P_l^m{ cos m phi } is given by          C
C                   SHC( INDSHC( L, M, 1 )                           C
C                                                                    C
C                  Coeff for P_l^m{ sin m phi } is given by          C
C                   SHC( INDSHC( L, M, 2 )                           C
C                                                                    C
C  Character                                                         C
C  ---------                                                         C
C                                                                    C
C     CHBND     : Boundary flag (*). Either 'Inner' or 'Outer'       C
C                                                                    C
C____________________________________________________________________C
C
C*********************************************************************
      SUBROUTINE ITMRBE( N2, NR, INARR, MHTO, MHLO, MHMO,
     1                   RHS, CHBND, LH, ZCOEF, SHC )
      IMPLICIT NONE
C____________________________________________________________________C
C Variable declarations - Parameters ................................C
      INTEGER N2, NR, INARR( * ), MHTO( * ), MHLO( * ), MHMO( * ), LH
      DOUBLE PRECISION RHS( N2 ), ZCOEF, SHC( LH*(LH+2) )
      CHARACTER *(*)   CHBND
C____________________________________________________________________C
C Variable declarations - Working variables .........................C
      INTEGER NRR, NH, IHR, IRR, NRC, INDFUN, L, M, ICS, IHARM,
     1        INDSHC
      DOUBLE PRECISION FAC
      EXTERNAL MHDBCS
C____________________________________________________________________C
C START OF PROGRAM **************************************************C
C____________________________________________________________________C
C     .
C     . Check value of NR
C     .
      NRR    = INARR( 2 )
      NH     = INARR( 3 )
C     .
      IF ( NRR.NE.NR ) THEN
        PRINT *,' Subroutine ITMRBE.'
        PRINT *,' NR  = ', NR
        PRINT *,' NRR = ', NRR
        PRINT *,' Program aborted.'
        STOP
      ENDIF
C     .
C     . Check the validity of CHBND
C     .
      IF ( CHBND(1:1).EQ.'I' .OR. CHBND(1:1).EQ.'i' ) THEN
        IRR   = 1
        GOTO 50
      ENDIF
C     .
      IF ( CHBND(1:1).EQ.'O' .OR. CHBND(1:1).EQ.'o' ) THEN
        IRR   = NR
        GOTO 50
      ENDIF
C     .
      PRINT *,' Subroutine ITMRBE. CHBND = ', CHBND
      PRINT *,' Program aborted.'
      STOP
C     .
 50   CONTINUE
C     .
      DO IHR = 1, NH
        IF ( MHTO( IHR ).NE.3 ) GOTO 60
C       .
C       . Harmonic is a temperature function
C       .
        L   = MHLO( IHR )
        IF ( MHMO( IHR ).GE.0 ) THEN
          ICS = 1
          M   = MHMO( IHR )
        ELSE
          ICS = 2
          M   = MHMO( IHR )*(-1)
        ENDIF
C       .
        IHARM = INDSHC( L, M, ICS )
        IF ( IHARM.EQ.0 ) THEN
          FAC = ZCOEF
        ELSE
          FAC = SHC( IHARM )
        ENDIF
C       .
        NRC = INDFUN( IRR, IHR, INARR )
        RHS( NRC ) = FAC
C       .
 60   CONTINUE
      ENDDO
C     .
      RETURN
      END
C*********************************************************************

