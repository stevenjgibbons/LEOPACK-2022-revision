C*********************************************************************
C subroutine VELocity Rhs Boundary cond. Enforce *********************
C            ---      -   -              -       *********************
C Steve Gibbons Mon Oct 18 16:04:01 BST 1999                         C
C____________________________________________________________________C
C                                                                    C
C Imposes the boundary conditions upon the right hand side of        C
C the equation A x = rhs, once x has been treated by VELMBE.         C
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
C  Double Precision                                                  C
C  ----------------                                                  C
C                                                                    C
C     RHS       : Right hand vector. Dim ( N2 )                      C
C                                                                    C
C  Character                                                         C
C  ---------                                                         C
C                                                                    C
C     CHBND     : Boundary flag (*). Either 'Inner' or 'Outer'       C
C                                                                    C
C____________________________________________________________________C
C
C*********************************************************************
      SUBROUTINE VELRBE( N2, NR, INARR, MHTO, RHS, CHBND)
      IMPLICIT NONE
C____________________________________________________________________C
C Variable declarations - Parameters ................................C
      INTEGER N2, NR, INARR( * ), MHTO( * )
      DOUBLE PRECISION RHS( N2 )
      CHARACTER *(*)   CHBND
C____________________________________________________________________C
C Variable declarations - Working variables .........................C
      INTEGER NRR, NH, IHR, IRR, IRR2, NRC, INDFUN
      DOUBLE PRECISION ZERO
      PARAMETER ( ZERO = 0.0d0 )
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
        PRINT *,' Subroutine VELRBE.'
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
        IRR2  = 2
        GOTO 50
      ENDIF
C     .
      IF ( CHBND(1:1).EQ.'O' .OR. CHBND(1:1).EQ.'o' ) THEN
        IRR2  = NR - 1
        IRR   = NR
        GOTO 50
      ENDIF
C     .
      PRINT *,' Subroutine VELRBE. CHBND = ', CHBND
      PRINT *,' Program aborted.'
      STOP
C     .
 50   CONTINUE
C     .
      DO IHR = 1, NH
        IF ( MHTO( IHR ).NE.1 .AND. MHTO( IHR ).NE.2 ) GOTO 60
C       .
C       . Harmonic is either poloidal or toroidal velocity
C       .
        NRC = INDFUN( IRR, IHR, INARR )
        RHS( NRC ) = ZERO
C       .
        IF ( MHTO( IHR ).EQ.2 ) THEN
          NRC = INDFUN( IRR2, IHR, INARR )
          RHS( NRC ) = ZERO
        ENDIF
C       .
 60   CONTINUE
      ENDDO
C     .
      RETURN
      END
C*********************************************************************

