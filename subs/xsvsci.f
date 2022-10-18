C*********************************************************************
C subroutine Xtra Special Vector Single Component Integrate **********
C            -    -       -      -      -         -         **********
C Steve Gibbons Mon Apr  9 12:18:19 MET DST 2001                     C
C____________________________________________________________________C
C                                                                    C
C The function XSV( NCMX, NPHP, NTHP, NR ) defines, in real space,   C
C the values of a function. SF( ICMP, iphi, ithe, ir ) is the value  C
C at  RAD = xarr( ir )                                               C
C                                                                    C
C  THETA = ACOS[ GAUX( ithe ) ]                                      C
C                                 and                                C
C  PHI   = (iphi-1)*DELTAP                                           C
C                                 with                               C
C  deltap = 2*pi/(NPHP*M0).                                          C
C                                                                    C
C XSVSCI will return the value DINT, the volume integral (over the   C
C full spherical shell or sphere) of the function stored in the      C
C component ICMP.                                                    C
C                                                                    C
C____________________________________________________________________C
C                                                                    C
C Input variables :-                                                 C
C ===============                                                    C
C                                                                    C
C  Integer                                                           C
C  -------                                                           C
C                                                                    C
C     NCMX      : Maximum number of components stored in XSV         C
C     NPHP      : Number of phi points.                              C
C     NTHP      : Number of theta points.                            C
C     NR        : Total number of radial grid nodes.                 C
C     ICM       : Component of XSV which stores scalar function.     C
C                                                                    C
C  Double Precision                                                  C
C  ----------------                                                  C
C                                                                    C
C     XARR      : Array of radial spacing values. Dim ( NR )         C
C                                                                    C
C     GAUW      : Gauss weights as calculated by GAUWTS. ( NTHP )    C
C                                                                    C
C     XSV       : Xtra Special Function: is an array containing a    C
C                  function over a set of theta, phi and r points.   C
C                    Dimensions are                                  C
C                      ( NCMX, NPHP, NTHP, NR )                      C
C                                                                    C
C     DINT      : Spherical integral of function.                    C
C                                                                    C
C____________________________________________________________________C
C
C*********************************************************************
      SUBROUTINE XSVSCI( NCMX, NPHP, NTHP, NR, ICM, XARR, GAUW,
     1                   XSV, DINT )
      IMPLICIT NONE
C____________________________________________________________________C
C Variable declarations - Parameters ................................C
      INTEGER          NCMX, NPHP, NTHP, NR, ICM
      DOUBLE PRECISION GAUW( NTHP ), DINT, XARR( NR ),
     1                 XSV( NCMX, NPHP, NTHP, NR )
C____________________________________________________________________C
C Variable declarations - Working variables .........................C
      INTEGER          ITHE, IPHI, IR
      DOUBLE PRECISION PI, ZERO, ZCOEF, WEIGHT, RAD, FAC1, FAC2,
     1                 VFAC, DRAD
      PARAMETER        ( ZERO = 0.0d0, PI = 3.14159265358979312D0 )
C____________________________________________________________________C
C START OF PROGRAM **************************************************C
C____________________________________________________________________C
C
C Check value of ICM
C
      IF ( ICM.LT.1 .OR. ICM.GT.NCMX ) THEN
        PRINT *,' Subroutine XSVSCI.'
        PRINT *,' ICM = ', ICM,' NCMX = ', NCMX
        PRINT *,' Program aborted.'
        STOP
      ENDIF
C
      DINT = ZERO
C
      FAC1 = 0.5d0/DBLE( NPHP )
C
C ............... begin looping around radial grid nodes .............
      DO IR = 1, NR
        RAD  = XARR( IR )
        VFAC = 4.0d0*PI*RAD**2
        FAC2 = ZERO
C
C ............... begin looping around theta points ..................
        DO ITHE = 1, NTHP
          WEIGHT = GAUW( ITHE )*FAC1
          ZCOEF  = ZERO
          DO IPHI = 1, NPHP
            ZCOEF = ZCOEF + XSV( ICM, IPHI, ITHE, IR )
          ENDDO
          FAC2 = FAC2 + WEIGHT*ZCOEF
        ENDDO
C       .
C       . End looping around theta points
C       .
        IF ( IR.EQ.1 ) THEN
          DRAD = 0.5d0*(  XARR( 2 ) - XARR( 1 )  )
        ENDIF
C       .
        IF ( IR.GT.1 .AND. IR.LT.NR ) THEN
          DRAD = 0.5d0*(  XARR( IR+1 ) - XARR( IR-1 )  )
        ENDIF
C       .
        IF ( IR.EQ.NR ) THEN
          DRAD = 0.5d0*(  XARR( NR ) - XARR( NR-1 )  )
        ENDIF
C       .
        DINT = DINT + DRAD*FAC2*VFAC
C       .
      ENDDO
C     .
C     . End looping around radial grid nodes
C     .
      RETURN
      END
C*********************************************************************
