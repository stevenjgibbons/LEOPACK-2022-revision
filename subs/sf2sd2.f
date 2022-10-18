C*********************************************************************
C subroutine Scalar Function 2 Spectral Decomposition add Xsv ********
C            -      -        - -        -                 -   ********
C Steve Gibbons Thu Nov 30 15:24:44 MET 2000                         C
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
C                                                                    C
C SV3 is a solution vector containing a spectral decomposition       C
C of a scalar function, psi, at NR radial grid nodes such that       C
C                                                                    C
C  psi = \sum_{ih3 = 1, nh3} psi_{ih3}(r) Y_{ih3}                    C
C                                                                    C
C  Y_{ih3} = P_L^M( \cos \theta ) cos ( M \phi ) for MM3( ih3 ) = M  C
C     or                                                             C
C  Y_{ih3} = P_L^M( \cos \theta ) sin ( M \phi ) for MM3( ih3 ) = -M C
C                                                                    C
C    where L = ML3( ih3 ) and                                        C
C    the P_L^M are Schmidt quasi-normalised associated               C
C    Legendre functions (see routine SCHNLA ).                       C
C                                                                    C
C The radial functions psi_{ih3}(r) are stored in the vector SV3     C
C with psi_{ih3}( ir ) being stored in the element                   C
C                                                                    C
C  SV3[ ( ih3 - 1 )*NR + IR ]                                        C
C                                                                    C
C If the array SF denotes values of a function at grid node IR       C
C then SF2SD2 will add FAC*SF to the appropriate elements in SV3.    C
C                                                                    C
C This routine requires the array SF2SA to be pre-calculated         C
C by a call to SF2SD1 (this array must incorporate the scaling FAC). C
C                                                                    C
C____________________________________________________________________C
C                                                                    C
C Input variables :-                                                 C
C ===============                                                    C
C                                                                    C
C  Integer                                                           C
C  -------                                                           C
C                                                                    C
C     ILNR      : First radial grid node to be acted upon.           C
C     IRNR      : Final radial grid node to be acted upon.           C
C     NR        : Total number of radial grid nodes.                 C
C     NH3       : Number of scalar functions                         C
C     ML3       : Array dim ( NH3 ). Sph. harm degree, L.            C
C     MM3       : Array dim ( NH3 ). Sph. harm order, M, or -M.      C
C     M0        : Smallest non-zero wavenumber in solution.          C
C     NTHP      : Number of theta points.                            C
C     NPHP      : Number of phi points.                              C
C     NCMX      : Maximum number of components stored in XSV         C
C     ICM       : Component of XSV which stores scalar function.     C
C                                                                    C
C  Double Precision                                                  C
C  ----------------                                                  C
C                                                                    C
C     XSV       : Xtra Special Function: is an array containing a    C
C                  function over a set of theta, phi and r points.   C
C                    Dimensions are                                  C
C                      ( NCMX, NPHP, NTHP, NR )                      C
C                                                                    C
C     SV3       : Dim ( NR*NH3 ) Radial functions psi_{ip3}(r)       C
C                                                                    C
C     FTF       : Work array: dim ( 2*NPHP )                         C
C                                                                    C
C     SF2SA     : Dim ( NH3, NTHP ). Array prepared by SF2SD1.       C
C____________________________________________________________________C
C
C*********************************************************************
      SUBROUTINE SF2SD2( ILNR, IRNR, NR, NH3, ML3, MM3, M0, NTHP,
     1                   NPHP, NCMX, ICM, SF, SV3, FTF, SF2SA )
      IMPLICIT NONE
C____________________________________________________________________C
C Variable declarations - Parameters ................................C
      INTEGER          ILNR, IRNR, NR, NH3, ML3( NH3 ),
     1                 MM3( NH3 ), M0, NTHP, NPHP, NCMX, ICM
      DOUBLE PRECISION SF( NCMX, NPHP, NTHP, NR ),
     1                 FTF( 2*NPHP ), SV3( * ), SF2SA( NH3, NTHP )
C____________________________________________________________________C
C Variable declarations - Working variables .........................C
      INTEGER          ISIGN, ITHETA, LENV, IPHI, INDLOC, IH3,
     1                 L, M, INDEX, IR, IRMNR
      DOUBLE PRECISION ZERO, ZCOEF, FAC2
      PARAMETER        ( ZERO = 0.0d0, ISIGN = 1 )
C____________________________________________________________________C
C START OF PROGRAM **************************************************C
C____________________________________________________________________C
C
C Check value of ICM
C
      IF ( ICM.LT.1 .OR. ICM.GT.NCMX ) THEN
        PRINT *,' Subroutine SF2SD2.'
        PRINT *,' ICM = ', ICM
        PRINT *,' Program aborted.'
        STOP
      ENDIF
C
      LENV = 2*NPHP
C ............... begin looping around radial grid nodes .............
      DO IR = ILNR, IRNR
       IRMNR = IR - NR
C
C ............... begin looping around theta points ..................
       DO ITHETA = 1, NTHP
        CALL DVECZ( FTF, LENV )
        ZCOEF  = ZERO
        INDLOC = -1
        DO IPHI = 1, NPHP
          FAC2   = SF( ICM, IPHI, ITHETA, IR )
          INDLOC = INDLOC + 2
          FTF( INDLOC ) = FAC2
          ZCOEF = ZCOEF + FAC2
        ENDDO
C       .
C       . Do Fast Fourier Transform
C       .
        CALL FFTRLV( FTF, NPHP, ISIGN )
C       .
C       . ftf now contains spectral coefficients
C       . Loop around harmonics
C       .
        INDLOC = IRMNR
        DO IH3 = 1, NH3
           INDLOC = INDLOC + NR
C
           L      = ML3( IH3 )
C
C First do monopole term
C
           IF ( L.EQ.0 ) THEN
             FAC2 = SF2SA( IH3, ITHETA )*ZCOEF
             SV3( INDLOC ) = SV3( INDLOC ) + FAC2
             GOTO 70
           ENDIF
C          .
C          . Find "pseudo M"
C          .
           M      = MM3( IH3 )
C          .
           IF ( M.GE.0 ) THEN
             INDEX = 1 + 2*M/M0
           ELSE
             INDEX = 2 - 2*M/M0
           ENDIF
C          .
           FAC2 = SF2SA( IH3, ITHETA )*FTF( INDEX )
           SV3( INDLOC ) = SV3( INDLOC ) + FAC2
C          .
 70     CONTINUE
        ENDDO
C       .      Ended loop ih3 = 1, nh3
C       .
       ENDDO
C ............. ended looping around theta points ...................
C     .
      ENDDO
C ............. ended looping around radial grid nodes ..............
C     .
      RETURN
      END
C*********************************************************************
