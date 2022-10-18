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
C then SF2SDX will add FAC*SF to the appropriate elements in SV3.    C
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
C     LH        : Maximum spherical harmonic degree, l.              C
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
C     PA        : Schmidt Normalised Legendre Functions              C
C                  Dim { ( LH + 1 )*( LH + 2 )/2 , NTHP }            C
C                  P_l^m(X_j) = PA( L*(L+1)/2 + M + 1 , j )          C
C                                                                    C
C     GAUW      : Gauss weights as calculated by GAUWTS. ( NTHP )    C
C                                                                    C
C     XSV       : Xtra Special Function: is an array containing a    C
C                  function over a set of theta, phi and r points.   C
C                    Dimensions are                                  C
C                      ( NCMX, NPHP, NTHP, NR )                      C
C                                                                    C
C     SV3       : Dim ( NR*NH3 ) Radial functions psi_{ip3}(r)       C
C                                                                    C
C     FAC       : Coefficient of function to be added.               C
C                                                                    C
C     FTF       : Work array: dim ( 2*NPHP )                         C
C                                                                    C
C____________________________________________________________________C
C
C*********************************************************************
      SUBROUTINE SF2SDX( ILNR, IRNR, NR, LH, NH3, ML3, MM3, M0, NTHP,
     1                  NPHP, NCMX, ICM, PA, GAUW, SF, SV3, FAC, FTF )
      IMPLICIT NONE
C____________________________________________________________________C
C Variable declarations - Parameters ................................C
      INTEGER          ILNR, IRNR, NR, LH, NH3, ML3( NH3 ),
     1                 MM3( NH3 ), M0, NTHP, NPHP, NCMX, ICM
      DOUBLE PRECISION PA( ( LH + 1 )*( LH + 2 )/2, NTHP ),
     1                 GAUW( NTHP ), SF( NCMX, NPHP, NTHP, NR ),
     2                 FTF( 2*NPHP ), SV3( * ), FAC
C____________________________________________________________________C
C Variable declarations - Working variables .........................C
      INTEGER          ISIGN, ITHETA, IOP, LENV, IPHI, INDLOC, IH3,
     1                 ICS, L, M, INDEX, IP, MPS, IR
      DOUBLE PRECISION ZERO, ZCOEF, FAC2, FACL, DALF, WEIGHT
      PARAMETER        ( ZERO = 0.0d0, IOP = 0, ISIGN = 1 )
C____________________________________________________________________C
C START OF PROGRAM **************************************************C
C____________________________________________________________________C
C
C Check value of ICM
C
      IF ( ICM.LT.1 .OR. ICM.GT.NCMX ) THEN
        PRINT *,' Subroutine SF2SDX.'
        PRINT *,' ICM = ', ICM
        PRINT *,' Program aborted.'
        STOP
      ENDIF
C
C ............... begin looping around radial grid nodes .............
      DO IR = ILNR, IRNR
C
C ............... begin looping around theta points ..................
       DO ITHETA = 1, NTHP
        LENV = 2*NPHP
        CALL VECOP( FTF, ZERO, LENV, IOP )
        ZCOEF  = ZERO
        WEIGHT = GAUW( ITHETA )
        DO IPHI = 1, NPHP
          FAC2   = SF( ICM, IPHI, ITHETA, IR )
          INDLOC = 2*IPHI - 1
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
        DO IH3 = 1, NH3
           INDLOC = IH3*NR - NR + IR
C
           L      = ML3( IH3 )
C
C First do monopole term
C
           IF ( L.EQ.0 ) THEN
             FAC2 = 0.5d0*FAC*ZCOEF*WEIGHT/DBLE( NPHP )
             SV3( INDLOC ) = SV3( INDLOC ) + FAC2
             GOTO 70
           ENDIF
C
           FACL   = 0.25d0*DBLE( 2*L + 1 )
           M      = MM3( IH3 )
           ICS    = 1
C          .
C          . Modify M and ICS if the harmonic is sin (m phi) dep.
C          .
           IF ( M.LT.0   ) ICS = 2
           IF ( ICS.EQ.2 ) M   = -M
C          .
C          . Store associated Legendre Functions
C          .
           IP     = L*(L+1)/2+M+1
           DALF   = PA( IP, ITHETA )
C          .
C          . Make sure that harmonic is
C          . compatible with the symmetry imposition
C          . MPS (pseudo M) = M/M0
C          .
           MPS = M/M0
C Uncomment following lines for this safety check:
c          IF ( MPS*M0.NE.M ) THEN
c            PRINT *,' Subroutine SF2SDX.'
c            PRINT *,' MM3(',IH3,') = ', MM3( IH3 )
c            PRINT *,' Program aborted.'
c            STOP
c          ENDIF
C          .
           INDEX          = 2*MPS + ICS
           FAC2 = FAC*WEIGHT*DALF*FTF( INDEX )*FACL
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
