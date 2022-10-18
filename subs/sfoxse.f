C*********************************************************************
C subroutine Shear Flow Only Xtra Special Evaluation *****************
C            -     -    -    -    -       -          *****************
C Steve Gibbons Wed Dec 20 10:29:24 WET 2000                         C
C____________________________________________________________________C
C                                                                    C
C We have received an array XSV( NCMX, NPHP, NTHP, NR ) into which   C
C we have placed the following:-                                     C
C                                                                    C
C  ICM             Function                                          C
C  ---             --------                                          C
C                                                                    C
C    1             v_r                                               C
C    2             v_{theta}                                         C
C    3             v_{phi}                                           C
C    4             (curl v)_r                                        C
C    5             (curl v)_{theta}                                  C
C    6             (curl v)_{phi}                                    C
C                                                                    C
C we put the r, theta and phi components of                          C
C              (-CG)*( k x V ) + CF*( V x curl V )                   C
C                    respectively into ICM  =  7,  8 and 9.          C
C                                                                    C
C                                                                    C
C____________________________________________________________________C
C                                                                    C
C Input variables :-                                                 C
C ===============                                                    C
C  Integer                                                           C
C  -------                                                           C
C     NTHP      : Number of theta points.                            C
C     NPHP      : Number of phi points.                              C
C     NR        : Number of radial grid nodes.                       C
C     ILNR      : Lowest radial grid node.                           C
C     IRNR      : Highest radial grid node.                          C
C     NCMX      : Maximum number of components stored in XSV         C
C                                                                    C
C  Double Precision                                                  C
C  ----------------                                                  C
C                                                                    C
C     XSV       : eXtra Special Vector. An array of dimensions       C
C                ( NCMX, NPHP, NTHP, NR)                             C
C                 At the radial grid node, IR, and theta point       C
C                 ITHE and phi point IPHI, the components listed     C
C                 above are stored in                                C
C                   XSV( ICM, IPHI, ITHE, IR )                       C
C                 where ICM corresponds to the number above.         C
C                                                                    C
C     GAUX      : Cosines of the NTHPTS evaluated by the routine     C
C                  gauwts. Dimension ( NTHPTS ).                     C
C     CG        : Coefficient of Coriolis term.                      C
C     CF        : Coefficient of ( V x curl V ) term.                C
C                                                                    C
C____________________________________________________________________C
C
C*********************************************************************
      SUBROUTINE SFOXSE( NTHP, NPHP, NR, ILNR, IRNR, NCMX, XSV, GAUX,
     1                   CG, CF )
      IMPLICIT NONE
C____________________________________________________________________C
C Variable declarations - Parameters ................................C
      INTEGER          NTHP, NPHP, NR, ILNR, IRNR, NCMX
      DOUBLE PRECISION XSV( NCMX, NPHP, NTHP, NR ), GAUX( NTHP ),
     1                 CG, CF
C____________________________________________________________________C
C Variable Declarations - Working Variables .........................C
      INTEGER          ITHE, IPHI, IR
      DOUBLE PRECISION VRAD, VTHE, VPHI, CRAD, CTHE, CPHI, OUTPHI,
     1                 OUTRAD, OUTTHE, COSTH, SINTH, COSTHG, SINTHG,
     2                 COSTHM
C____________________________________________________________________C
C START OF PROGRAM **************************************************C
C____________________________________________________________________C
C
C Check value of NCMX
C
      IF ( NCMX.LT.9 ) THEN
        PRINT *,' Subroutine SFOXSE.'
        PRINT *,' NCMX = ', NCMX
        PRINT *,' Program aborted.'
        STOP
      ENDIF
C
C Loop around radial grid nodes
C
      DO IR = ILNR, IRNR
C      ............. loop around theta and phi points ................
       DO ITHE = 1, NTHP
         COSTH  = GAUX( ITHE )
         SINTH  = DSQRT( 1.0d0 - COSTH*COSTH )
         COSTHG = CG*COSTH
         COSTHM = COSTHG*(-1.0d0)
         SINTHG = CG*SINTH
         DO IPHI = 1, NPHP
C          .
C          . Fill all the scalar variables
C          . First, the velocity components
C          .
           VRAD   = XSV( 1, IPHI, ITHE, IR )
           VTHE   = XSV( 2, IPHI, ITHE, IR )
           VPHI   = XSV( 3, IPHI, ITHE, IR )
C          .
C          . Now curl v components
C          .
           CRAD   = XSV( 4, IPHI, ITHE, IR )*CF
           CTHE   = XSV( 5, IPHI, ITHE, IR )*CF
           CPHI   = XSV( 6, IPHI, ITHE, IR )*CF
C          .
C          . Evaluate Coriolis terms
C          .
           OUTRAD = SINTHG*VPHI
           OUTTHE = COSTHG*VPHI
           OUTPHI = COSTHM*VTHE - SINTHG*VRAD
C          .
C          . Finally evaluate forcing terms
C          .
           XSV(  7, IPHI, ITHE, IR ) = OUTRAD +
     1           VTHE*CPHI - VPHI*CTHE
C
           XSV(  8, IPHI, ITHE, IR ) = OUTTHE +
     1           VPHI*CRAD - VRAD*CPHI
C
           XSV(  9, IPHI, ITHE, IR ) = OUTPHI +
     1           VRAD*CTHE - VTHE*CRAD
C          .
         ENDDO
       ENDDO
C ............. ended looping around theta, phi points ...............
      ENDDO
C
C End loop around radial grid nodes
C
      RETURN
      END
C*********************************************************************
