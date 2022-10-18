C*********************************************************************
C subroutine Magnetic Dynamo and Convection Xtra Special Evaluation **
C            -        -          -          -    -       -          **
C Steve Gibbons Tue Dec 12 10:19:32 WET 2000                         C
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
C    4             (Grad Theta)_r                                    C
C    5             (Grad Theta)_{theta}                              C
C    6             (Grad Theta)_{phi}                                C
C    7             (curl v)_r                                        C
C    8             (curl v)_{theta}                                  C
C    9             (curl v)_{phi}                                    C
C    10            B_r                                               C
C    11            B_{theta}                                         C
C    12            B_{phi}                                           C
C    13            (curl B)_r                                        C
C    14            (curl B)_{theta}                                  C
C    15            (curl B)_{phi}                                    C
C                                                                    C
C  Into ICM  = 16, we put v.Grad( theta ) and we put the r, theta    C
C and phi components of   (-CG)*( k x V ) + CF*( V x curl V )        C
C                         + CJ*( curl B x B )                        C
C                    respectively into ICM  = 17, 18 and 19.         C
C r, theta and phi components of ( v x B ) are put into              C
C  ICM = 20, 21 and 22.                                              C
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
C     CJ        : Coefficient of ( B x curl B ) term.                C
C                                                                    C
C____________________________________________________________________C
C
C*********************************************************************
      SUBROUTINE MDCXSE( NTHP, NPHP, NR, ILNR, IRNR, NCMX, XSV, GAUX,
     1                   CG, CF, CJ )
      IMPLICIT NONE
C____________________________________________________________________C
C Variable declarations - Parameters ................................C
      INTEGER          NTHP, NPHP, NR, ILNR, IRNR, NCMX
      DOUBLE PRECISION XSV( NCMX, NPHP, NTHP, NR ), GAUX( NTHP ),
     1                 CG, CF, CJ
C____________________________________________________________________C
C Variable Declarations - Working Variables .........................C
      INTEGER          ITHE, IPHI, IR
      DOUBLE PRECISION VRAD, VTHE, VPHI, GRAD, GTHE, GPHI,
     1                 CVRAD, CVTHE, CVPHI, SCALF, OUTPHI, OUTRAD,
     2                 OUTTHE, COSTH, SINTH, BRAD, BTHE, BPHI,
     3                 CBRAD, CBTHE, CBPHI, COSTHG, SINTHG, COSTHM
C____________________________________________________________________C
C START OF PROGRAM **************************************************C
C____________________________________________________________________C
C
C Check value of NCMX
C
      IF ( NCMX.LT.22 ) THEN
        PRINT *,' Subroutine MDCXSE.'
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
         COSTHG = COSTH*CG
         COSTHM = COSTHG*(-1.0d0)
         SINTHG = SINTH*CG
         DO IPHI = 1, NPHP
C          .
C          . Fill all the scalar variables
C          . First, the velocity components
C          .
           VRAD   = XSV(  1, IPHI, ITHE, IR )
           VTHE   = XSV(  2, IPHI, ITHE, IR )
           VPHI   = XSV(  3, IPHI, ITHE, IR )
C          .
C          . The Grad( theta ) components
C          .
           GRAD   = XSV(  4, IPHI, ITHE, IR )
           GTHE   = XSV(  5, IPHI, ITHE, IR )
           GPHI   = XSV(  6, IPHI, ITHE, IR )
C          .
C          . Now curl v components
C          .
           CVRAD  = XSV(  7, IPHI, ITHE, IR )*CF
           CVTHE  = XSV(  8, IPHI, ITHE, IR )*CF
           CVPHI  = XSV(  9, IPHI, ITHE, IR )*CF
C          .
C          . Now B components
C          .
           BRAD   = XSV( 10, IPHI, ITHE, IR )
           BTHE   = XSV( 11, IPHI, ITHE, IR )
           BPHI   = XSV( 12, IPHI, ITHE, IR )
C          .
C          . Now curl B components
C          .
           CBRAD  = XSV( 13, IPHI, ITHE, IR )*CJ
           CBTHE  = XSV( 14, IPHI, ITHE, IR )*CJ
           CBPHI  = XSV( 15, IPHI, ITHE, IR )*CJ
C          .
           SCALF  = VRAD*GRAD + VTHE*GTHE + VPHI*GPHI
           XSV( 16, IPHI, ITHE, IR ) = SCALF
C          .
C          . Evaluate Coriolis terms
C          .
           OUTRAD = SINTHG*VPHI
           OUTTHE = COSTHG*VPHI
           OUTPHI = COSTHM*VTHE - SINTHG*VRAD
C          .
C          . Now evaluate forcing terms for
C          . vorticity equation.
C          .
           XSV( 17, IPHI, ITHE, IR ) = OUTRAD +
     1           VTHE*CVPHI - VPHI*CVTHE
     2         - BTHE*CBPHI + BPHI*CBTHE
C
           XSV( 18, IPHI, ITHE, IR ) = OUTTHE +
     1           VPHI*CVRAD - VRAD*CVPHI
     2         - BPHI*CBRAD + BRAD*CBPHI
C
           XSV( 19, IPHI, ITHE, IR ) = OUTPHI +
     1           VRAD*CVTHE - VTHE*CVRAD
     2         - BRAD*CBTHE + BTHE*CBRAD
C          .
C          . Now evaluate forcing terms for
C          . induction equation.
C          .
           XSV( 20, IPHI, ITHE, IR ) = 
     1           VTHE*BPHI - VPHI*BTHE
C
           XSV( 21, IPHI, ITHE, IR ) = 
     1           VPHI*BRAD - VRAD*BPHI
C
           XSV( 22, IPHI, ITHE, IR ) = 
     1           VRAD*BTHE - VTHE*BRAD
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
