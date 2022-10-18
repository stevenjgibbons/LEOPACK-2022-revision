C*********************************************************************
C subroutine Xtra Special V Cross B Evaluate *************************
C            -    -       - -     - -        *************************
C Steve Gibbons Mon Dec 18 10:28:51 WET 2000                         C
C____________________________________________________________________C
C                                                                    C
C We have received an array XSV( NCMX, NPHP, NTHP, NR ) into which   C
C we have placed the following:-                                     C
C                                                                    C
C  ICM             Function                                          C
C  ---             --------                                          C
C                                                                    C
C    ICMRAD        v_r           (6.lt.ICMRAD .and. ICMRAD.le.NCMX)  C
C    ICMTHE        v_{theta}     (6.lt.ICMTHE .and. ICMTHE.le.NCMX)  C
C    ICMPHI        v_{phi}       (6.lt.ICMPHI .and. ICMPHI.le.NCMX)  C
C                                                                    C
C (Also, ofcourse ICMRAD, ICMTHE and ICMPHI are distinct) ...        C
C                                                                    C
C    4             B_r                                               C
C    5             B_{theta}                                         C
C    6             B_{phi}                                           C
C                                                                    C
C Into ICM = 1, 2 and 3, we place the r, theta and phi components    C
C respectively of ( v x B ).                                         C
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
C     ICMRAD    : Component of XSV containing the velocity rad comp. C
C     ICMTHE    : Component of XSV containing the velocity the comp. C
C     ICMPHI    : Component of XSV containing the velocity phi comp. C
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
C____________________________________________________________________C
C
C*********************************************************************
      SUBROUTINE XSVCBE( NTHP, NPHP, NR, ILNR, IRNR, NCMX, ICMRAD,
     1                   ICMTHE, ICMPHI, XSV )
      IMPLICIT NONE
C____________________________________________________________________C
C Variable declarations - Parameters ................................C
      INTEGER          NTHP, NPHP, NR, ILNR, IRNR, NCMX, ICMRAD,
     1                 ICMTHE, ICMPHI
      DOUBLE PRECISION XSV( NCMX, NPHP, NTHP, NR )
C____________________________________________________________________C
C Variable Declarations - Working Variables .........................C
      INTEGER          ITHE, IPHI, IR
      DOUBLE PRECISION VRAD, VTHE, VPHI, BRAD, BTHE, BPHI
C____________________________________________________________________C
C START OF PROGRAM **************************************************C
C____________________________________________________________________C
C
C Check value of NCMX
C
      IF ( NCMX.LT.9  ) THEN
        PRINT *,' Subroutine XSVCBE.'
        PRINT *,' NCMX = ', NCMX
        PRINT *,' Program aborted.'
        STOP
      ENDIF
C
C Check values of ICMRAD etc.
C
      IF ( ICMRAD.LE.6 .OR. ICMRAD.GT.NCMX  .OR. ICMTHE.LE.6 .OR.
     1     ICMTHE.GT.NCMX .OR. ICMPHI.LE.6 .OR. ICMPHI.GT.NCMX .OR.
     2     ICMRAD.EQ.ICMTHE .OR. ICMTHE.EQ.ICMPHI .OR.
     3     ICMPHI.EQ.ICMRAD  ) THEN
        PRINT *,' Subroutine XSVCBE.'
        PRINT *,' ICMRAD = ', ICMRAD
        PRINT *,' ICMTHE = ', ICMTHE
        PRINT *,' ICMPHI = ', ICMPHI
        PRINT *,' NCMX   = ', NCMX
        PRINT *,' Program aborted.'
        STOP
      ENDIF
C
C Loop around radial grid nodes
C
      DO IR = ILNR, IRNR
C      ............. loop around theta and phi points ................
       DO ITHE = 1, NTHP
         DO IPHI = 1, NPHP
C          .
C          . Fill all the scalar variables
C          . First, the velocity components
C          .
           VRAD   = XSV( ICMRAD, IPHI, ITHE, IR )
           VTHE   = XSV( ICMTHE, IPHI, ITHE, IR )
           VPHI   = XSV( ICMPHI, IPHI, ITHE, IR )
C          .
C          . Now B components
C          .
           BRAD   = XSV( 4, IPHI, ITHE, IR )
           BTHE   = XSV( 5, IPHI, ITHE, IR )
           BPHI   = XSV( 6, IPHI, ITHE, IR )
C          .
C          . Now evaluate forcing terms for
C          . induction equation.
C          .
           XSV( 1, IPHI, ITHE, IR ) = VTHE*BPHI - VPHI*BTHE
C
           XSV( 2, IPHI, ITHE, IR ) = VPHI*BRAD - VRAD*BPHI
C
           XSV( 3, IPHI, ITHE, IR ) = VRAD*BTHE - VTHE*BRAD
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

