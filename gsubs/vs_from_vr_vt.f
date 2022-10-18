C*********************************************************************
C subroutine VS FROM VR and V_Theta **********************************
C            -- ---- --     - -     **********************************
C Steve Gibbons Wed Jan 24 11:31:11 WET 2001                         C
C____________________________________________________________________C
C                                                                    C
C When doing a constant z plot of velocities (for arrows etc.)       C
C it is the s component of the flow which is important. This routine C
C takes real arrays of dimension (NS, NPHI) for VR and VTheta and    C
C fills the array with VS with vs = vr sin theta + v_theta cos theta C
C                                                                    C
C The height above the equator is given by the double precision ZED. C
C                                                                    C
C The coordinates are defined by a COMMON block:                     C
C                                                                    C
C      COMMON  / PARAMP / NRAD, NTHE, RFIRST, RLAST, TFIRST, TLAST   C
C                                                                    C
C It consists of two INTEGERs                                        C
C NRAD (number of equally spaced radial grid nodes) and              C
C NTHE (number of equally spaced phi grid nodes), and four REAL      C
C variables RFIRST, RLAST, TFIRST, TLAST                             C
C                                                                    C
C RFIRST is the inner displacement from axis.                        C
C RLAST is the outer displacement from axis.                         C
C                                                                    C
C (Since we are doing constant z, these may not actually             C
C correspond with the inner and outer radii (they only will in the   C
C case z = 0 infact).                                                C
C                                                                    C
C TFIRST is the first phi point, in radians.                         C
C TLAST is the last phi point, in radians.                           C
C____________________________________________________________________C
C                                                                    C
C Input variables :-                                                 C
C ===============                                                    C
C                                                                    C
C  Integer                                                           C
C  -------                                                           C
C                                                                    C
C     NS        : Number of s grid nodes.                            C
C     NPHI      : Number of phi grid nodes.                          C
C                                                                    C
C  Real                                                              C
C  ----                                                              C
C                                                                    C
C     VR        : Dim ( NS, NPHI ). Radial velocity component        C
C     VT        : Dim ( NS, NPHI ). Theta velocity component         C
C                                                                    C
C     VS        : Dim ( NS, NPHI ). S velocity component (output)    C
C                                                                    C
C  Double Precision                                                  C
C  ----------------                                                  C
C                                                                    C
C     ZED       : Height above equator.                              C
C                                                                    C
C____________________________________________________________________C
C
C*********************************************************************
      SUBROUTINE VS_FROM_VR_VT( NS, NPHI, VR, VT, VS, ZED )
      IMPLICIT NONE
C____________________________________________________________________C
C Common block variables ............................................C
      INTEGER NRAD, NTHE
      REAL    RFIRST, RLAST, TFIRST, TLAST
      COMMON  / PARAMP / NRAD, NTHE, RFIRST, RLAST, TFIRST, TLAST
C____________________________________________________________________C
C Variable declarations - Parameters ................................C
      INTEGER NS, NPHI
      REAL             VR( NRAD, NTHE ), VT( NRAD, NTHE ),
     1                 VS( NRAD, NTHE )
      DOUBLE PRECISION ZED
C____________________________________________________________________C
C Variable declarations - Working variables .........................C
      INTEGER          IS, IPHI
      REAL             REALZ, R1, THE, COSTH, SINTH
C____________________________________________________________________C
C START OF PROGRAM **************************************************C
C____________________________________________________________________C
C
      IF ( NS.NE.NRAD .OR. NPHI.NE.NTHE ) THEN
         PRINT *,' Subroutine VS_FROM_VR_VT.'
         PRINT *,' NS   = ', NS  ,' NRAD = ', NRAD
         PRINT *,' NPHI = ', NPHI,' NTHE = ', NTHE
         PRINT *,' Program aborted.'
         STOP
      ENDIF
C
      REALZ = REAL( ZED )
C
      DO IS = 1, NS
        R1    = RFIRST + (RLAST-RFIRST)*REAL(IS-1)/REAL(NS-1)
        THE   = ATAN2( R1, REALZ )
        COSTH = COS( THE )
        SINTH = SIN( THE )
        DO IPHI = 1, NPHI
          VS( IS, IPHI ) =
     1                       VR( IS, IPHI )*SINTH +
     2                       VT( IS, IPHI )*COSTH
        ENDDO
      ENDDO
C
      RETURN
      END
C*********************************************************************
