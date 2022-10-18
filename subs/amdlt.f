C*********************************************************************
C subroutine Adapted Matrix DL Terms *********************************
C            -       -      -- -     *********************************
C Steve Gibbons Tue Oct 26 15:24:17 GMT 1999                         C
C____________________________________________________________________C
C                                                                    C
C Supplies coefficients to the routine AMLICA for terms with the     C
C function D_l to the power 0, 1 or 2. Must be declared EXTERNAL     C
C in the (sub)program which calls AMLICA and must replave SUB1 in    C
C the calling sequence.                                              C
C                                                                    C
C Set IPARS( 1 ) to IFLAG                                            C
C                                                                    C
C Now if IFLAG = 1, D_l^0 is added to the matrix - i.e. the identity C
C                   IHD can be zero.                                 C
C                                                                    C
C        IFLAG = 2, D_l is added to the matrix.                      C
C                   This returns                                     C
C                                                                    C
C                     CVEC( 1 ) = (-1.0d0)*DBLE( L*L+L )/(RAD*RAD)   C
C                     CVEC( 2 ) =   2.0d0/RAD                        C
C                     CVEC( 3 ) =   1.0d0                            C
C                                                                    C
C                   IHD must be atleast 2.                           C
C                                                                    C
C        IFLAG = 3, D_l^2 is added to the matrix                     C
C                                                                    C
C                     CVEC( 1 ) = (L+2)(L+1)L(L-1)/RAD**4            C
C                     CVEC( 2 ) =  0.0d0                             C
C                     CVEC( 3 ) = -2L(L+1)/RAD**2                    C
C                     CVEC( 4 ) = 4/RAD                              C
C                     CVEC( 5 ) =  1.0d0                             C
C                                                                    C
C  The integer L is passed in as IPARS( 2 ).                         C
C                                                                    C
C  All other parameters are looked after by AMLICA and most are      C
C  dummy parameters.                                                 C
C                                                                    C
C____________________________________________________________________C
C
C*********************************************************************
      SUBROUTINE AMDLT( CVEC, RAD, IPARS, DPARS, IHD )
      IMPLICIT NONE
C____________________________________________________________________C
C Variable declarations - Parameters ................................C
      INTEGER IPARS( * ), IHD
      DOUBLE PRECISION CVEC( * ), RAD, DPARS( * )
C____________________________________________________________________C
C Variable declarations - Working variables .........................C
      INTEGER IHDMIN( 3 ), IFLAG, L, ND
      DOUBLE PRECISION TOL
      PARAMETER ( TOL = 1.0d-8 )
C____________________________________________________________________C
C START OF PROGRAM **************************************************C
C____________________________________________________________________C
C Check on value of RAD
C     .
      IF ( ABS( RAD ).LT.TOL ) THEN
         PRINT *,' Subroutine AMDLT.'
         PRINT *,' RAD = ', RAD
         PRINT *,' Program aborted.'
         STOP
      ENDIF
C     .
C
C Just apply the minimum values of IHD to the IHDMIN array.
C
C     .
      IHDMIN( 1 ) = 0
      IHDMIN( 2 ) = 2
      IHDMIN( 3 ) = 4
C     .
      IFLAG = IPARS( 1 )
      L     = IPARS( 2 )
C     .
C     . Check for valid value of IFLAG
C     .
      IF ( IFLAG.LT.1 .OR. IFLAG.GT.3 ) THEN
         PRINT *,' Subroutine AMDLT.'
         PRINT *,' IFLAG = ', IFLAG
         PRINT *,' Program aborted.'
         STOP
      ENDIF
C     .
C     . ok - so our chosen boundary condition is valid
C     . now check that IHD is large enough
C     .
      IF (     IHD.NE.IHDMIN( IFLAG )    ) THEN
         PRINT *,' Subroutine AMDLT.'
         PRINT *,' IHD = ',IHD,' and must be atleast'
         PRINT *,  IHDMIN( IFLAG ),' for option ',IFLAG
         PRINT *,' Program aborted.'
         STOP
      ENDIF
C     .
C     . Zero all coefficients up to IHD
C     .
      DO ND = 1, IHD + 1
        CVEC( ND ) = 0.0d0
      ENDDO
C     .
C     . Option 1:
C     . D_l^0
C     .
      IF ( IFLAG.EQ.1 ) THEN
         CVEC( 1 ) = 1.0d0
         RETURN
      ENDIF
C     .
C     . Option 2:
C     . D_l
C     .
      IF ( IFLAG.EQ.2 ) THEN
         CVEC( 1 ) = (-1.0d0)*DBLE( L*L+L )/(RAD*RAD)
         CVEC( 2 ) = 2.0d0/RAD
         CVEC( 3 ) = 1.0d0
         RETURN
      ENDIF
C     .
C     . Option 3:
C     . D_l^2
C     .
      IF ( IFLAG.EQ.3 ) THEN
         CVEC( 1 ) = DBLE( (L+2)*(L+1)*L*(L-1) )/RAD**4
         CVEC( 3 ) = (-2.0d0)*DBLE( L*L + L )/RAD**2
         CVEC( 4 ) = 4.0d0/RAD
         CVEC( 5 ) = 1.0d0
         RETURN
      ENDIF
C     .
      END
C*********************************************************************

