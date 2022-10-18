C*********************************************************************
C subroutine Magneto Hydro Dynamic Boundary Condition Subroutine *****
C            -       -     -       -        -         -          *****
C Steve Gibbons Wed Oct  6 06:39:09 BST 1999                         C
C____________________________________________________________________C
C                                                                    C
C Allows the enforcement of some of the most common boundary cond.s  C
C to matrices by use of the NGMBCE routine.                          C
C                                                                    C
C In the (sub)program which calls the routine NGMBCE, MHDBCS must be C
C declared EXTERNAL and MHDBCS must replace SUB1 in the calling      C
C sequence.                                                          C
C                                                                    C
C The function of MHDBCS is determined by IFLAG, the value of        C
C INARR( 5 ). Current values of IFLAG are                            C
C                                                                    C
C  IFLAG = 1. Fixed value. Simply will fix the value of a function   C
C             to the corresponding right hand side value.            C
C             Examples are poloidal velocity at impenetrable bndry,  C
C             temperature when fixed, toroidal velocity with rigid   C
C             boundary cond.s                                        C
C                                                                    C
C  IFLAG = 2. Fixed derivative. Simply fix the first deriv. of a     C
C             function to the value in right hand side.              C
C             Examples are poloidal velocity at rigid boundary and   C
C             temperature when subject to a fixed heat flux boundary C
C             condition.                                             C
C                                                                    C
C  IFLAG = 3. Fixed second deriv. Fix 2nd deriv. of function to      C
C             RHS value. Example: poloidal velocity with stress free C
C             boundary.                                              C
C                                                                    C
C  IFLAG = 4. t'(r) = t(r)/r  or r t'(r) - t(r) = 0                  C
C             Example is toroidal velocity with stress free b.c.     C
C             A zero value fo r is checked for, in which case t(r)   C
C             is simply set to zero.                                 C
C                                                                    C
C  IFLAG = 5. p'(r) = l p(r)/r or r p'(r) - l p(r) = 0               C
C                                                                    C
C             This is the boundary condition for the poloidal        C
C             magnetic field at an insulating inner boundary.        C
C             If r is zero then there is no inner core and we        C
C             merely enforce the regularity condition, p(r) = 0.     C
C                                                                    C
C             For the toroidal field at insulating boundaries,       C
C             always use IFLAG = 1.                                  C
C                                                                    C
C             The spherical harmonic degree, l, is passed in as      C
C             the element INARR( 6 ).                                C
C                                                                    C
C  IFLAG = 6. r p'(r) + (l+1) p(r) = 0                               C
C                                                                    C
C             This is the boundary condition for poloidal            C
C             magnetic field at an insulating outer boundary.        C
C                                                                    C
C             The spherical harmonic degree, l, is passed in as      C
C             the element INARR( 6 ).                                C
C                                                                    C
C  The only input parameters which need discussing are INARR( 5 )    C
C  and INARR( 6 ), The others are all supplied by NGMBCE and in most C
C  cases are dummy parameters.                                       C
C                                                                    C
C____________________________________________________________________C
C
C*********************************************************************
      SUBROUTINE MHDBCS( CVEC, IRAD, RAD, NLN, NRN, NBN, FDCM,
     1                   NFDCM, NR, INARR, DPARR, VEC0, IHD, XARR )
      IMPLICIT NONE
C____________________________________________________________________C
C Variable declarations - Parameters ................................C
      INTEGER IRAD, NLN, NRN, NBN, NFDCM, NR, INARR( * ), IHD
      DOUBLE PRECISION CVEC( * ), RAD, DPARR( * ), VEC0( *),
     1                 XARR( * ), FDCM( NFDCM, NR, * )
C____________________________________________________________________C
C Variable declarations - Working variables .........................C
      INTEGER IHDMIN( 6 ), IFLAG, L, ND
C____________________________________________________________________C
C START OF PROGRAM **************************************************C
C____________________________________________________________________C
C
C Just apply the minimum values of IHD to the IHDMIN array.
C
C     .
      IHDMIN( 1 ) = 0
      IHDMIN( 2 ) = 1
      IHDMIN( 3 ) = 2
      IHDMIN( 4 ) = 1
      IHDMIN( 5 ) = 1
      IHDMIN( 6 ) = 1
C     .
      IFLAG = INARR( 5 )
      L     = INARR( 6 )
C     .
C     . Check for valid value of IFLAG
C     .
      IF ( IFLAG.LT.1 .OR. IFLAG.GT.6 ) THEN
         PRINT *,' Subroutine MHDBCS.'
         PRINT *,' IFLAG = ', IFLAG
         PRINT *,' Program aborted.'
         STOP
      ENDIF
C     .
C     . ok - so our chosen boundary condition is valid
C     . now check that IHD is large enough
C     .
      IF (     IHD.LT.IHDMIN( IFLAG )    ) THEN
         PRINT *,' Subroutine MHDBCS.'
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
C     . Now go through the boundary conditions
C     .
C     . Option 1:
C     . Fixed value (zero^{th} derivative).
C     .
      IF ( IFLAG.EQ.1 ) THEN
         CVEC( 1 ) = 1.0d0
         RETURN
      ENDIF
C     .
C     . Option 2:
C     . Fixed first derivative.
C     .
      IF ( IFLAG.EQ.2 ) THEN
         CVEC( 2 ) = 1.0d0
         RETURN
      ENDIF
C     .
C     . Option 3:
C     . Fixed second derivative.
C     .
      IF ( IFLAG.EQ.3 ) THEN
         CVEC( 3 ) = 1.0d0
         RETURN
      ENDIF
C     .
C     . Option 4:
C     . r t'(r) - t(r) = 0
C     .
      IF ( IFLAG.EQ.4 ) THEN
         CVEC( 1 ) = -1.0d0
         CVEC( 2 ) = RAD
         RETURN
      ENDIF
C     .
C     . Option 5:
C     . r p'(r) - l p(r) = 0
C     .
      IF ( IFLAG.EQ.5 ) THEN
         CVEC( 1 ) = (-1.0d0)*DBLE( L )
         CVEC( 2 ) = RAD
         RETURN
      ENDIF
C     .
C     . Option 6:
C     . r p'(r) + (l+1) p(r) = 0
C     .
      IF ( IFLAG.EQ.6 ) THEN
         CVEC( 1 ) = DBLE( L ) + 1.0d0
         CVEC( 2 ) = RAD
         RETURN
      ENDIF
C     .
      END
C*********************************************************************
