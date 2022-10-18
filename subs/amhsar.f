C*********************************************************************
C subroutine Adapted Matrix Heat Source Auxilliary Routine ***********
C            -       -      -    -      -          -       ***********
C Steve Gibbons Thu Nov 18 17:54:58 GMT 1999                         C
C____________________________________________________________________C
C                                                                    C
C If \nabla^2 T_0( r ) = C (C is a constant), then T_0 has the       C
C general solution                                                   C
C                                                                    C
C            CB1 r^2                                                 C
C T_0( r ) = ------- - CB2 r^{-1}      (cb1 and cb2 are constants)   C
C               2                                                    C
C                                                                    C
C                                                                    C
C and so                                                             C
C                                                                    C
C              [                                     ]               C
C \nabla T_0 = [  CB1 r  +  CB2 r^{-2}     0     0   ]               C
C              [                       ,      ,      ]               C
C                                                                    C
C                                                                    C
C                                          [         CB2  ]          C
C and so v . \nabla T_0 = l(l+1)P(r) Y_l^m [ CB1 +  ----- ]          C
C                                          [         r^3  ]          C
C                                                                    C
C AMHSAR is an auxilliary routine to AMLICA for adding these terms.  C
C CB1 must be stored in DPARS( 1 )                                   C
C CB2 must be stored in DPARS( 2 )                                   C
C IPARS( 1 ) must contain L (degree of velocity harmonic).           C
C All other terms are looked after by AMLICA.                        C
C____________________________________________________________________C
C
C*********************************************************************
      SUBROUTINE AMHSAR( CVEC, RAD, IPARS, DPARS, IHD )
C____________________________________________________________________C
C Variable declarations - Parameters ................................C
      INTEGER IPARS( * ), IHD
      DOUBLE PRECISION CVEC( * ), RAD, DPARS( * )
C____________________________________________________________________C
C Variable declarations - Working variables .........................C
      INTEGER L, ND
      DOUBLE PRECISION TOL, CB1, CB2, R3
      PARAMETER ( TOL = 1.0d-8 )
C____________________________________________________________________C
C START OF PROGRAM **************************************************C
C____________________________________________________________________C
C Check on value of RAD
C     .
      IF ( ABS( RAD ).LT.TOL ) THEN
         PRINT *,' Subroutine AMHSAR.'
         PRINT *,' RAD = ', RAD
         PRINT *,' Program aborted.'
         STOP
      ENDIF
C     .
      R3 = RAD*RAD*RAD
C     .
      L     = IPARS( 1 )
C     .
      CB1   = DPARS( 1 )
      CB2   = DPARS( 2 )
C     .
C     . ok - so our chosen boundary condition is valid
C     . now check that IHD is large enough
C     .
      IF (     IHD.LT.0    ) THEN
         PRINT *,' Subroutine AMHSAR.'
         PRINT *,' IHD = ',IHD,' and must be atleast 0.'
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
      CVEC( 1 ) = DBLE( L*L + L )*( CB1 + CB2/R3 )
C     .
      RETURN
      END
C*********************************************************************

