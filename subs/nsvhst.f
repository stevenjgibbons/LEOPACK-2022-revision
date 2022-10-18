C*********************************************************************
C subroutine New Solution Vector Heat Source Term ********************
C            -   -        -      -    -      -    ********************
C Steve Gibbons Thu Nov  9 15:36:15 MET 2000                         C
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
C NSVHST adds this term to the appropriate temperature radial funcs. C
C                                                                    C
C The poloidal velocity terms are stored in the vector SV1.          C
C The radial function IH, at grid node IR is stored at the location  C
C    ind = ( ih - 1 )*nr + ir                                        C
C SV1 must contain ONLY poloidal velocity scalars and the l and m    C
C indices are stored in the arrays ML1 and MM1. There are NH1 pol.   C
C vel. harmonics.                                                    C
C                                                                    C
C The vector FT3 contains the forcing terms to the temperature       C
C equation and is stored in the same way. There are NH3 radial       C
C functions with l and m indices given by ML3 and MM3.               C
C                                                                    C
C The terms are added for nodes IR = 2 to IR = nr - 1.               C
C The solved vector SV1 contains only values for nodes 3 to NR - 2   C
C and so the routine PVVCPL must be called in advance to complete    C
C the poloidal velocity solution vector.                             C
C                                                                    C
C____________________________________________________________________C
C                                                                    C
C Input variables :-                                                 C
C ===============                                                    C
C  Integer                                                           C
C  -------                                                           C
C                                                                    C
C     NR        : Number of radial grid nodes.                       C
C                  See above for key. (corresponds to input vec.)    C
C     NH1       : Number of poloidal velocity radial functions.      C
C     ML1       : Array length ( * ) - atleast length NH1            C
C                  Sph. harm. degree, l, for poloidal velocity.      C
C     MM1       : Array length ( * ) - atleast length NH1            C
C                  Sph. harm. order, m, for cos m phi dep. (p.vel.)  C
C                 -Sph. harm. order, m, for sin m phi dep. (p.vel.)  C
C                                                                    C
C     NH3       : Number of temperature radial functions.            C
C     ML3       : Array length ( * ) - atleast length NH3            C
C                  Sph. harm. degree, l, for temperature.            C
C     MM3       : Array length ( * ) - atleast length NH3            C
C                  Sph. harm. order, m, for cos m phi dep. (temp.)   C
C                 -Sph. harm. order, m, for sin m phi dep. (temp.)   C
C                                                                    C
C  Double Precision                                                  C
C  ----------------                                                  C
C                                                                    C
C     SV1       : Solution vector. Dim ( * ) ( poloidal velocity )   C
C                 Length must be atleast NR*NH1                      C
C     FT3       : Forcing term for heat eqn. Dim ( * )               C
C                 Length must be atleast NR*NH3                      C
C     FAC       : Multiplier of term to be added.                    C
C     XARR      : Dim ( NR ). i^{th} element gives x_i.              C
C                                                                    C
C     CB1       : Constant - see above.                              C
C     CB2       : Constant - see above.                              C
C                                                                    C
C____________________________________________________________________C
C
C*********************************************************************
      SUBROUTINE NSVHST( NR, NH1, ML1, MM1, NH3, ML3, MM3, SV1, FT3,
     1                   FAC, XARR, CB1, CB2 )
      IMPLICIT NONE
C____________________________________________________________________C
C Variable declarations - Parameters ................................C
      INTEGER NR, NH1, ML1(NH1), MM1(NH1), NH3, ML3(NH3), MM3(NH3)
      DOUBLE PRECISION SV1( * ), FT3( * ), FAC, XARR( NR ), CB1, CB2
C____________________________________________________________________C
C Variable declarations - Working variables .........................C
C
      INTEGER IR, L1, MMM1, IH1, IH3, IND1, IND3, IBEG1, IBEG3, INDL
      DOUBLE PRECISION DLOW, RAD, D0F, F1, R3, ZERO, DLL1
      PARAMETER ( ZERO = 0.0d0, DLOW = 1.0d-9 )
C____________________________________________________________________C
C START OF PROGRAM **************************************************C
C____________________________________________________________________C
C
C     . early exit?
C
      IF ( FAC.EQ.ZERO ) RETURN
C     .
C     . Loop around poloidal velocity harmonics
C     .
      DO IH1 = 1, NH1
        L1   = ML1( IH1 )
        MMM1 = MM1( IH1 )
        INDL = L1*L1 + L1
        DLL1 = DBLE( INDL )
C
C so let's look for the corresponding temperature harmonic in FT3
C
        DO IH3 = 1, NH3
          IF (   ML3( IH3 ).EQ.L1 .AND. MM3( IH3 ).EQ.MMM1  ) THEN
C
C o.k. we've found the corresponding harmonic
C
            IBEG1 = (IH1 - 1)*NR
            IBEG3 = (IH3 - 1)*NR
            DO IR = 2, NR - 1
C
C Find locations of source and destination vectors ...
C
              IND1 = IBEG1 + IR
              IND3 = IBEG3 + IR
C
              RAD = XARR( IR )
              IF ( ABS( RAD ).LT.DLOW ) THEN
                PRINT *,' Subroutine NSVHST.'
                PRINT *,' Rad at node ',IR,' is ',RAD
                PRINT *,' Division by zero imminent.'
                PRINT *,' Program aborted.'
                STOP
              ENDIF
              R3  = RAD*RAD*RAD
              D0F = SV1( IND1 )
              F1  = D0F*DLL1
              FT3( IND3 ) = FT3( IND3 ) + FAC*( CB1 + CB2/R3 )*F1
C
            ENDDO
C           . End of loop ir = 2, nr-1
C
C Jump out of ih3 loop if we have found our harmonic
C
            GOTO 46
C
          ENDIF
C         . check ml3(ih3).eq.l1 and mm3(ih3).eq.mm1
        ENDDO
C       . End of loop ih3 = 1, nh3
C
 46   CONTINUE
      ENDDO
C     . End of loop ih1 = 1, nh1
C
      RETURN
      END
C*********************************************************************
