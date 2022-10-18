C*********************************************************************
C subroutine New Solution Vector Buoyancy Term Add *******************
C            -   -        -      -        -    -   *******************
C Steve Gibbons Fri Nov 10 07:44:00 WET 2000                         C
C____________________________________________________________________C
C                                                                    C
C Adds the buoyancy term from SV3, a solution vector containing the  C
C temperature function, to the toroidal part of the vorticity        C
C equation in FT1. CH multiplies the terms added.                    C
C                                                                    C
C The temperature terms are stored in the vector SV3.                C
C The radial function IH, at grid node IR is stored at the location  C
C    ind = ( ih - 1 )*nr + ir                                        C
C SV3 must contain ONLY temperature scalars and the l and m indices  C
C are stored in the arrays ML3 and MM3. There are NH3 temp. harm.s   C
C                                                                    C
C The vector FT1 contains the forcing terms to the toroidal part     C
C of the vorticity equation and is stored in the same way. There     C
C are NH1 functions with l and m indices given by ML1 and MM1.       C
C                                                                    C
C The terms are added for nodes IR = 3 to IR = nr - 2.               C
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
C     NH3       : Number of poloidal velocity radial functions.      C
C     ML3       : Array length ( * ) - atleast length NH3            C
C                  Sph. harm. degree, l, for poloidal velocity.      C
C     MM3       : Array length ( * ) - atleast length NH3            C
C                  Sph. harm. order, m, for cos m phi dep. (temp.)   C
C                 -Sph. harm. order, m, for sin m phi dep. (temp.)   C
C                                                                    C
C     NH1       : Number of func.s in tor. vorticity eqn.            C
C     ML1       : Array length ( * ) - atleast length NH1            C
C                  Sph. harm. degree, l, for tor. vort. forcing t.   C
C     MM1       : Array length ( * ) - atleast length NH3            C
C                  Sph. harm. order, m, for cos m phi dep. (t. vort) C
C                 -Sph. harm. order, m, for sin m phi dep. (t. vort) C
C                                                                    C
C  Double Precision                                                  C
C  ----------------                                                  C
C                                                                    C
C     SV3       : Solution vector. Dim ( * ) ( temperature )         C
C                 Length must be atleast NR*NH3                      C
C     FT1       : Forcing term for toroidal part of vorticity eqn.   C
C                  Dim ( * ). Length must be atleast NR*NH1          C
C     CH        : Multiplier of term to be added.                    C
C                                                                    C
C____________________________________________________________________C
C
C*********************************************************************
      SUBROUTINE NSVBTA( NR, NH3, ML3, MM3, NH1, ML1, MM1, SV3, FT1,
     1                   CH )
      IMPLICIT NONE
C____________________________________________________________________C
C Variable declarations - Parameters ................................C
      INTEGER NR, NH3, ML3(NH3), MM3(NH3), NH1, ML1(NH1), MM1(NH1)
      DOUBLE PRECISION SV3( * ), FT1( * ), CH
C____________________________________________________________________C
C Variable declarations - Working variables .........................C
C
      INTEGER IR, L3, MMM3, IH3, IH1, IND3, IND1, IBEG3, IBEG1
      DOUBLE PRECISION ZERO
      PARAMETER ( ZERO = 0.0d0 )
C____________________________________________________________________C
C START OF PROGRAM **************************************************C
C____________________________________________________________________C
C
C     . early exit?
C
      IF ( CH.EQ.ZERO ) RETURN
C     .
C     . Loop around temperature harmonics
C     .
      DO IH3 = 1, NH3
        L3   = ML3( IH3 )
        MMM3 = MM3( IH3 )
C
C so let's look for the corresponding tor. vort. harmonic in FT1
C
        DO IH1 = 1, NH1
          IF (   ML1( IH1 ).EQ.L3 .AND. MM1( IH1 ).EQ.MMM3  ) THEN
C
C o.k. we've found the corresponding harmonic
C
            IBEG3 = (IH3 - 1)*NR
            IBEG1 = (IH1 - 1)*NR
            DO IR = 3, NR - 2
C
C Find locations of source and destination vectors ...
C
              IND3 = IBEG3 + IR
              IND1 = IBEG1 + IR
C
              FT1( IND1 ) = FT1( IND1 ) + CH*SV3( IND3 )
C
            ENDDO
C           . End of loop ir = 3, nr-2
C
C Jump out of ih1 loop if we have found our harmonic
C
            GOTO 46
C
          ENDIF
C         . check ml1(ih1).eq.l3 and mm1(ih1).eq.mm3
        ENDDO
C       . End of loop ih1 = 1, nh1
C
 46   CONTINUE
      ENDDO
C     . End of loop ih3 = 1, nh3
C
      RETURN
      END
C*********************************************************************


