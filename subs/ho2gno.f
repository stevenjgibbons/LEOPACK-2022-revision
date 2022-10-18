C*********************************************************************
C subroutine Harmonic Order 2 Grid Node Order ************************
C            -        -     - -    -    -     ************************
C Steve Gibbons Fri Feb 16 12:40:06 WET 2001                         C
C____________________________________________________________________C
C                                                                    C
C Reads a vector VECIN whose elements are indexed with               C
C  IND = (IH-1)*NR + IR                                              C
C and puts the values into a vector VECOUT, whose elements are       C
C indexed with                                                       C
C  IND = (IR-1)*NH + IH                                              C
C                                                                    C
C____________________________________________________________________C
C                                                                    C
C  Integer                                                           C
C  -------                                                           C
C                                                                    C
C     NR        : Number of radial grid nodes.                       C
C     NH        : Number of spherical harmonics radial functions.    C
C                                                                    C
C  Double Precision                                                  C
C  ----------------                                                  C
C                                                                    C
C     VECIN     : Dim (NR*NH). Vector stored in "harmonic order"     C
C                 i.e. harmonic 1 at all grid nodes comes before     C
C                 harmonic 2 at all nodes  etc.                      C
C                    [ IND = (IH-1)*NR + IR ]                        C
C                                                                    C
C     VECOUT    : Dim (NR*NH). Vector stored in "grid node order"    C
C                 i.e. every harmonic at node 1 comes before every   C
C                 harmonic at node 2, etc. [ IND = (IR-1)*NH + IH ]  C
C                                                                    C
C____________________________________________________________________C
C
C*********************************************************************
      SUBROUTINE HO2GNO( NR, NH, VECIN, VECOUT )
      IMPLICIT NONE
C____________________________________________________________________C
C Variable declarations - Parameters ................................C
      INTEGER          NR, NH
      DOUBLE PRECISION VECIN( * ), VECOUT( * )
C____________________________________________________________________C
C Variable declarations - Working variables .........................C
C
      INTEGER IH, IR, INDOUT, INDIN
C____________________________________________________________________C
C START OF PROGRAM **************************************************C
C____________________________________________________________________C
C
      INDIN = 0
      DO IH = 1, NH
        DO IR = 1, NR
C         .
C         . The following incrementation of indin
C         . should be equivalent to saying
C         . INDIN  = (IH-1)*NR + IR
C         .
          INDIN  = INDIN + 1
          INDOUT = (IR-1)*NH + IH
          VECOUT( INDOUT ) = VECIN( INDIN )
        ENDDO
      ENDDO
C
      RETURN
      END
C*********************************************************************
