C*********************************************************************
C subroutine Poloidal Velocity Vector ComPLete ***********************
C            -        -        -      -  --    ***********************
C Steve Gibbons Thu Nov  9 10:47:05 WET 2000                         C
C____________________________________________________________________C
C                                                                    C
C Because of the additional boundary conditions on the poloidal      C
C velocity radial function, the values at grid nodes 2 and NR-1 are  C
C not solved for in the solution vector.                             C
C PVVCPL completes these two values in the solution vector.          C
C                                                                    C
C The solution vector must consist only of poloidal velocity         C
C harmonics, of which there are NH1, and the values of these         C
C radial functions, at NR grid nodes, must be given by the index     C
C                                                                    C
C   ind = ( ih - 1 )*nr + ir                                         C
C                                                                    C
C The arrays PVLC and PVRC have been prepared by a call to PVCCF.    C
C                                                                    C
C____________________________________________________________________C
C                                                                    C
C Input variables :-                                                 C
C ===============                                                    C
C                                                                    C
C  Integer                                                           C
C  -------                                                           C
C                                                                    C
C     NR        : Number of grid nodes.                              C
C     NH1       : Number of poloidal harmonic functions.             C
C     NBN       : Maximum number of nodes on either side for         C
C                  central differencing.                             C
C                                                                    C
C  Double Precision                                                  C
C  ----------------                                                  C
C                                                                    C
C     SV1       : Solution vector for poloidal velocity.             C
C                  Dim (NR*NH1)                                      C
C     PVLC      : Dim (NBN). See PVCCF.                              C
C     PVRC      : Dim (NBN). See PVCCF.                              C
C____________________________________________________________________C
C
C*********************************************************************
      SUBROUTINE PVVCPL( NR, NH1, NBN, SV1, PVLC, PVRC )
      IMPLICIT NONE
C____________________________________________________________________C
C Variable declarations - Parameters ................................C
      INTEGER          NR, NH1, NBN
      DOUBLE PRECISION SV1( NR*NH1 ), PVLC( NBN ), PVRC( NBN )
C____________________________________________________________________C
C Variable declarations - Working variables .........................C
      INTEGER          I, IND, IH, IR, IND2
      DOUBLE PRECISION TEMP
C____________________________________________________________________C
C START OF PROGRAM **************************************************C
C____________________________________________________________________C
C
      DO IH = 1, NH1
C
        IR   = 2
        IND2 = ( IH - 1 )*NR + IR
        TEMP = 0.0d0
        DO I = 1, NBN
          IND = IND2 + I
          TEMP = TEMP + PVLC( I )*SV1( IND )
        ENDDO
        SV1( IND2 ) = TEMP
C
        IR   = NR - 1
        IND2 = ( IH - 1 )*NR + IR
        TEMP = 0.0d0
        DO I = 1, NBN
          IND = IND2 - I
          TEMP = TEMP + PVRC( I )*SV1( IND )
        ENDDO
        SV1( IND2 ) = TEMP
C
      ENDDO
C
      RETURN
      END
C*********************************************************************
