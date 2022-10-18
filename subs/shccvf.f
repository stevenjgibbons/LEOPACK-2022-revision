C*********************************************************************
C subroutine Single Harmonic Chebyshev Coefficient Vector Form *******
C            -      -        -         -           -      -    *******
C Steve Gibbons Fri Apr 28 15:40:38 BST 2000                         C
C____________________________________________________________________C
C                                                                    C
C SV is a standard solution vector arranged according to INDFUN and  C
C INARR with XARR containing the NR grid points, such that           C
C XARR( i + 1 ) is the i^{th} zero of the Chebyshev polynomial T_N,  C
C where N = NR - 2.                                                  C
C                                                                    C
C THIS ***MUST*** be the case! Otherwise it is meaningless to run    C
C this subroutine!!! It is ***NOT*** checked for.                    C
C                                                                    C
C Now if IH is the number of a harmonic between 1 and NH, then       C
C SHCCVF evaluates the coefficients c_k, where the IH^{th} radial    C
C function f(r) is given by                                          C
C                                                                    C
C f(r) = \sum_{k=0}^{N-1} c_k T_k( x ) where x is in the interval    C
C                                      [-1,1] (see CPCOVR).          C
C                                                                    C
C The coefficient c_k is put in the array CCV in the location given  C
C by the integer function CHINDF.                                    C
C                                                                    C
C The procedure is directly based on the subroutine COCHGA by        C
C Daniele Funaro (see my subroutine COCHGV).                         C
C                                                                    C
C____________________________________________________________________C
C                                                                    C
C Input variables :-                                                 C
C ===============                                                    C
C  Integer                                                           C
C  -------                                                           C
C                                                                    C
C     INARR     : Integer array dimension ( * ).                     C
C                  Elements may be arbitrary except for              C
C                  INARR( 1 ) = IFORMF - flag for vector format.     C
C                                        See INDFUN                  C
C                  INARR( 2 ) = NR. Number of radial grid nodes.     C
C                  INARR( 3 ) = NH = total number of radial func.s   C
C     IH        : Number of radial function (harmonic).              C
C                                                                    C
C  Double Precision                                                  C
C  ----------------                                                  C
C                                                                    C
C     SV        : Solution vector. Dim ( * ) but length atleast      C
C                  NR*NH.                                            C
C                                                                    C
C     CCV       : Chebyshev coefficient vector. Dim ( * ) but        C
C                 length atleast NR*NH. For each radial function IH, C
C                 the element CCV( NR-2, IH, INARR ) gives RI and    C
C                 CCV( NR-1, IH, INARR ) gives RO.                   C
C                                                                    C
C     XARR      : Dim (*). Radial values. No check is made on them.  C
C____________________________________________________________________C
C
C*********************************************************************
      SUBROUTINE SHCCVF( INARR, IH, SV, CCV, XARR )
      IMPLICIT NONE
C____________________________________________________________________C
C Variable declarations - Parameters ................................C
      INTEGER INARR( * ), IH
      DOUBLE PRECISION SV( * ), CCV( * ), XARR( * )
C____________________________________________________________________C
C Variable declarations - Working variables .........................C
      INTEGER IR, K, J, N, INDSV, INDCCV, INDFUN, CHINDF
      DOUBLE PRECISION PH, DN, SU, SK, DJ, DK, RI, RO
C____________________________________________________________________C
C START OF PROGRAM **************************************************C
C____________________________________________________________________C
C
      N  = INARR( 2 ) - 2
      RI = XARR( 1 )
      RO = XARR( N + 2 )
      PH = 1.57079632679489661923D0
      DN = DFLOAT(N)
      SU = 0.0D0
      DO 1 J = 1, N
        IR    = J + 1
        INDSV = INDFUN( IR, IH, INARR )
        SU    = SU+SV( INDSV )
 1    CONTINUE
      K             = 0
      INDCCV        = CHINDF( K, IH, INARR )
      CCV( INDCCV ) = SU/DN
      IF (N .EQ. 1) RETURN
          SK = -2.0D0
      DO 2 K=1,N-1
          DK = DFLOAT(K)
          SU = 0.0D0
      DO 3 J = 1, N
          IR    = J + 1
          INDSV = INDFUN( IR, IH, INARR )
          DJ    = 2.0D0*DFLOAT(J)-1.0D0
          SU    = SU+SV( INDSV )*DCOS(DK*DJ*PH/DN)
 3    CONTINUE
      INDCCV        = CHINDF( K, IH, INARR )
      CCV( INDCCV ) = SK*SU/DN
      SK            = -SK
 2    CONTINUE
C
C Now enter the inner and outer radii in the spare locations
C
      K             = N
      INDCCV        = CHINDF( K, IH, INARR )
      CCV( INDCCV ) = RI
C
      K             = N+1
      INDCCV        = CHINDF( K, IH, INARR )
      CCV( INDCCV ) = RO
C
      RETURN
      END
C*********************************************************************
