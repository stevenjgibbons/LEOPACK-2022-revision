C*********************************************************************
C integer function CHebyshev INDex Function **************************
C                  --        ---   -        **************************
C Steve Gibbons Fri Apr 28 09:00:05 BST 2000                         C
C____________________________________________________________________C
C                                                                    C
C Our solutions are normally held at NR grid points.                 C
C If these radial grid points are the zeros of the Chebshev          C
C polynomial T_{N} where N = NR - 2 (with the exception of grid      C
C nodes 1 and NR which are the inner and outer boundaries) then      C
C we can transform the radial functions f(x) into a Chebyshev basis  C
C                                                                    C
C  f( x ) = \sum_{k=0}^{N-1} c_k T_k( y )                            C
C                                                                    C
C where the transformed variable y is given by the subroutine        C
C CPCOVR.                                                            C
C                                                                    C
C Returns the location in the solution vector of the value of the    C
C coefficient c_k (k=IK) for the IH^{th} harmonic.                   C
C                                                                    C
C If IK = N or N+1, then we cannot be referring to a coefficient     C
C in the above expansion. However, we wish to store the values of    C
C the inner radius and outer radius as a reference.                  C
C                                                                    C
C Hence, CHINDF( N, IH, INARR ) is the location of RI and            C
C        CHINDF( N+1, IH, INARR ) is the location of RO.             C
C____________________________________________________________________C
C Input variables :-                                                 C
C ===============                                                    C
C  Integer                                                           C
C  -------                                                           C
C     IK        : Chebyshev coefficient index, k, for c_k.           C
C     IH        : Number of harmonic.                                C
C     INARR     : Array of integers. Dimension ( * )                 C
C                 At present the elements of INARR are as follows.   C
C                                                                    C
C                 INARR( 1 ) = IFORMF. Format flag.                  C
C                  This flag decides how the elements are arranged   C
C                  in the solution vector.                           C
C                                                                    C
C                  The current options are:-                         C
C                                                                    C
C                   IFORMF = 1/3. CHINDF = IK*NH + IH                C
C                   IFORMF = 2/4. CHINDF = ( IH - 1 )*NR + IK + 1    C
C                                                                    C
C                   (NR = N + 2, where nodes are zeros of T_N)       C
C                                                                    C
C  Here, NH is the TOTAL number of harmonics in the solution         C
C  vector - or atleast the part of which is visible to that part     C
C  of the program. NR is the number of radial grid nodes             C
C  corresponding to that harmonic (which for cases IFORMF = 1 and    C
C  IFORMF = 2 is identical for all harmonics - more complicated      C
C  options (for example magnetic fields which have to be resolved    C
C  beyond the region of fluid flow) may be added later and this      C
C  routine should be flexible to all possibilities with extra        C
C  constraints being added in other elements of INARR.               C
C                                                                    C
C  IFORMF = 1/3 is the option likely to be used for the purpose of   C
C  solution as it allows the banded formation of a matrix.           C
C                                                                    C
C  IFORMF = 2/4 is the option likely to be used for the purpose of   C
C  displaying solutions as it stores adjacent nodes for each         C
C  harmonic together.                                                C
C                                                                    C
C                 INARR( 2 ) = NR. See above.                        C
C                 INARR( 3 ) = NH. See above.                        C
C                                                                    C
C____________________________________________________________________C
C
C*********************************************************************
      FUNCTION CHINDF ( IK, IH, INARR )
      IMPLICIT NONE
C____________________________________________________________________C
C Variable declarations - Parameters ................................C
      INTEGER CHINDF, IK, IH, INARR( * )
C____________________________________________________________________C
C Variable declarations - Working variables .........................C
      INTEGER IFORMF, NR, NH, N
C____________________________________________________________________C
C START OF PROGRAM **************************************************C
C____________________________________________________________________C
C
      IFORMF = INARR( 1 )
      NR     = INARR( 2 )
      NH     = INARR( 3 )
      N      = NR - 2
C
      IF ( IK.LT.0 .OR. IK.GE.(N+2) ) THEN
        PRINT *,' Function CHINDF. IK = ', IK
        PRINT *,' N = ', N,'. Program aborted.'
        STOP
      ENDIF
C
      IF ( IH.LT.1 .OR. IH.GT.NH ) THEN
        PRINT *,' Function CHINDF. IH = ', IH
        PRINT *,' NH = ', NH,'. Program aborted.'
        STOP
      ENDIF
C
      IF ( IFORMF.EQ.1 .OR. IFORMF.EQ.3 ) THEN
        CHINDF = IK*NH + IH
        RETURN
      ENDIF
C
      IF ( IFORMF.EQ.2 .OR. IFORMF.EQ.4 ) THEN
        CHINDF = ( IH - 1 )*NR + IK + 1
        RETURN
      ENDIF
C
      PRINT *,' Function CHINDF. IFORMF = ', IFORMF
      PRINT *,' Not current option. Program aborted.'
      STOP
      END
C*********************************************************************
