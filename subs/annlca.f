C*********************************************************************
C subroutine Arbitrary Node Non-Linear Contribution Add **************
C            -         -    -   -      -            -   **************
C Steve Gibbons Sun Nov 21 14:54:41 GMT 1999                         C
C____________________________________________________________________C
C                                                                    C
C  Adds to the double precision matrix, A, the terms which a single  C
C harmonic IHC contributes to a single harmonic IHR as a result of   C
C interaction with a harmonic IH0.                                   C
C                                                                    C
C  The radial functions to harmonics IHC and IHR, which are          C
C represented in the solution vector are represented at NR grid      C
C nodes whose spacing is dictated by the array XARR.                 C
C                                                                    C
C  The harmonic IH0 need not necessarily appear in this solution     C
C vector and it is stored in the vector VEC0 which is indexed by     C
C the array INARR0 (NR0 not necessarily equal to NR) and the         C
C spacings are given by XARR0.                                       C
C                                                                    C
C Let the radial function of the IHR harmonic be denoted f_{IHR)(r)  C
C Let the radial function of the IHC harmonic be denoted f_{IHC)(r)  C
C                                                                    C
C If the contribution for a given radial node, $r_j$, is of the form C
C                                                                    C
C f_{IHR)(r_j) = \sum_{nd=0,IHD} f_{IHC}^{nd}(r_j) c_{IHC}^{nd}(r_j) C
C                                                                    C
C  (where f^{nd} denotes the (nd)^{th} derivative of f with respect  C
C  to r, evaluated at r_j)                                           C
C                                                                    C
C then the coefficient c_{IHC}^{nd}(r_j) must be supplied by the     C
C subroutine SUB1 (declared EXTERNAL in the calling (sub)program )   C
C which returns c_{IHC}^{nd}(r_j) in the array element CVEC(nd + 1). C
C                                                                    C
C Unlike in AMLICA, where c_{IHC}^{nd} depends only upon r and       C
C some constants (which are passed in through the arrays IPARS and   C
C DPARS), c_{IHC}^{nd} will be a function of VEC0.                   C
C                                                                    C
C IHD0 is the highest derivative which is required of VEC0.          C
C                                                                    C
C The derivative of VEC0 is taken by the routine SVRINT.             C
C                                                                    C
C NNDS is the number of nodes which may be used to calculate this    C
C derivative. NNDS must be atleast (IHD0+1) and for any accuracy     C
C to speak of should be atleast (IHD0+2).                            C
C                                                                    C
C This can allow for a far greater accuracy in taking the deriv.s    C
C of the functions of VEC0 than may be possible in the matrix where  C
C the equations are to be solved.                                    C
C                                                                    C
C The resolution may ofcourse be identical.                          C
C                                                                    C
C The subroutine SUB1 *MUST* have the calling sequence               C
C                                                                    C
C SUB1( CVEC, RAD, IPARS, DPARS, IHD, IH0, IHD0, NNDS,               C
C       XARR0, INARR0, VEC0, WORK1, WORK2, IWORK, WORKM )            C
C                                                                    C
C where  RAD is the double precision value of the radius,            C
C        IPARS is an integer array to provide SUB1 with parameters.  C
C          IPARS is not referred to by ANNLCA other than to pass     C
C          this information.                                         C
C        DPARS is like IPARS but contains double precision elements. C
C        IHD is the highest derivative which will be refered to by   C
C        ANNLCA. So all CVEC( i ) must be zero with 1.le.i.le.ihd+1  C
C        if it is not assigned another value.                        C
C                                                                    C
C        Other variables are described above.                        C
C                                                                    C
C____________________________________________________________________C
C Input variables :-                                                 C
C ===============                                                    C
C  Integer                                                           C
C  -------                                                           C
C     N1        : Leading dimension of the matrix.                   C
C     N2        : Second dimension of the matrix.                    C
C     KL        : Number of lower diagonals in banded matrix.        C
C     KU        : Number of upper diagonals in banded matrix.        C
C     KLE       : Number of extra lower diagonals in banded          C
C                  matrix. This is to make the matrix solvable by    C
C                   LAPACK routines.                                 C
C                                                                    C
C     IMF       : Matrix format flag.                                C
C                                                                    C
C          imf = 1; Matrix is in LAPACK banded format                C
C                   ie element a_{i,j} is stored in                  C
C                   A( kle + ku + 1 + i - j , j )                    C
C                                                                    C
C          imf = 2; Matrix is banded but with element a_{i,j}        C
C                   stored in A( kl + 1 + j - i , i ).               C
C                                                                    C
C          imf = 3; Matrix is square - ie a_{i,j} is stored          C
C                   in A( i, j ).                                    C
C                                                                    C
C     IHC       : Number of the input (column) harmonic              C
C     IHR       : Number of the output (row) harmonic                C
C                                                                    C
C     INARR     : Array of integers. Dimension ( * )                 C
C                 At present the elements of INARR are as follows.   C
C                                                                    C
C                 INARR( 1 ) = IFORMF. Format flag.                  C
C                  This flag decides how the elements are arranged   C
C                  in the solution vector.                           C
C                                                                    C
C                  The current options are:-                         C
C                                                                    C
C                   IFORMF = 1. INDFUN = ( IR - 1 )*NH + IH          C
C                   IFORMF = 2. INDFUN = ( IH - 1 )*NR + IR          C
C                                                                    C
C  where IR and IH are the current grid node and harmonic resp.      C
C  and NR and NH are the total numbers of nodes and harmonics        C
C  in the solution vector.                                           C
C                                                                    C
C                 INARR( 2 ) = NR. Number of radial grid nodes.      C
C                 INARR( 3 ) = NH. Number of harmonics in sol. vect. C
C                                                                    C
C     IHD       : Highest derivative involved.                       C
C     NBN       : Number of diagonal elements in radius.             C
C     ILNR      : Left-most node at which to start adding rows to    C
C                   matrix.                                          C
C     IRNR      : Right-most node at which to start adding rows to   C
C                   matrix.                                          C
C                                                                    C
C     NFDCM     : Leading dim of SVFDC. See SVFDCF.                  C
C                  (Must be atleast 2*NBN + 1 )                      C
C                                                                    C
C     NR        : Number of radial grid nodes.                       C
C     NDCS      : Maximum distinct finite difference schemes         C
C                  stored by the array SVFDC.                        C
C     IS        : Number of required finite difference scheme.       C
C                  (This points to the fourth dim. element in        C
C                   the coefficient array SVFDC.)                    C
C                                                                    C
C     NDRVS     : Highest derivative stored in SVFDC.                C
C                 (Must be atleast 1).                               C
C     NDRVM     : Limit on NDRVS. Array bound for SVFDC.             C
C                                                                    C
C     IPARS     : Information array for SUB1. Dimension ( * )        C
C                                                                    C
C     IH0       : Harmonic with which IHC interacts to produce IHR.  C
C     IHD0      : Highest deriv. required from VEC0.                 C
C     NNDS      : Number of nodes which may be used to calculate     C
C                  derivatives in array VEC0.                        C
C                  (MUST be atleast IHD0+1 - and should be no less   C
C                   than IHD0+2 ....                                 C
C                                                                    C
C     INARR0    : See INDFUN. Format of VEC0.                        C
C                                                                    C
C     IWORK     : Work array. Dim ( NNDS ).                          C
C                                                                    C
C  Subroutines                                                       C
C  -----------                                                       C
C     SUB1      : Determines what multiplies each derivative in      C
C                  the matrix. Must have calling sequence ...        C
C                                                                    C
C     CALL SUB1 ( CVEC, RAD, IPARS, DPARS, IHD, IH0, IHD0, NNDS,     C
C                 XARR0, INARR0, VEC0, WORK1, WORK2, IWORK, WORKM )  C
C                                                                    C
C  Double Precision                                                  C
C  ----------------                                                  C
C                                                                    C
C     AMAT      : Matrix. Dimensions ( N1, N2 )                      C
C              Will generally be banded due to the nature of the     C
C              numerical scheme. KL, KU and KLE parameterise this.   C
C                                                                    C
C     FAC       : Multiplication factor.                             C
C                                                                    C
C     XARR      : Array of dimension ( NR )                          C
C                 XARR( j ) = element x_j                            C
C                                                                    C
C     WORK      : Working array of dimension ( NDRVS + 1 )           C
C                                                                    C
C     DPARS     : Information array for SUB1. Dimension ( * )        C
C                                                                    C
C     SVFDC     : Finite difference coefficient matrix.              C
C                  Dimension ( NFDCM, NR, NDRVM+1, NDCS ).           C
C                   Array is generated by the routine svfdcf         C
C                 See documentation for SVFDCF for details.          C
C                                                                    C
C     VEC0      : Vector containing 'old function'                   C
C                                                                    C
C     XARR0     : Array of dimension ( NR0 )                         C
C                 XARR0( j ) = element x_j in storage of VEC0.       C
C                                                                    C
C     WORK1     : Work array dim ( NNDS )                            C
C     WORK2     : Work array dim ( NNDS )                            C
C     WORKM     : Work array dim ( NNDS, NNDS )                      C
C                                                                    C
C____________________________________________________________________C
C
C*********************************************************************
      SUBROUTINE ANNLCA( N1, N2, KL, KU, KLE, IMF, IHC, IHR, INARR,
     1      IHD, NBN, ILNR, IRNR, NFDCM, NR, NDCS, IS, NDRVS, NDRVM,
     2      IPARS, SUB1, AMAT, FAC, XARR, WORK, DPARS, SVFDC, IH0,
     3      IHD0, NNDS, INARR0, VEC0, XARR0, IWORK, WORK1, WORK2,
     4      WORKM )
      IMPLICIT NONE
C____________________________________________________________________C
C Variable declarations - Parameters ................................C
      INTEGER N1, N2, KL, KU, KLE, IMF, IHC, IHR, INARR( * ), IHD,
     1        NBN, ILNR, IRNR, NFDCM, NR, NDCS, IS, NDRVS, NDRVM,
     2        IPARS( * ), IH0, IHD0, NNDS, INARR0( * ), IWORK( NNDS )
      EXTERNAL SUB1
      DOUBLE PRECISION AMAT( N1, N2 ), FAC, XARR( NR ),
     1                 WORK( NDRVS + 1 ), DPARS( * ),
     2                 SVFDC( NFDCM, NR, NDRVM+1, NDCS ), VEC0( * ),
     3                 WORK1( NNDS ), WORK2( NNDS ),
     4                 WORKM( NNDS, NNDS ), XARR0( * )
C____________________________________________________________________C
C Variable declarations - Working variables .........................C
C
      INTEGER IRAD, INODE, NLN, NRN, INDFUN, INDR, INDC,
     1        IROW, ICOL, NH, IFORMF, NLCS, NRCS, I,
     2        NDER, ICOROW, ICOCOL, NRR
      DOUBLE PRECISION RAD, LOW
      PARAMETER ( LOW = 1.0d-9 )
C____________________________________________________________________C
C START OF PROGRAM **************************************************C
C____________________________________________________________________C
C
C First we check for an easy exit - 
C if FAC = 0.0d0 there is no point in doing this  ...
C
C     .
      IF ( ABS( FAC ).LT.LOW ) RETURN
C     .
      IFORMF = INARR( 1 )
      NRR    = INARR( 2 )
      NH     = INARR( 3 )
C     .
C     . Check validity of NRR. Because it is an array
C     . dimension, NR has to be passed in the parameter
C     . list and so we must ensure we are consistent
C     .
      IF ( NRR.NE.NR ) THEN
         PRINT *,' Subroutine ANNLCA.'
         PRINT *,' INARR( 2 ) = ', NRR
         PRINT *,' NR = ', NR
         PRINT *,' Program aborted.'
         STOP
      ENDIF
C     .
C     . Check that we really are dealing with
C     . a non-uniform grid.
C     .
      IF ( IFORMF.NE.3 .AND. IFORMF.NE.4 ) THEN
        PRINT *,' Subroutine ANNLCA. IFORMF = ',IFORMF
        PRINT *,' This is not the uniform grid code.'
        PRINT *,' Program aborted.'
        STOP
      ENDIF
C     .
C     . Check size of IHD
C     .
      IF ( IHD.GT.NDRVS ) THEN
        PRINT *,' Subroutine ANNLCA. '
        PRINT *,' IHD   = ', IHD
        PRINT *,' NDRVS = ', NDRVS
        STOP
      ENDIF
C     .
C     . Check N2 for case of IFORMF = 3 and 4
C     .
      IF ( IFORMF.EQ.3 .OR. IFORMF.EQ.4 ) THEN
        IF ( N2.NE.NH*NR ) THEN
          PRINT *,' Subroutine ANNLCA. N2 = ', N2
          PRINT *,' NH = ', NH,' NR = ', NR
          PRINT *,' Program aborted.'
          STOP
        ENDIF
      ENDIF
C     .
C     . Check N1 for IMF.EQ.1 and IMF.EQ.3
C     .
      IF ( ( IMF.EQ.1 .AND. N1.NE.(KLE+KL+KU+1) ) .OR.
     1     ( IMF.EQ.3 .AND. N1.NE.N2           )      ) THEN
        PRINT *,' Subroutine ANNLCA. IMF = ',IMF
        PRINT *,' N1 = ', N1,' N2 = ',N2
        PRINT *,' KL = ', KL,' KU = ',KU
        PRINT *,' KLE = ', KLE
        PRINT *,' Program aborted.'
        STOP
      ENDIF
C     .
      IF ( NNDS.LT.(IHD0+1) ) THEN
        PRINT *,' Subroutine ANNLCA.'
        PRINT *,' NNDS  = ', NNDS
        PRINT *,' IHD0  = ', IHD0
        PRINT *,' Program aborted.'
        STOP
      ENDIF
C     .
C     . Now we will check that we have sufficient
C     . band width to be able to calculate the
C     . IHD^{th} derivative of a function.
C     . IHD must be less than or equal to the
C     . minimum value of (NLN + NRN)
C     . 
C     . To do this, we calculate NLCS and NRCS
C     . which are respectively the minimum number
C     . of points you will ever find on the left
C     . and the right of your node, IRAD.
C     .
      NLCS = ILNR - 1
      NRCS = NR   - IRNR
C
C Just check that we have sufficiently many grid nodes
C to be able to enforce all equations.
C
      IF ( NLCS.LT.0 ) THEN
        PRINT *,' Subroutine ANNLCA.'
        PRINT *,' ILNR = ', ILNR
        PRINT *,' Program aborted.'
        STOP
      ENDIF
C
      IF ( NRCS.LT.0 ) THEN
        PRINT *,' Subroutine ANNLCA.'
        PRINT *,' IRNR = ', IRNR
        PRINT *,' Program aborted.'
        STOP
      ENDIF
C
      I = MIN( NLCS, NRCS) + NBN
      IF ( IHD.GT.I ) THEN
        PRINT *,' Subroutine ANNLCA.'
        PRINT *,' You want a derivative of order ', IHD
        PRINT *,' However, there is a node at which you '
        PRINT *,' only have ',I,' points for a derivative.'
        PRINT *,' Program aborted.'
        STOP
      ENDIF
C     .
C     . Ok - we now loop around the radial nodes for the
C     . rows ... 
C     .
      DO IRAD = ILNR, IRNR
C       .
C       . Calculate the corresponding index and radius
C       .
        INDR = INDFUN( IRAD, IHR, INARR )
        RAD = XARR( IRAD )
C       .
C       . Calculate NLN and NRN
C       .
        NLN = MIN( NBN, IRAD - 1 )
        NRN = MIN( NBN, NR - IRAD )
C       .
        CALL SUB1( WORK, RAD, IPARS, DPARS, IHD, IH0, IHD0, NNDS,
     1             XARR0, INARR0, VEC0, WORK1, WORK2, IWORK, WORKM )
C       .
C       . So now loop around the nodes from (IRAD - NLN)
C       . to (IRAD + NRN)
C       .
        DO INODE = IRAD - NLN, IRAD + NRN
C          .
C          . This is the row of the finite
C          . difference coefficients matrix
C          . which will contain the correct values.
C          .
           ICOCOL = INODE - IRAD + NBN + 1
C          .
C          . find the corresponding matrix column
C          .
           INDC = INDFUN( INODE, IHC, INARR )
C          .
C          . Find the actual location in the AMAT array
C          .
           CALL MATIND (INDR,INDC,IMF,KL,KU,KLE,N1,N2,IROW,ICOL)
C          .
           DO NDER = 0, IHD
             ICOROW = NDER + 1
             AMAT( IROW, ICOL ) = AMAT( IROW, ICOL ) + 
     1         FAC*WORK( ICOROW )*SVFDC( ICOCOL, IRAD, ICOROW, IS )
           ENDDO
C          .
        ENDDO
C       .
      ENDDO
C     .
      RETURN
      END
C*********************************************************************
