C*********************************************************************
C subroutine Uniform Grid Matrix Interaction Contribution Add ********
C            -       -    -      -           -            -   ********
C Steve Gibbons Thu Sep 16 15:53:52 BST 1999                         C
C____________________________________________________________________C
C                                                                    C
C  Adds to the double precision matrix, A, the terms which a single  C
C harmonic IHC contributes to a single harmonic IHR in a single term C
C of a coupled o.d.e. when the grid in radius is uniform.            C
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
C The subroutine SUB1 *MUST* have the calling sequence               C
C                                                                    C
C SUB1( CVEC, IRAD, RAD, COEFM, NCFM, NLN, NRN, INARR, DPARR, VEC0,  C
C       IHD )                                                        C
C                                                                    C
C where  IRAD is the number of the radial grid node of the           C
C        harmonic which is being contributed to.                     C
C        RAD is the double precision value of the radius,            C
C        COEFM is the coefficient matrix of leading order NCFM,      C
C        NLN is the number of nodes to the left with which to        C
C          calculate derivatives.                                    C
C        NLR is the number of nodes to the right with which to       C
C          calculate derivatives.                                    C
C        INARR is an integer array as used by                        C
C          RADVLF and INDFUN although other elements of the array    C
C          may be used arbitrarily by the routine SUB1.              C
C        DPARR is a double precision array as used by                C
C          RADVLF and INDFUN although other elements of the array    C
C          may be used arbitrarily by the routine SUB1.              C
C        VEC0 is a solution vector whose elements may be required.   C
C        IHD is the highest derivative which will be refered to by   C
C        UGMICA. So all CVEC( i ) must be zero with 1.le.i.le.ihd+1  C
C        if it is not assigned another value.                        C
C                                                                    C
C        The parameters IRAD, COEFM, NCFM, NLN, NLR and VEC0 are     C
C        all included into the calling sequence incase this is       C
C        a non-linear term.                                          C
C                                                                    C
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
C     NCFM      : Dimension of square finite diff. coefficient       C
C                  matrix.                                           C
C                                                                    C
C     IPCM      : Work array for routine FDCINV                      C
C                 Dimension ( NCFM )                                 C
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
C                   IFORMF = 1,3. INDFUN = ( IR - 1 )*NH + IH        C
C                   IFORMF = 2,4. INDFUN = ( IH - 1 )*NR + IR        C
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
C     ILNC      : Left-most node which may be used to form deriv.s   C
C     IRNC      : Right-most node which may be used to form deriv.s  C
C                                                                    C
C  ILNC must never be greater than ILNR, and similarly,              C
C  IRNC must never be less than IRNR.                                C
C                                                                    C
C  Short word of explanation about ILNR, IRNR, ILNC, IRNC ...        C
C  If ILNR = ILNC = 1 and IRNR = IRNC = NR, then the contribution    C
C  to every node of IHR is calculated using every node of harmonic   C
C  IHC. However, normally one or two boundary conditions need to be  C
C  enforced at each boundary and so it makes sense not to spend time C
C  calculating matrix terms which will consequently be deleted.      C
C  Say 2 boundary conditions need to be enforced at IRAD = 1 and     C
C  1 at IRAD = NR, then set ILNR = 3 and                             C
C  IRNR = NR - 1. (This is also a way to get around a possible       C
C  division by zero error if grid point 1 is at the origin).         C
C  I can't think of any situations where you wouldn't want to        C
C  include either boundary point when taking a derivative ... unless C
C  the column is to be removed from the solution vector for some     C
C  reason. For generality, these are flexible terms.                 C
C                                                                    C
C  Subroutines                                                       C
C  -----------                                                       C
C     SUB1      : Determines what multiplies each derivative in      C
C                  the matrix. Must have calling sequence ...        C
C                                                                    C
C     CALL SUB1 ( CVEC, IRAD, RAD, COEFM, NCFM, NLN, NRN,            C
C                 INARR, DPARR, VEC0, IHD )                          C
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
C     WORK      : Work array for routine FDCINV                      C
C                 Dimension ( NCFM )                                 C
C                                                                    C
C     COEFM     : Dimension ( NCFM, NCFM ). Work array.              C
C                                                                    C
C     DPARR     : Array of dimension ( * ).                          C
C      DPARR( 1 ) = RI                                               C
C      DPARR( 2 ) = RO                                               C
C                                                                    C
C     VEC0      : Solution vector of dimension ( * ).                C
C                Only referred to if done so by subroutine SUB1.     C
C                                                                    C
C                                                                    C
C                                                                    C
C                                                                    C
C____________________________________________________________________C
C
C*********************************************************************
      SUBROUTINE UGMICA( N1, N2, KL, KU, KLE, IMF, IHC, IHR, NCFM,
     1                   IPCM, INARR, IHD, NBN, ILNR, IRNR, ILNC,
     2                   IRNC, SUB1, AMAT, FAC, WORK, COEFM, DPARR,
     3                   VEC0 )
      IMPLICIT NONE
C____________________________________________________________________C
C Variable declarations - Parameters ................................C
      INTEGER N1, N2, KL, KU, KLE, IMF, IHC, IHR, NCFM,
     1        IPCM( NCFM ), INARR( * ), IHD, NBN, ILNR, IRNR, ILNC,
     2        IRNC
      EXTERNAL SUB1
      DOUBLE PRECISION AMAT( N1, N2 ), FAC, WORK( NCFM ), 
     1                 COEFM( NCFM, NCFM ), DPARR( * ), VEC0( * )
C____________________________________________________________________C
C Variable declarations - Working variables .........................C
C
      INTEGER IRAD, INODE, NLN, NRN, INDFUN, INDR, INDC,
     1        IROW, ICOL, NR, NH, IFORMF, NLCS, NRCS, I,
     2        NDER, ICOROW, ICOCOL
      DOUBLE PRECISION RAD, H, LOW, XARR( 1 )
      PARAMETER ( LOW = 1.0d-9 )
      LOGICAL ONRC
C onrc is the 'no repeat calculation flag'.
C if onrc is .TRUE. then we do not calculate the
C finite difference coefficients again as they are
C identical to the previous ones.
C xarr is not actually referenced as this is a fixed grid routine
C____________________________________________________________________C
C START OF PROGRAM **************************************************C
C____________________________________________________________________C
C
      XARR( 1 ) = 0.0d0
C
C First we check for an easy exit - 
C if FAC = 0.0d0 there is no point in doing this  ...
C
      IFORMF = INARR( 1 )
C     .
C     . Check uniform format
C     .
      IF ( IFORMF.EQ.3 .OR. IFORMF.EQ.4 ) THEN
        PRINT *,' Subroutine UGMICA.'
        PRINT *,' IFORMF = ', IFORMF
        PRINT *,' This is a uniform grid routine.'
        PRINT *,' Program aborted.'
        STOP
      ENDIF
C     .
      IF ( ABS( FAC ).LT.LOW ) RETURN
C     .
      NR     = INARR( 2 )
      NH     = INARR( 3 )
C     .
C     . Check N2 for case of IFORMF = 1,2,3 and 4
C     .
      IF ( IFORMF.EQ.1 .OR. IFORMF.EQ.2 .OR. IFORMF.EQ.3 .OR.
     1                     IFORMF.EQ.4            ) THEN
        IF ( N2.NE.NH*NR ) THEN
          PRINT *,' Subroutine UGMICA. N2 = ', N2
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
         PRINT *,' Subroutine UGMICA. IMF = ',IMF
         PRINT *,' N1 = ', N1,' N2 = ',N2
         PRINT *,' KL = ', KL,' KU = ',KU
         PRINT *,' KLE = ', KLE
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
      NLCS = ILNR - ILNC
      NRCS = IRNC - IRNR
C
C Just check that we have sufficiently many grid nodes
C to be able to enforce all equations.
C
      IF ( NLCS.LT.0 ) THEN
        PRINT *,' Subroutine UGMICA.'
        PRINT *,' ILNR = ', ILNR
        PRINT *,' ILNC = ', ILNC
        PRINT *,' Program aborted.'
        STOP
      ENDIF
C
      IF ( NRCS.LT.0 ) THEN
        PRINT *,' Subroutine UGMICA.'
        PRINT *,' IRNR = ', IRNR
        PRINT *,' IRNC = ', IRNC
        PRINT *,' Program aborted.'
        STOP
      ENDIF
C
      I = MIN( NLCS, NRCS) + NBN
      IF ( IHD.GT.I ) THEN
        PRINT *,' Subroutine UGMICA.'
        PRINT *,' You want a derivative of order ', IHD
        PRINT *,' However, there is a node at which you '
        PRINT *,' only have ',I,' points for a derivative.'
        PRINT *,' Program aborted.'
        STOP
      ENDIF
C     .
C     . Ok - we now loop around the radial nodes for the
C     . rows ... first set our logical flag to .FALSE.
C     . i.e. calculate differencing coefficients each
C     . time unless told otherwise
C     .
      ONRC = .FALSE.
      DO IRAD = ILNR, IRNR
C       .
C       . Calculate the corresponding index and radius
C       .
        INDR = INDFUN( IRAD, IHR, INARR )
        CALL RADVLF( RAD, IRAD, INARR, DPARR, XARR, H )
C       .
C       . Calculate NLN and NRN
C       .
        NLN = MIN( NBN, IRAD - ILNC )
        NRN = MIN( NBN, IRNC - IRAD )
C       .
C       . Next line ensures that coefficients are
C       . calculated if it is not a central difference
C       . or if this is the first call to the central
C       . nodes ...
C       .
        IF ( ONRC .AND. NLN.EQ.NBN .AND. NRN.EQ.NBN ) GOTO 50
C       .
C       . Calculate the coefficients
C       .
        CALL FDCINV( H, NLN, NRN, NBN, NCFM, COEFM, IPCM, WORK )
C       .
        IF ( NLN.EQ.NBN .AND. NRN.EQ.NBN ) ONRC = .TRUE.
 50     CONTINUE
C       .
C       . ok - the finite difference coefficients
C       . are calculated so we can now call SUB1 and
C       . see what we have to multiply each of the deriv.s by.
C       . Now the array WORK will not be required until
C       . the next call of FDCINV and so we can store the
C       . coefficients in it that SUB1 returns ...
C       .
        CALL SUB1( WORK, IRAD, RAD, COEFM, NCFM, NLN, NRN,
     1             INARR, DPARR, VEC0, IHD )
C       .
C       . So now loop around the nodes from (IRAD - NLN)
C       . to (IRAD + NRN)
C       .
        DO INODE = IRAD - NLN, IRAD + NRN
C          .
C          . This is the column of the finite
C          . difference coefficients matrix
C          . which will contain the correct values.
C          .
           ICOCOL = INODE - IRAD + NLN + 1
C          .
C          . find the corresponding matrix column
C          .
           INDC = INDFUN( INODE, IHC, INARR )
C          .
C          . Find the actual location in the AMAT array
C          .
           CALL MATIND (INDR,INDC,IMF,KL,KU,KLE,N1,N2,IROW,ICOL)
C          .
C          .
C          .
           DO NDER = 0, IHD
             ICOROW = NDER + 1
             AMAT( IROW, ICOL ) = AMAT( IROW, ICOL ) + 
     1           FAC*COEFM( ICOROW, ICOCOL )*WORK( ICOROW )
           ENDDO
C          .
        ENDDO
C       .
      ENDDO
C     .
      RETURN
      END
C*********************************************************************
