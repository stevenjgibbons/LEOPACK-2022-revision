C*********************************************************************
C subroutine Uniform Grid Matrix Boundary Condition Enforce **********
C            -       -    -      -        -         -       **********
C Steve Gibbons Fri Sep 17 08:58:45 BST 1999                         C
C____________________________________________________________________C
C                                                                    C
C  Replaces the equation for grid node IRR of harmonic IHR with      C
C  an equation for harm. IHC at grid node IRAD.                      C
C                                                                    C
C  Now UGMBCE does NOT zero the row of the matrix - it will assume   C
C  this has already been done. If this row is not zero then a call   C
C  must be made to BMRCOP with N2 = NDIM, IRC = 1,                   C
C  NRC = INDFUN( IRR, IHR, INARR ), ALPHA = BETA = 0.0d0.            C
C                                                                    C
C  If the condition is of the form                                   C
C                                                                    C
C    \sum_{nd = 0,IHD} f_{IHC}^{nd}( node IRAD ) c_{IHC}^{nd}(IRAD)  C
C    = BVALUE                                                        C
C                                                                    C
C  Then call UGMBCE with an EXTERNAL subroutine SUB1 which will      C
C  return c_{IHC}^{nd}(IRAD) in the d.p. array CVEC( nd + 1 )        C
C  and when solving the system make sure that                        C
C  RHS( IND ) = BVALUE ... where IND = INDFUN( IRR, IHR, INARR )     C
C                                                                    C
C  For example if the boundary condition to be enforced on           C
C  harmonic IHC at node NR is                                        C
C                                                                    C
C  dP/dr + (l-1)P/r = 0 (where P represents the radial function      C
C  to harmonic IHC) then SUB1 must return                            C
C                                                                    C
C   CVEC( 1 ) = (L-1)/RAD                                            C
C   CVEC( 2 ) = 1.0d0                                                C
C   CVEC( i ) = 0.0d0  for all i greater than 2                      C
C                                                                    C
C  and BVALUE is set to zero.                                        C
C  In almost all circumstances, IHC will be equal to IHR as all      C
C  boundary conditions should apply to harmonics independent of      C
C  other harmonics - but should the occasion arise ...               C
C                                                                    C
C The subroutine SUB1 *MUST* have the calling sequence               C
C                                                                    C
C SUB1( CVEC, IRAD, RAD, COEFM, NCFM, NLN, NRN, INARR, DPARR, VEC0,  C
C       IHD )                                                        C
C                                                                    C
C where  IRAD is the number of the grid nodes at which the           C
C          condition must be enforced.                               C
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
C        UGMBCE. So all CVEC( i ) must be zero with 1.le.i.le.ihd+1  C
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
C     IRR       : Grid node of harmonic IHR whose equation is to be  C
C                  replaced by the new constraint.                   C
C                                                                    C
C     IRAD      : Grid node at which the condition applies.          C
C                 This is not necessarily the same as IRR. For ex.,  C
C                 a poloidal velocity harmonic with stress free      C
C                 boundaries must satisfy p( r_o ) = p''( r_o ) = 0. C
C                 Therefore two boundary cond.s must be enforced and C
C                 so two of the equations for harmonic IHR at this   C
C                 boundary must be replaced with b.c. equations.     C
C                 The first can be set by calling UGMBCE with IRAD   C
C                 = NR ( where the b.c. must be enforced) and IRR    C
C                 = NR (which node's equation must be replaced ).    C
C                 The second b.c. must also be applied with IRAD =   C
C                 NR but IRR must be set equal to (say) (NR - 1).    C
C                                                                    C
C     Care should be taken to ensure IRAD is not too far from IRR    C
C     although the routine MATIND will trap any attempt to write     C
C     outside the bounds of the matrix.                              C
C                                                                    C
C     ILNC      : Left-most node which may be used to form deriv.s   C
C     IRNC      : Right-most node which may be used to form deriv.s  C
C                                                                    C
C  Short word of explanation about ILNC and ILNC ...                 C
C  We want to enforce a value at the node IRR of the harmonic        C
C  IHR and so this will require some other nodes around IRR,         C
C  depending upon the number of the higest derivative this condition C
C  requires.                                                         C
C                                                                    C
C  Being for boundary conditions, this routine will almost           C
C  inevitably wish to use the nodes corresponding to 1 and NR.       C
C  The only exception is for o.d.e.s where a condition is enforced   C
C  in the middle of an interval - still, there are few circumstances C
C  under which you would not wish the end nodes.                     C
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
C____________________________________________________________________C
C
C*********************************************************************
      SUBROUTINE UGMBCE( N1, N2, KL, KU, KLE, IMF, IHC, IHR, NCFM,
     1                   IPCM, INARR, IHD, NBN, IRR, IRAD, ILNC, IRNC,
     2                   SUB1, AMAT, FAC, WORK, COEFM, DPARR, VEC0 )
      IMPLICIT NONE
C____________________________________________________________________C
C Variable declarations - Parameters ................................C
      INTEGER N1, N2, KL, KU, KLE, IMF, IHC, IHR, NCFM,
     1        IPCM( NCFM ), INARR( * ), IHD, NBN, IRR, IRAD, ILNC,
     2        IRNC
      EXTERNAL SUB1
      DOUBLE PRECISION AMAT( N1, N2 ), FAC, WORK( NCFM ), 
     1                 COEFM( NCFM, NCFM ), DPARR( * ), VEC0( * )
C____________________________________________________________________C
C Variable declarations - Working variables .........................C
C
      INTEGER INODE, NLN, NRN, INDFUN, INDR, INDC,
     1        IROW, ICOL, NR, NH, IFORMF, I,
     2        NDER, ICOROW, ICOCOL
      DOUBLE PRECISION RAD, H, LOW, XARR( 1 )
      PARAMETER ( LOW = 1.0d-9 )
C xarr is not required as this is a uniform grid routine
C____________________________________________________________________C
C START OF PROGRAM **************************************************C
C____________________________________________________________________C
C
      XARR( 1 ) = 0.0d0
C
C First we check for an easy exit - 
C if FAC = 0.0d0 there is no point in doing this  ...
C     .
      IF ( ABS( FAC ).LT.LOW ) RETURN
C     .
      IFORMF = INARR( 1 )
      NR     = INARR( 2 )
      NH     = INARR( 3 )
C     .
C     . Check that we really are dealing with
C     . a uniform grid.
C     .
      IF ( IFORMF.NE.1 .AND. IFORMF.NE.2 ) THEN
        PRINT *,' Subroutine UGMBCE. IFORMF = ',IFORMF
        PRINT *,' Grid must be uniform for this routine.'
        PRINT *,' Program aborted.'
        STOP
      ENDIF
C     .
C     . Check N2 for case of IFORMF = 1,2,3 and 4
C     .
      IF ( IFORMF.EQ.1 .OR. IFORMF.EQ.2 .OR. IFORMF.EQ.3 .OR.
     1                     IFORMF.EQ.4            ) THEN
        IF ( N2.NE.NH*NR ) THEN
          PRINT *,' Subroutine UGMBCE. N2 = ', N2
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
         PRINT *,' Subroutine UGMBCE. IMF = ',IMF
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
C     . value of (NLN + NRN)
C     . 
C     . Calculate the corresponding index and radius
C     .
      INDR = INDFUN( IRR, IHR, INARR )
      CALL RADVLF( RAD, IRAD, INARR, DPARR, XARR, H )
C     .
C     . Calculate NLN and NRN
C     . As we are applying a boundary condition, one of these
C     . is likely to be zero ...
C     .
      NLN = MIN( NBN, IRAD - ILNC )
      NRN = MIN( NBN, IRNC - IRAD )
      I = ( NLN + NRN )
      IF ( IHD.GT.I ) THEN
        PRINT *,' Subroutine UGMBCE.'
        PRINT *,' You want a derivative of order ', IHD
        PRINT *,' However, there are only '
        PRINT *, I,' points for a derivative.'
        PRINT *,' Program aborted.'
        STOP
      ENDIF
C
      IF ( NLN.LT.0 ) THEN
        PRINT *,' Subroutine UGMBCE.'
        PRINT *,' IRAD = ', IRAD
        PRINT *,' NLN  = ', NLN
        PRINT *,' Program aborted.'
        STOP
      ENDIF
C
      IF ( NRN.LT.0 ) THEN
        PRINT *,' Subroutine UGMBCE.'
        PRINT *,' IRAD = ', IRAD
        PRINT *,' NRN  = ', NRN
        PRINT *,' Program aborted.'
        STOP
      ENDIF
C
C     .
C     . Calculate the coefficients
C     .
      CALL FDCINV( H, NLN, NRN, NBN, NCFM, COEFM, IPCM, WORK )
C     .
C     . ok - the finite difference coefficients
C     . are calculated so we can now call SUB1 and
C     . see what we have to multiply each of the deriv.s by.
C     . Now the array WORK will not be required until
C     . the next call of FDCINV and so we can store the
C     . coefficients in it that SUB1 returns ...
C     .
      CALL SUB1( WORK, IRAD, RAD, COEFM, NCFM, NLN, NRN,
     1           INARR, DPARR, VEC0, IHD )
C     .
C     . So now loop around the nodes from (IRAD - NLN)
C     . to (IRAD + NRN)
C     .
      DO INODE = IRAD - NLN, IRAD + NRN
C        .
C        . This is the column of the finite
C        . difference coefficients matrix
C        . which will contain the correct values.
C        .
         ICOCOL = INODE - IRAD + NLN + 1
C        .
C        . find the corresponding matrix column
C        .
         INDC = INDFUN( INODE, IHC, INARR )
C        .
C        . Find the actual location in the AMAT array
C        .
         CALL MATIND (INDR,INDC,IMF,KL,KU,KLE,N1,N2,IROW,ICOL)
C        .
C        .
C        .
         DO NDER = 0, IHD
           ICOROW = NDER + 1
           AMAT( IROW, ICOL ) = AMAT( IROW, ICOL ) + 
     1         FAC*COEFM( ICOROW, ICOCOL )*WORK( ICOROW )
         ENDDO
C        .
      ENDDO
C     .
      RETURN
      END
C*********************************************************************
