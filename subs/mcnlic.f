C*********************************************************************
C subroutine Magnetic Convection Non-Linear Interaction Calculate ****
C            -        -          -   -      -           -         ****
C Steve Gibbons Fri Feb 11 13:52:06 GMT 2000                         *
C____________________________________________________________________C
C                                                                    C
C Given a harmonic set defined by MHT, MHL, MHM, MHTR (see CINDSW)   C
C for a single vector containing Velocity, Temperature AND Magnetic  C
C Field (although the routine is quite valid without mag. field),    C
C MCNLIC will calculate the interactions for the terms               C
C                                                                    C
C  v . Grad( THETA )        (using VSPCC)                            C
C                                                                    C
C  curl ( k x v )           (using CFVICC)                           C
C                                                                    C
C  curl ( v . Grad ) v      (using VCCPCG)                           C
C                                                                    C
C  curl ( B . Grad ) B      (using VCCPCG)                           C
C                                                                    C
C  curl ( v x B )           (using VCPCC)                            C
C                                                                    C
C____________________________________________________________________C
C                                                                    C
C Input variables :-                                                 C
C ===============                                                    C
C  Integer                                                           C
C  -------                                                           C
C                                                                    C
C     INARR    : Format of the solution vector. (see INDFUN)         C
C                                                                    C
C     MHT      : Type of harmonics (Standard system)                 C
C     MHL      : Sph. harm degree, l. (Standard system)              C
C     MHM      : Order, m or (-m) (Standard system)                  C
C     MHTR     : Type of curl harmonics. (Standard system)           C
C                                                                    C
C     LH        : Maximum spherical harmonic degree, l.              C
C     NTHP      : The number of theta points.                        C
C     NPHP      : The number of phi points.                          C
C     MMAX      : Maximum sph. harmonic order, m.                    C
C                                                                    C
C     IHNA     : Number of alpha harmonics. Dim ( MAXNVI )           C
C     IHNB     : Number of beta harmonics. Dim ( MAXNVI )            C
C     IHNG     : Number of gamma harmonics. Dim ( MAXNVI )           C
C                                                                    C
C     MAXNVI   : Maximum number of vector interactions.              C
C                                                                    C
C     ISLARR   : Dimension (5). Starting element for different       C
C                categories of vector interaction.                   C
C                                                                    C
C          ISLARR( 1 ): First element for    v . Grad( THETA )       C
C          ISLARR( 2 ): First element for    curl ( k x v )          C
C          ISLARR( 3 ): First element for    curl ( v . Grad ) v     C
C          ISLARR( 4 ): First element for    curl ( B . Grad ) B     C
C          ISLARR( 5 ): First element for    curl ( v x B )          C
C                                                                    C
C     NVIARR   : Dimension (5). Number of interactions for different C
C                categories of vector interaction.                   C
C                                                                    C
C          NVIARR( 1 ): Number of elmnts for   v . Grad( THETA )     C
C          NVIARR( 2 ): Number of elmnts for   curl ( k x v )        C
C          NVIARR( 3 ): Number of elmnts for   curl ( v . Grad ) v   C
C          NVIARR( 4 ): Number of elmnts for   curl ( B . Grad ) B   C
C          NVIARR( 5 ): Number of elmnts for   curl ( v x B )        C
C                                                                    C
C                                                                    C
C  Character                                                         C
C  ---------                                                         C
C                                                                    C
C     TVHI      : *(4) Type of vector interaction. Dim. (MAXNVI).    C
C                 = 'CQSS', 'CQST' etc. according to the corresp.    C
C                 vector interaction.                                C
C                                                                    C
C  Double Precision                                                  C
C  ----------------                                                  C
C                                                                    C
C     GAUX      : Cosines of the NTHP evaluated by the routine       C
C                  gauwts. Dimension ( NTHP ).                       C
C     GAUW      : Corresponding weights. As above.                   C
C     PA        : Schmidt Normalised Legendre Functions              C
C                  Dimension ( ( LH + 1 )*( LH + 2 )/2, NTHP )       C
C                  P_l^m(X_j) = PA( L*(L+1)/2 + M + 1 , j )          C
C     DPA       : Derivatives of the above.                          C
C                                                                    C
C     FTF1      : Work array - dim. (2*NPHP)                         C
C     FTF2      : Work array - dim. (2*NPHP)                         C
C     FTF3      : Work array - dim. (2*NPHP)                         C
C     VF1       : Work array - dim. ( NPHP, NTHP, 3)                 C
C     VF2       : Work array - dim. ( NPHP, NTHP, 3)                 C
C     VF3       : Work array - dim. ( NPHP, NTHP, 3)                 C
C     QST       : Work array - dim. ( LH * ( LH + 2 ), 3 )           C
C     SF        : Work array - dim. ( NPHP, NTHP )                   C
C     SHC       : Work array - dim. ( LH * ( LH + 2 ) )              C
C                                                                    C
C     CVI       : Coefficient of vector interaction. Dim ( MAXNVI )  C
C                                                                    C
C____________________________________________________________________C
C
C*********************************************************************
      SUBROUTINE MCNLIC( INARR, MHT, MHL, MHM, MHTR, LH, NTHP, NPHP,
     1       MMAX, GAUX, GAUW, PA, DPA, FTF1, FTF2, FTF3, VF1, VF2,
     2       VF3, QST, SF, SHC, IHNA, IHNB, IHNG, CVI, TVHI, MAXNVI,
     3       ISLARR, NVIARR )
      IMPLICIT NONE
C____________________________________________________________________C
C Variable declarations - Parameters ................................C
      INTEGER INARR( * ), MHT( * ), MHL( * ), MHM( * ), MHTR( * ),
     1        LH, NTHP, NPHP, MMAX, IHNA( * ), IHNB( * ), IHNG( * ),
     2        MAXNVI, ISLARR( 5 ), NVIARR( 5 )
      CHARACTER *(4) TVHI( * )
      DOUBLE PRECISION VF1( NPHP, NTHP, 3), VF2( NPHP, NTHP, 3),
     1                 SF( NPHP, NTHP ), SHC( LH*( LH + 2) )
      DOUBLE PRECISION 
     1                 VF3( NPHP, NTHP, 3), FTF1( 2*NPHP ),
     2                 FTF2( 2*NPHP ), FTF3( 2*NPHP ),
     3                 CVI( * ), QST( LH*(LH+2), 3)
      DOUBLE PRECISION
     1                 GAUX( NTHP ), GAUW( NTHP ),
     2                 PA ( ( LH + 1 )*( LH + 2 )/2 , NTHP),
     3                 DPA ( ( LH + 1 )*( LH + 2 )/2 , NTHP)
C____________________________________________________________________C
C Variable declarations - Working variables .........................C
C
      INTEGER NVI, ISL, NH, MAXH
      CHARACTER *(3) CHVMFF
C____________________________________________________________________C
C START OF PROGRAM **************************************************C
C____________________________________________________________________C
C     .
      NH = INARR( 3 )
C     .
C     . First do v . Grad( THETA ) terms
C     . this uses routine VSPCC ...
C     .
      ISL    = 1
      MAXH   = MAXNVI
      NVI    = 0
      CALL VSPCC( NVI, MAXH, IHNA( ISL ), IHNB( ISL ), IHNG( ISL ),
     1            NH, MHT, MHL, MHM, NH, MHT, MHL, MHM, NH, MHTR,
     2            MHL, MHM, LH, NTHP, NPHP, MMAX, TVHI( ISL ),
     3            CVI( ISL ), GAUX, GAUW, PA, DPA, FTF1, VF1, VF2,
     4            SF, SHC )
C     .
C     . Register outcome of VSPCC
C     .
      ISLARR( 1 ) = ISL
      NVIARR( 1 ) = NVI
C     .
C     . Now do curl ( k x v ) terms
C     . this uses routine CFVICC ...
C     .
      ISL    = ISL  + NVI
      MAXH   = MAXH - NVI
      NVI    = 0
      CALL CFVICC( NVI, MAXH, IHNA( ISL ), IHNG( ISL ), NH, MHT, MHL,
     1             MHM, MHTR, MHL, MHM, LH, NTHP, NPHP, MMAX, 
     2             TVHI( ISL ), CVI( ISL ), GAUX, GAUW, PA, DPA,
     3             FTF1, FTF2, FTF3, VF1, QST )
C     .
C     . Register outcome of CFVICC
C     .
      ISLARR( 2 ) = ISL
      NVIARR( 2 ) = NVI
C     .
C     . Now do curl ( v . Grad ) v terms
C     . this uses routine VCCPCG ...
C     .
      ISL    = ISL  + NVI
      MAXH   = MAXH - NVI
      NVI    = 0
      CHVMFF = 'VEL'
      CALL VCCPCG( NVI, MAXH, IHNA( ISL ), IHNB( ISL ), IHNG( ISL ),
     1     NH, MHT, MHL, MHM, NH, MHT, MHL, MHM, NH, MHTR, MHL, MHM,
     2     LH, NTHP, NPHP, MMAX, TVHI( ISL ), CVI( ISL ), GAUX, GAUW,
     3     PA, DPA, FTF1, FTF2, FTF3, VF1, VF2, VF3, QST, CHVMFF )
C     .
C     . Register outcome of VCCPCG
C     .
      ISLARR( 3 ) = ISL
      NVIARR( 3 ) = NVI
C     .
C     . Now do curl ( B . Grad ) B terms
C     . this uses routine VCCPCG ...
C     .
      ISL    = ISL  + NVI
      MAXH   = MAXH - NVI
      NVI    = 0
      CHVMFF = 'MAG'
      CALL VCCPCG( NVI, MAXH, IHNA( ISL ), IHNB( ISL ), IHNG( ISL ),
     1     NH, MHT, MHL, MHM, NH, MHT, MHL, MHM, NH, MHTR, MHL, MHM,
     2     LH, NTHP, NPHP, MMAX, TVHI( ISL ), CVI( ISL ), GAUX, GAUW,
     3     PA, DPA, FTF1, FTF2, FTF3, VF1, VF2, VF3, QST, CHVMFF )
C     .
C     . Register outcome of VCCPCG
C     .
      ISLARR( 4 ) = ISL
      NVIARR( 4 ) = NVI
C     .
C     . Now do curl ( u x B ) terms
C     . this uses routine VCPCC ...
C     .
      ISL    = ISL  + NVI
      MAXH   = MAXH - NVI
      NVI    = 0
      CALL VCPCC( NVI, MAXH, IHNA( ISL ), IHNB( ISL ), IHNG( ISL ),
     1            NH, MHT, MHL, MHM, NH, MHT, MHL, MHM, NH, MHTR,
     2            MHL, MHM, LH, NTHP, NPHP, MMAX, TVHI( ISL ),
     3            CVI( ISL ), GAUX, GAUW, PA, DPA, FTF1, FTF2, FTF3,
     4            VF1, VF2, VF3, QST )
C     .
C     . Register outcome of VCCPCG
C     .
      ISLARR( 5 ) = ISL
      NVIARR( 5 ) = NVI
C     .
      RETURN
      END
C*********************************************************************
