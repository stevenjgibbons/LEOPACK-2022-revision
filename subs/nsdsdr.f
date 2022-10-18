C*********************************************************************
C subroutine Number of Solid Body Singularities Determination Rout. **
C            -         -     -    -             -             -     **
C Steve Gibbons Wed Jul  5 13:42:12 BST 2000                         C
C____________________________________________________________________C
C                                                                    C
C If we are solving the momentum equation when the velocity has      C
C stress-free boundary conditions then we are likely to encounter    C
C (a) singularit(y)ies due to the fact that the system is not        C
C determined to within a solid body rotation about a given           C
C co-ordinate axis.                                                  C
C                                                                    C
C This routine looks at the harmonic sets present and the boundary   C
C conditions and determines the number of toroidal singularities,    C
C NTS, present in the system.                                        C
C                                                                    C
C All singularities are for toroidal velocity harmonics with l = 1.  C
C In all cases where axi-symmetric terms are present along with      C
C stress-free velocity boundary conditions, the t_1^0 harmonic is    C
C undetermined.                                                      C
C                                                                    C
C The t_1^{1s} and t_1^{1c} harmonics will only be undetermined if   C
C if the rotation vector is not explicitly included in the matrix.   C
C If the Coriols force is explicitly included the OCF logical flag   C
C is set to .TRUE. Otherwise, for instance when the matrix contains  C
C only diffusive terms, OCF = .FALSE.                                C
C                                                                    C
C NTSM, the maximum number of allowed toroidal singularities,        C
C defines the size of the integer array NUHA (number of undetermined C
C harmonics)                                                         C
C                                                                    C
C                                                                    C
C____________________________________________________________________C
C                                                                    C
C Input variables :-                                                 C
C ===============                                                    C
C  Integer                                                           C
C  -------                                                           C
C                                                                    C
C     NTSM      : Maximum number of toroidal singularities.          C
C     NTS       : Number of toroidal singularities.                  C
C     NUMA      : Numbers of undetermined harmonics. Dim (NTSM)      C
C     NH        : Number of harmonics.                               C
C     MHT       : Harmonic types, dim (NH).                          C
C     MHL       : Harmonic degrees, dim (NH).                        C
C     MHM       : Harmonic orders, dim (NH).                         C
C     MHP       : Pointers to finite difference scheme, dim (NH).    C
C     NDCS      : Maximum distinct finite difference schemes         C
C                  stored by the array SVFDC.                        C
C     MHIBC     : Inner boundary condition. Dim (NDCS).              C
C                                                                    C
C  MHIBC( is ) = 1 --> No condition is to be assumed at the bndry.   C
C  MHIBC( is ) = 2 --> Function must vanish at the bndry.            C
C  MHIBC( is ) = 3 --> First derivative must vanish at the bndry.    C
C  MHIBC( is ) = 4 --> Function must vanish at the bndry AND         C
C                       first derivative must vanish at the bndry.   C
C  MHIBC( is ) = 5 --> Function must vanish at the bndry AND         C
C                       second derivative must vanish at the bndry.  C
C  MHIBC( is ) = 6 --> r df/dr - f(r) = 0   at the bndry.            C
C  MHIBC( is ) = 7 --> r df/dr - l f(r) = 0 at the bndry.            C
C                        where L = LARR( is )                        C
C                                                                    C
C     MHOBC     : Outer boundary condition. Dim (NDCS).              C
C                                                                    C
C  MHOBC( is ) = 1 --> No condition is to be assumed at the bndry.   C
C  MHOBC( is ) = 2 --> Function must vanish at the bndry.            C
C  MHOBC( is ) = 3 --> First derivative must vanish at the bndry.    C
C  MHOBC( is ) = 4 --> Function must vanish at the bndry AND         C
C                       first derivative must vanish at the bndry.   C
C  MHOBC( is ) = 5 --> Function must vanish at the bndry AND         C
C                       second derivative must vanish at the bndry.  C
C  MHOBC( is ) = 6 --> r df/dr - f(r) = 0   at the bndry.            C
C  MHOBC( is ) = 7 --> r df/dr + (l+1) f(r) = 0 at the bndry.        C
C                        where L = LARR( is )                        C
C                                                                    C
C  Logical                                                           C
C  -------                                                           C
C                                                                    C
C     OCF       : Set to .TRUE. if matrix will included explicit     C
C                  Coriolis force terms. .FALSE. otherwise.          C
C                                                                    C
C____________________________________________________________________C
C
C*********************************************************************
      SUBROUTINE NSDSDR( NTSM, NTS, NUMA, NH, MHT, MHL, MHM, MHP,
     1                   NDCS, MHIBC, MHOBC, OCF )
      IMPLICIT NONE
C____________________________________________________________________C
C Variable declarations - Parameters ................................C
      INTEGER NTSM, NTS, NUMA( NTSM ), NH, MHT( NH ), MHL( NH ),
     1        MHM( NH ), MHP( NH ), NDCS, MHIBC( NDCS ),
     2        MHOBC( NDCS )
      LOGICAL OCF
C____________________________________________________________________C
C Variable declarations - Working variables .........................C
C
      INTEGER IH, IS
C____________________________________________________________________C
C START OF PROGRAM **************************************************C
C____________________________________________________________________C
C
      NTS = 0
C
      DO IH = 1, NH
C       .
        IF (   MHT( IH ).NE.2   ) GOTO 60
        IF (   MHL( IH ).NE.1   ) GOTO 60
C       .
C       . OK, so harmonic IH is toroidal velocity with l=1.
C       .
        IS = MHP( IH )
C       .
C       . Case of m=1 toroidal velocity harmonic with
C       . stress free boundaries
C       .
        IF ( MHIBC( IS ).EQ.6 .AND. MHOBC( IS ).EQ.6 .AND.
     1      ( MHM( IH ).EQ.0   .OR.   MHM( IH ).EQ.1   .OR.
     2        MHM( IH ).EQ.-1 )     ) THEN
          IF (  MHM( IH ).NE.0 .AND. OCF ) GOTO 60
          NTS = NTS + 1
          IF ( NTS.GT.NTSM ) THEN
            PRINT *,' Subroutine NSDSDR.'
            PRINT *,' NTS = ', NTS,' NTSM = ', NTSM
            PRINT *,' Program aborted.'
            STOP
          ENDIF
          NUMA( NTS ) = IH
          GOTO 60
        ENDIF
C       .
 60   CONTINUE
      ENDDO
C     .
      RETURN
      END
C*********************************************************************
