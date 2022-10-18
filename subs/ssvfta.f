C*********************************************************************
C subroutine Stored derivative Solution Vector Forcing Term Add ******
C            -                 -        -      -       -    -   ******
C Steve Gibbons Sun Feb 13 14:39:14 GMT 2000                         *
C____________________________________________________________________C
C                                                                    C
C Adds to a vector RHS the following terms in the MHD equations:     C
C                                                                    C
C F_Theta =  u . ( CB1 r + CB2 r^{-2} , 0 , 0 ) : SSVHST             C
C            - CC u . Grad ( Theta )            : SRVGTA             C
C                                                                    C
C F_vort. =     CH curl ( Theta )               : SSVTA              C
C             - CG curl ( K x v )               : SRKCVA             C
C             - CF curl ( v. Grad) v            : SRCVVA             C
C             + CJ curl ( B. Grad) B            : SRCVVA             C
C                                                                    C
C                                                                    C
C F_B     =   CM curl ( v x B )                 : SRMATA             C
C                                                                    C
C____________________________________________________________________C
C                                                                    C
C Input variables :-                                                 C
C ===============                                                    C
C  Integer                                                           C
C  -------                                                           C
C                                                                    C
C     INARR     : Format array for solution vectors. Dim (3).        C
C                                                                    C
C     MHT       : Dim ( NH ). See above.                             C
C     MHL       : Dim ( NH ). See above.                             C
C     MHM       : Dim ( NH ). See above.                             C
C                                                                    C
C     MHTR      : Dim ( NH ). MHT after calling CINDSW.              C
C                                                                    C
C     ISLARR    : Starting locations from MCNLIC. Dim ( 5 ).         C
C     NVIARR    : Number of coefficients from MCNLIC. Dim ( 5 ).     C
C                                                                    C
C     IHNA      : Number of alpha harmonics. Dim ( * )               C
C     IHNB      : Number of beta  harmonics. Dim ( * )               C
C     IHNG      : Number of gamma harmonics. Dim ( * )               C
C                                                                    C
C  Double Precision                                                  C
C  ----------------                                                  C
C                                                                    C
C     CB1       : Coefficient of v_r.r term in heat equation.        C
C     CB2       : Coefficient of v_r/(rr) term in heat equation.     C
C     CC        : Coefficient of v.Grad(Theta) term in heat eqn.     C
C     CH        : Coefficient of buoyancy term in vorticity eqn.     C
C     CG        : Coefficient of Coriolis term in vorticity eqn.     C
C     CF        : Coefficient of Inertial term in vorticity eqn.     C
C     CJ        : Coefficient of Lorentz  term in vorticity eqn.     C
C     CM        : Coefficient of Advection term in induction eqn.    C
C                                                                    C
C     V0        : Vector at time step i, zero^th derivatives.        C
C     V1        : Vector at time step i, first   derivatives.        C
C     V2        : Vector at time step i, second  derivatives.        C
C     V3        : Vector at time step i, third   derivatives.        C
C                                                                    C
C  (Pre-calculate these with CASVDR)                                 C
C                                                                    C
C     RHS       : Right hand side vector.                            C
C     XARR      : Dim ( NR ). i^{th} element gives x_i.              C
C     CVI       : Dim ( * ). Coefficients calculated by MCNLIC.      C
C                                                                    C
C  Character                                                         C
C  ---------                                                         C
C                                                                    C
C     TVHI      : *(4) Type of vector interaction. Dim. (MAXNVI).    C
C                 = 'CQSS', 'CQST' etc. according to the corresp.    C
C                 vector interaction.                                C
C                                                                    C
C____________________________________________________________________C
C
C*********************************************************************
      SUBROUTINE SSVFTA( INARR, MHT, MHL, MHM, MHTR, ISLARR, NVIARR,
     1                   IHNA, IHNB, IHNG, CB1, CB2, CC, CH, CG, CF,
     2                   CJ, CM, V0, V1, V2, V3, RHS, XARR, CVI, TVHI)
      IMPLICIT NONE
C____________________________________________________________________C
C Variable declarations - Parameters ................................C
      INTEGER INARR( * ), MHT( * ), MHL( * ), MHM( * ), MHTR( * ),
     1        ISLARR( 5 ), NVIARR( 5 ), IHNA( * ), IHNB( * ),
     2        IHNG( * )
      DOUBLE PRECISION CB1, CB2, CC, CH, CG, CF, CJ, CM, V0( * ),
     1                 V1( * ), V2( * ), V3( * ), RHS( * ),
     2                 XARR( * ), CVI( * )
      CHARACTER *(4) TVHI( * )
C____________________________________________________________________C
C Variable declarations - Working variables .........................C
C
      INTEGER NR, ILNR, IRNR, ILNT, IRNT, ISL, NVI, ICMPI, ICMPO
      DOUBLE PRECISION FAC
      CHARACTER *(3) CHVMFF
C____________________________________________________________________C
C START OF PROGRAM **************************************************C
C____________________________________________________________________C
C
      NR = INARR( 2 )
C     .
C     . First add on heat source terms
C     .
      ILNR  = 2
      IRNR  = NR - 1
      FAC   = 1.0d0
      CALL SSVHST( NR, V0, INARR, MHT, MHL, MHM, RHS,
     1             INARR, MHTR, MHL, MHM, FAC, ILNR, IRNR,
     2             XARR, CB1, CB2 )
C     .
C     . add on v . Grad( Theta ) terms
C     .
      NVI   = NVIARR( 1 )
      ISL   = ISLARR( 1 )
      ILNR  = 2
      IRNR  = NR - 1
      FAC   = -1.0d0*CC
      CALL SRVGTA( NVI, IHNA( ISL ), IHNB( ISL ), IHNG( ISL ),
     1             INARR, MHT, MHL, MHM, INARR, MHT, MHL, MHM,
     2             INARR, MHTR, MHL, MHM, ILNR, IRNR, NR,
     3             TVHI( ISL ), CVI( ISL ), V0, V1, V0, V1,
     4             RHS, XARR, FAC )
C     .
C     . add on Buoyancy terms
C     .
      ILNR  = 3
      IRNR  = NR - 2
      FAC   = CH
      ICMPI = 3
      ICMPO = 2
      CALL SSVTA( NR, V0, INARR, MHT, MHL, MHM, ICMPI,
     1            RHS, INARR, MHTR, MHL, MHM, ICMPO,
     2            FAC, ILNR, IRNR )
C     .
C     . add on Coriolis terms
C     .
      NVI   = NVIARR( 2 )
      ISL   = ISLARR( 2 )
      ILNR  = 2
      IRNR  = NR - 1
      ILNT  = 3
      IRNT  = NR - 2
      FAC   = -1.0d0*CG
      CALL SRKCVA( NVI, IHNA( ISL ), IHNG( ISL ), INARR, MHT,
     1             MHL, MHM, INARR, MHTR, MHL, MHM, ILNR, IRNR,
     2             ILNT, IRNT, NR, TVHI( ISL ), CVI( ISL ), V0,
     3             V1, V2, XARR, FAC, RHS)
C     .
C     . add on inertial term
C     .
      NVI    = NVIARR( 3 )
      ISL    = ISLARR( 3 )
      ILNR   = 2
      IRNR   = NR - 1
      ILNT   = 3
      IRNT   = NR - 2
      FAC    = -1.0d0*CF
      CHVMFF = 'VEL'
      CALL SRCVVA( NVI, IHNA( ISL ), IHNB( ISL ), IHNG( ISL ),
     1             INARR, MHT, MHL, MHM, INARR, MHT, MHL, MHM,
     2             INARR, MHTR, MHL, MHM, ILNR, IRNR, ILNT, IRNT,
     3             NR, TVHI( ISL ), CVI( ISL ), V0, V1, V2, V0,
     4             V1, V2, V3, RHS, XARR, FAC, CHVMFF)
C     .
C     . add on Lorentz term
C     .
      NVI    = NVIARR( 4 )
      ISL    = ISLARR( 4 )
      ILNR   = 2
      IRNR   = NR - 1
      ILNT   = 3
      IRNT   = NR - 2
      FAC    = CJ
      CHVMFF = 'MAG'
      CALL SRCVVA( NVI, IHNA( ISL ), IHNB( ISL ), IHNG( ISL ),
     1             INARR, MHT, MHL, MHM, INARR, MHT, MHL, MHM,
     2             INARR, MHTR, MHL, MHM, ILNR, IRNR, ILNT, IRNT,
     3             NR, TVHI( ISL ), CVI( ISL ), V0, V1, V2, V0,
     4             V1, V2, V3, RHS, XARR, FAC, CHVMFF)
C     .
C     . add on Induction term
C     .
      NVI    = NVIARR( 5 )
      ISL    = ISLARR( 5 )
      ILNR   = 2
      IRNR   = NR - 1
      ILNT   = 3
      IRNT   = NR - 2
      FAC    = CM
      CALL SRMATA( NVI, IHNA( ISL ), IHNB( ISL ), IHNG( ISL ),
     1             INARR, MHT, MHL, MHM, INARR, MHT, MHL, MHM,
     2             INARR, MHTR, MHL, MHM, ILNR, IRNR, ILNT, IRNT,
     3             NR, TVHI( ISL ), CVI( ISL ), V0, V1, V2, V0,
     4             V1, V2, RHS, XARR, FAC )
C     .
      RETURN
      END
C*********************************************************************
