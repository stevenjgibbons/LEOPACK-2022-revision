C*********************************************************************
C subroutine Fixed Velocity Timestep Growth Rate Find ****************
C            -     -        -        -      -    -    ****************
C Steve Gibbons Thu May 25 15:15:50 BST 2000                         C
C____________________________________________________________________C
C                                                                    C
C Takes an initial vector VEC (non-zero) and time-steps linearly     C
C with the equations below (linear in the velocity, non-linear terms C
C only with imposed velocity field, v_0) until either MNTS timesteps C
C are taken, or the growth rate for all wavenumbers has converged    C
C such that the rate of change in the growth rate is less than       C
C GRT for NTSC timesteps.                                            C
C The resulting growth rates are returned in SIGM.                   C
C                                                                    C
C____________________________________________________________________C
C                                                                    C
C  Let F_Theta and F_Vort the forcing terms in the                   C
C  heat and vorticity equations respectively with                    C
C                                                                    C
C  F_Theta = v ( CB1 r + CB2 r^{-2} , 0, 0 ) - CC v_0 . Grad (Theta) C
C                                                                    C
C  F_Vort  =   CH curl ( Theta r ) - CG curl ( k x v )               C
C             - CF curl ( v . Grad ) v_0                             C
C             - CF curl ( v_0 . Grad ) v                             C
C                                                                    C
C____________________________________________________________________C
C                                                                    C
C  Let v^{i} and Theta^{i} be velocity and temperature               C
C respectively at time step i.                                       C
C                                                                    C
C  Using MCDRHF and SSVFTA, we form the forcing terms                C
C                                                                    C
C      F1 =  [CA+DELTAT.CD.CFAC Lap] Theta^{i} +                     C
C                                DELTAT. F_Theta( v, The ),          C
C                                                                    C
C      F2 =  [CE+DELTAT.CI.CFAC Lap] curl v^{i} +                    C
C                                DELTAT.F_Vort( v, The ),            C
C                                                                    C
C A predictor ( Theta_(0)^{i+1}, v_(0)^{i+1} )                       C
C is then formed by solving                                          C
C                                                                    C
C [ CA + DELTAT.CD.(CFAC-1) Lap] Theta_(0)^{i+1} = F1                C
C [ CE + DELTAT.CI.(CFAC-1) Lap] v_(0)^{i+1} = F2                    C
C                                                                    C
C and then functional iteration is used to semi-implicitly solve     C
C the equations                                                      C
C                                                                    C
C [ CA + DELTAT.CD.(CFAC-1) Lap] Theta_(j)^{i+1} = F1                C
C [ CE + DELTAT.CI.(CFAC-1) Lap] v_(j)^{i+1} = F2                    C
C                                                                    C
C  where                                                             C
C                                                                    C
C F1 = [CA+DELTAT.CD.CFAC Lap] Theta^{i}                             C
C      + CFAC.DELTAT.F_Theta( v,The^{i} )                            C
C      + (1-CFAC).DELTAT.F_Theta( v,The_{j-1}^{i+1} )                C
C                                                                    C
C F2 = [CE+DELTAT.CI.CFAC Lap] Vort^{i}                              C
C      + CFAC.DELTAT.F_Vort( v,The^{i} )                             C
C      + (1-CFAC).DELTAT.F_Vort( v,The_{j-1}^{i+1} )                 C
C                                                                    C
C This functional iteration continues until                          C
C                                                                    C
C DABS( |V_(j)^{i+1}| - |V_{j-1}^{i+1}| ) < TOL                      C
C                                                                    C
C where V represents the full solution vector.                       C
C                                                                    C
C INARR( 1 ) for VEC = 4.                                            C
C____________________________________________________________________C
C                                                                    C
C Input variables :-                                                 C
C ===============                                                    C
C  Integer                                                           C
C  -------                                                           C
C                                                                    C
C     KL        : Number of diagonals in matrix.                     C
C     NBN       : Number of bounding nodes.                          C
C     NCFM      : Leading dimension of SVFDC and FDCM arrays.        C
C     NR        : Number of radial grid nodes.                       C
C     NDRVM     : Max deriv.s stored in SVFDC array.                 C
C     NH        : Number of spherical harmonics in solution vec.     C
C     NDCS      : Number of finite difference schemes.               C
C                                                                    C
C     MHT       : Dim ( * ). Harmonic type.                          C
C     MHL       : Dim ( * ). Harmonic degree, l.                     C
C     MHM       : Dim ( * ). Harmonic order, m.                      C
C     MHP       : Pointer array for harmonics. If HMP( ih ) = is     C
C                  then 'is' is the finite difference scheme used    C
C                   to take derivatives of that harm. radial func.   C
C                   If MHP is negative, the harmonic is avoided and  C
C                   CASVDR moves on to the next harmonic.            C
C     MHTR      : Dim ( * ). Harmonic type for curl.                 C
C     LH        : Maximum spherical harmonic degree, l.              C
C     MMAX      : Maximum value of m.                                C
C     M0        : Minimum value of non-zero m such that              C
C                  MOD( m, M0 ) = 0 for all m.                       C
C     N2        : Dimension of vectors. (2nd dim of DMAT).           C
C                                                                    C
C     NOIT      : On input. NOIT = Maximum number of iterations      C
C                 permitted in the functional iteration part of the  C
C                 process in the routine.                            C
C                 On output. NOIT is number of iterations taken      C
C                 UNLESS - maximum number is exceeded (NOIT = -1)    C
C                 OR the iteration norm,                             C
C                 DABS( |V_(j)^{i+1}| - |V_{j-1}^{i+1}| )            C
C                 is actually growing (NOIT = -2).                   C
C                                                                    C
C     LULOG     : Set to zero to suppress output of information.     C
C                 Otherwise, set to output chanel wth this number.   C
C                                                                    C
C     IPIV      : Dim (N2). Output from SFDDMF.                      C
C                                                                    C
C     MHIBC     : MHIBC( is ) describes the inner boundary           C
C                  condition for scheme IS.                          C
C                                                                    C
C  MHIBC( is ) = 1 --> No condition is to be assumed at the bndry.   C
C  MHIBC( is ) = 2 --> Function must vanish at the bndry.            C
C  MHIBC( is ) = 3 --> First derivative must vanish at the bndry.    C
C  MHIBC( is ) = 4 --> Function must vanish at the bndry AND         C
C                       first derivative must vanish at the bndry.   C
C  MHIBC( is ) = 5 --> Function must vanish at the bndry AND         C
C                       second derivative must vanish at the bndry.  C
C  MHIBC( is ) = 6 --> r df/dr - f(r) = 0   at the bndry.            C
C  MHIBC( is ) = 7 --> insulating magnetic field.                    C
C                                                                    C
C     MHOBC     : MHOBC( is ) describes the outer boundary           C
C                  condition for scheme IS. (See above for key.)     C
C                                                                    C
C     NTS       : Number of toroidal singularities.                  C
C                 (See SFDDMF, DMWMF )                               C
C                                                                    C
C     NNDS      : Number of nodes (between 2 and NR)                 C
C     IW1       : Dim (NNDS). Work array.                            C
C                                                                    C
C     NWM       : Number of wavenumbers present.                     C
C     MNTS      : Maximum number of timesteps permitted.             C
C                                                                    C
C     IEF       : On output contains -1 if max number of time-steps  C
C                 was exceeded. On a successful output, IEF is the   C
C                 number of time-steps taken.                        C
C                 If -2 on exit, FVITSR has failed.                  C
C                                                                    C
C                                                                    C
C     MA        : Dim (NWM ). Actual wavenumbers.                    C
C     NTSC      : Number of timesteps for which convergence must     C
C                   be maintained.                                   C
C                                                                    C
C  Double Precision                                                  C
C  ----------------                                                  C
C                                                                    C
C     VEC       : Dim(*). Solution vector at time step (i)           C
C     VEC1      : Dim(*). Solution vector at time step (i+1)         C
C     SVFDC     : Finite difference coefficient matrix.              C
C                  Dimension ( NFDCM, NR, NDRVM+1, NDCS ).           C
C                   Array is generated by the routine svfdcf         C
C                 See documentation for SVFDCF for details.          C
C                                                                    C
C     FDCM      : Finite difference coefficient matrix.              C
C                  Dimension ( NCFM, NR, 1 ).                        C
C                   Array is generated by the routine fdcmbd         C
C                 See documentation for FDCMBD for details.          C
C       MUST be calculated with:                                     C
C                                                                    C
C                     NDRVM = 1                                      C
C                     NLMN  = 2                                      C
C                     NRMN  = NR - 1                                 C
C                     NLMC  = 2                                      C
C                     NRMC  = NR - 1                                 C
C                                                                    C
C     PA        : Schmidt Normalised Legendre Functions              C
C                  Dim { ( LH + 1 )*( LH + 2 )/2 , NTHP }            C
C                  P_l^m(X_j) = PA( L*(L+1)/2 + M + 1 , j )          C
C     DPA       : Derivatives of the above.                          C
C     GAUX      : Cosines of the NTHP evaluated by the routine       C
C                  gauwts. Dimension ( NTHP ).                       C
C     GAUW      : Gauss weights computed bu GAUWTS. Dim ( NTHP )     C
C     SHC       : Work array. Dim ( LH*(LH+2) ).                     C
C     DSHC      : Work array. Dim ( LH*(LH+2) ).                     C
C     RVF1      : Work arr. dim. ( NR, NPHP, NTHP, 3 )               C
C     RVF2      : Work arr. dim. ( NR, NPHP, NTHP, 3 )               C
C     RVF3      : Work arr. dim. ( NR, NPHP, NTHP, 3 )               C
C     SF        : Work arr. dim. ( NPHP, NTHP )                      C
C     VF        : Work arr. dim. ( NPHP, NTHP, 3 )                   C
C                                                                    C
C     XARR      : Dim ( NR ). i^{th} element gives x_i.              C
C     DMAT      : Diffusion matrix. Dim ( 3*KL + 1, N2 ).            C
C     DPARS     : Dim( * ). Supplies double precision parameters.    C
C                                                                    C
C                  DPARS(  1 ) = DELTAT                              C
C                  DPARS(  2 ) = TOL                                 C
C                  DPARS(  3 ) = CFAC                                C
C                  DPARS(  4 ) = CA                                  C
C                  DPARS(  5 ) = CB1                                 C
C                  DPARS(  6 ) = CB2                                 C
C                  DPARS(  7 ) = CC                                  C
C                  DPARS(  8 ) = CD                                  C
C                  DPARS(  9 ) = CE                                  C
C                  DPARS( 10 ) = CF                                  C
C                  DPARS( 11 ) = CG                                  C
C                  DPARS( 12 ) = CH                                  C
C                  DPARS( 13 ) = CI                                  C
C                  DPARS( 14 ) = (not used)                          C
C                  DPARS( 15 ) = CK (probably not used)              C
C                  DPARS( 16 ) = CL (probably not used)              C
C                  DPARS( 17 ) = (not used)                          C
C                                                                    C
C     V0        : Dim ( * ). Work array.                             C
C     V1        : Dim ( * ). Work array.                             C
C     V2        : Dim ( * ). Work array.                             C
C     V3        : Dim ( * ). Work array.                             C
C     V4        : Dim ( * ). Work array.                             C
C     DTV       : Dim ( * ). Work array.                             C
C     RFI       : Dim ( * ). Work array.                             C
C     RFI1      : Dim ( * ). Work array.                             C
C                                                                    C
C     RQST1     : Work array. Dim (  LH*(LH+2) ,3, NR ).             C
C     RQST2     : Work array. Dim (  LH*(LH+2) ,3, NR ).             C
C     RQST3     : Work array. Dim (  LH*(LH+2) ,3, NR ).             C
C     RQSTA     : Work array. Dim (  LH*(LH+2) ,3, NR ).             C
C                                                                    C
C     W1        : Work array. Dim ( NNDS )                           C
C     W2        : Work array. Dim ( NNDS )                           C
C     W3        : Work array. Dim ( NNDS, NNDS )                     C
C     U         : U matrix :  Dim ( N2, NTS )                        C
C     V         : V matrix :  Dim ( N2, NTS )                        C
C                                                                    C
C     V0RVF     : Fixed velocity, v_0. Dim( NR, NPHP, NTHP, 3 )      C
C     CV0RVF    : Curl of v_0. Dim( NR, NPHP, NTHP, 3 )              C
C                                                                    C
C     SIGM      : Dim (NWM). Array of growth rates.                  C
C     GRT       : Tolerance for growth rate convergence.             C
C____________________________________________________________________C
C
C*********************************************************************
      SUBROUTINE FVTGRF( KL,NBN,NCFM,NR,NDRVM,NH,NDCS,MHT,MHL,MHM,MHP,
     1   MHTR,LH,MMAX,M0,NPHP,NTHP,VEC,SVFDC,FDCM,PA,DPA,GAUX,GAUW,
     2   SHC,DSHC,RVF1,RVF2,RVF3,SF,VF,XARR,N2,NOIT,LULOG,IPIV,MHIBC,
     3   MHOBC,DMAT,DPARS,V0,V1,V2,V3,V4,DTV,RFI,RFI1,VEC1,RQST1,
     4   RQST2,RQST3,RQSTA,NTS,NNDS,IW1,W1,W2,W3,U,V,V0RVF,CV0RVF,NWM,
     5   MNTS,SIGM,IEF,MA,GRT,NTSC )
      IMPLICIT NONE
C____________________________________________________________________C
C Variable declarations - Parameters ................................C
      INTEGER KL, NBN, NCFM, NR, NDRVM, NH, NDCS, MHT( * ), MHL( * ),
     1        MHM( * ), MHP( * ), MHTR( * ), LH, MMAX, M0, NPHP, NTHP,
     2        N2, NOIT, LULOG, IPIV( * ), MHIBC( * ), MHOBC( * ), NTS,
     3        NNDS, IW1( NNDS ), NWM, MNTS, IEF, MA( NWM ), NTSC
      DOUBLE PRECISION VEC( * ), SVFDC( NCFM, NR, NDRVM+1, NDCS ),
     1                 XARR( * ), DMAT( 3*KL + 1, N2 ),
     2                 DPARS( * ), V0( * ), V1( * ), V2( * ),
     3                 V3( * ), V4( * ), DTV( * ), RFI( * ),
     2                 RFI1( * ), VEC1( * )
      DOUBLE PRECISION V0RVF( NR, NPHP, NTHP, 3 ),
     1                 CV0RVF( NR, NPHP, NTHP, 3 ),
     2                 SIGM( NWM ), GRT
      DOUBLE PRECISION FDCM( NCFM, NR, 1 ), SF( NPHP, NTHP ),
     1                 RVF1( NR, NPHP, NTHP, 3 ), GAUX( NTHP ),
     2                 RVF2( NR, NPHP, NTHP, 3 ), GAUW( NTHP ),
     3                 RVF3( NR, NPHP, NTHP, 3 ), VF( NPHP, NTHP, 3 )
      DOUBLE PRECISION SHC( LH*(LH+2) ), DSHC( LH*(LH+2) ),
     1                 PA ( ( LH + 1 )*( LH + 2 )/2 , NTHP),
     2                 DPA ( ( LH + 1 )*( LH + 2 )/2 , NTHP)
      DOUBLE PRECISION RQST1( LH*(LH+2), 3, NR ),
     1                 RQST2( LH*(LH+2), 3, NR ),
     2                 RQST3( LH*(LH+2), 3, NR ),
     3                 RQSTA( LH*(LH+2), 3, NR )
      DOUBLE PRECISION W1( NNDS ), W2( NNDS ), W3( NNDS, NNDS ),
     1                 U( N2, NTS ), V( N2, NTS )
C____________________________________________________________________C
C Variable declarations - Working variables .........................C
      INTEGER NWMAX
      PARAMETER ( NWMAX = 20 )
      DOUBLE PRECISION WORKAR( 5, NWMAX ), DLOW, DKE( 2 ), FM, FP,
     1                 DELTAT, DTIME
      PARAMETER ( DLOW = 1.0d-8 )
      INTEGER NTSARR( NWMAX ), ITS, MM, IMM, NCW, INARR( 3 ),
     1        IH, NOIT2
C
C Explanation of work arrays.
C ===========================
C
C WORKAR( 1, IMM ) contains sqrt( ke ) in wavenumber MA( imm )
C                  at the previous time step.
C WORKAR( 2, IMM ) contains sqrt( ke ) in wavenumber MA( imm )
C                  at the current time step.
C WORKAR( 3, IMM ) contains growth rate for wavenumber MA( imm )
C                  at the previous time step.
C WORKAR( 4, IMM ) contains growth rate for wavenumber MA( imm )
C                  at the current time step.
C
C this is given by
C
C                     2 * [ WORKAR( 2, IMM ) - WORKAR( 1, IMM ) ]
C WORKAR( 4, IMM ) =  --------------------------------------------
C                    DELTAT * [WORKAR( 2, IMM ) + WORKAR( 1, IMM )]
C
C
C WORKAR( 5, IMM ) contains the rate of change of the growth rate
C  with
C
C                    
C WORKAR( 5, IMM ) =  WORKAR( 4, IMM ) - WORKAR( 3, IMM )
C
C NTSARR( imm ) is the number of time-steps where the 
C                growth rate qualifies.
C NCW = number of converged wavenumbers
C____________________________________________________________________C
C START OF PROGRAM **************************************************C
C____________________________________________________________________C
C
C Check that we have sufficient wavenumbers ...
C
      IF ( NWM.GT.NWMAX ) THEN
        PRINT *,' Subroutine FVTGRF.'
        PRINT *,' NWM = ', NWM,' NWMAX = ', NWMAX
        PRINT *,' Program aborted.'
        STOP
      ENDIF
C
      ITS   = 0
      DTIME = 0.0d0
C
C Check size of time-step
C
      DELTAT = DPARS( 1 )
      IF ( DELTAT.LT.DLOW ) THEN
        PRINT *,' Subroutine FVTGRF.'
        PRINT *,' DELTAT = ', DELTAT
        PRINT *,' Division by zero imminent.'
        PRINT *,' Program aborted.'
        STOP
      ENDIF
C
      INARR( 1 ) = 4
      INARR( 2 ) = NR
      INARR( 3 ) = NH
C
C Initialise key arrays to 'harmless' values
C
      DO IMM = 1, NWM
        WORKAR( 4, IMM ) = 400.0d0
        WORKAR( 5, IMM ) = 500.0d0
        NTSARR( IMM )    = 0
      ENDDO
C
C Start looping around time-steps
C
 50   CONTINUE
      ITS   = ITS + 1
      DTIME = DTIME + DELTAT
C
      DO IMM = 1, NWM
        WORKAR( 1, IMM ) = WORKAR( 2, IMM )
        WORKAR( 2, IMM ) = 0.0d0
        WORKAR( 3, IMM ) = WORKAR( 4, IMM )
      ENDDO
C
      NOIT2 = NOIT
      CALL FVITSR( KL, NBN, NCFM, NR, NDRVM, NH, NDCS, MHT, MHL,
     1   MHM, MHP, MHTR, LH, MMAX, M0, NPHP, NTHP, VEC, SVFDC, FDCM,
     2   PA, DPA, GAUX, GAUW, SHC, DSHC, RVF1, RVF2, RVF3, SF, VF,
     3   XARR, N2, NOIT2, LULOG, IPIV, MHIBC, MHOBC, DMAT, DPARS, V0,
     4   V1, V2, V3, V4, DTV, RFI, RFI1, VEC1, RQST1, RQST2, RQST3,
     5   RQSTA, NTS, NNDS, IW1, W1, W2, W3, U, V, V0RVF, CV0RVF )
C
C Return if too many iterations were taken
C
      IF ( NOIT2.LT.0 ) THEN
        IEF = -2
        RETURN
      ENDIF
C
C Calculate the kinetic energies of different components
C
      DO IH = 1, NH
        CALL SHKEER( IH, NDCS, NR, INARR, MHT, MHL, MHP, NBN,
     1             NDRVM, NDRVM, NCFM, VEC1, XARR, DKE, SVFDC )
        MM = MHM( IH )
        MM = IABS( MM )
        DO IMM = 1, NWM
          IF ( MA( IMM ).EQ.MM ) THEN
            FM = DKE( 1 )
            WORKAR( 2, IMM ) = WORKAR( 2, IMM ) + FM
          ENDIF
        ENDDO
      ENDDO
C
      NCW = 0
C
C Calculate growth rates
C
      DO IMM = 1, NWM
C       .
C       . First calculate growth rate at current time-step
C       .
        WORKAR( 2, IMM ) = DSQRT( WORKAR( 2, IMM ) )
        FP = WORKAR( 2, IMM ) + WORKAR( 1, IMM )
        FM = WORKAR( 2, IMM ) - WORKAR( 1, IMM )
        IF ( DABS( FP ).LT.DLOW ) THEN
          PRINT *,' Subroutine FVTGRF.'
          PRINT *,' FP = ', FP
          PRINT *,' Division by zero imminent.'
          PRINT *,' Program aborted.'
          STOP
        ENDIF
        WORKAR( 4, IMM ) = 2.0d0*FM/(DELTAT*FP)
        SIGM( IMM ) = WORKAR( 4, IMM )
C       .
C       . Calculate rate of change of growth rate
C       .
        FM = WORKAR( 4, IMM ) - WORKAR( 3, IMM )
        WORKAR( 5, IMM ) = DABS( FM )/DELTAT
C       .
        IF ( DABS( WORKAR( 5, IMM ) ).LT.GRT ) THEN
          NTSARR( IMM ) = NTSARR( IMM ) + 1
        ELSE
          NTSARR( IMM ) = 0
        ENDIF
C       .
        IF ( NTSARR( IMM ).GE.NTSC ) THEN
          NCW = NCW + 1
        ENDIF
C       .
        IF ( LULOG.NE.0 ) WRITE ( LULOG, 812 ) DTIME, MA( IMM ),
     1          WORKAR( 2, IMM ), WORKAR( 4, IMM ), WORKAR( 5, IMM )
 812    FORMAT('fvtgrf: t= ',1PD12.4,' m= ',I4,' k= ',1PD12.4,' s= ',
     1                 1PD12.4,' c= ',1PD12.4)
C       .
      ENDDO
C
C All of our growth rates have converged
C
      IF ( NCW.EQ.NWM ) THEN
        IEF = ITS
        RETURN
      ENDIF
C
C We have not converged and so return with error
C if we have done max. number of time-steps 
C
      IF ( ITS.EQ.MNTS ) THEN
        IEF = -1
        RETURN
      ENDIF
C
C Do next time-step
C Copy VEC1 into VEC:
C
      CALL VECCP( VEC1, VEC, N2 )
C
      GOTO 50
C
      END
C*********************************************************************

