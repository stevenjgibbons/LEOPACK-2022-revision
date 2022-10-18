C*********************************************************************
C subroutine Magnetic and Convection Diffusion Right Hand Form *******
C            -            -          -         -     -    -    *******
C Steve Gibbons Mon Feb  7 19:05:42 GMT 2000                         C
C____________________________________________________________________C
C                                                                    C
C If we are time stepping the induction equation, this routine       C
C adds the highlighted terms to the right hand sides of              C
C                                                                    C
C ( CA + dt.CD.(c-1) Lap ) Theta^{i+1} =                             C
C                            ( CA + dt.CD.c Lap ) Theta^{i} +        C
C                            ^^^^^^^^^^^^^^^^^^^^      Forc. Term    C
C                                                                    C
C ( CE + dt.CI.(c-1) Lap ) omega^{i+1} =                             C
C                            ( CE + dt.CI.c Lap ) omega^{i} +        C
C                            ^^^^^^^^^^^^^^^^^^^^      Forc. Term    C
C                                                                    C
C  (omega is the vorticity = curl u)                                 C
C                                                                    C
C ( CK + dt.CL.(c-1) Lap ) B^{i+1} =                                 C
C                            ( CK + dt.CL.c Lap ) B^{i} +            C
C                            ^^^^^^^^^^^^^^^^^^^^      Forc. Term    C
C                                                                    C
C It does NOT zero the array RHS.                                    C
C                                                                    C
C It zeros the matrix, assumes that all boundary conditions are      C
C insulating, velocity boundary conditions are rigid,                C
C                                                                    C
C For a harmonic ih, MHT( ih ) must be 4 for poloidal and 5 for      C
C toroidal. MHL( ih ) gives the spherical harmonic degree, l.        C
C MHM( ih ) gives the spherical harmonic order, m, when ih has a     C
C (cos m phi) dependence and -m for a (sin m phi) dependence.        C
C                                                                    C
C MHP( ih ) = is where MHBCS( is ) = 2 when MHT( ih ) = 2 or 5,      C
C MHBCS( is ) = 7 when MHT( ih ) = 4 and MHBCS( is ) = 4 when        C
C MHT( ih ) = 1.                                                     C
C                                                                    C
C The magnetic field solution vector must be stored with IFORMF = 4. C
C (See indfun) i.e. INDFUN = ( IH - 1 )*NR + IR                      C
C                                                                    C
C____________________________________________________________________C
C                                                                    C
C Input variables :-                                                 C
C ===============                                                    C
C  Integer                                                           C
C  -------                                                           C
C                                                                    C
C     NR        : Number of radial grid nodes.                       C
C     NH        : Number of spherical harmonic radial functions.     C
C                                                                    C
C     MHT       : Dim ( NH ). See above.                             C
C     MHL       : Dim ( NH ). See above.                             C
C     MHM       : Dim ( NH ). See above.                             C
C                                                                    C
C     MHTR      : Dim ( NH ). MHT after calling CINDSW.              C
C                                                                    C
C  Double Precision                                                  C
C  ----------------                                                  C
C                                                                    C
C     CA        : Coef. of time deriv. in heat equation.             C
C     CE        : Coef. of time deriv. in vorticity equation.        C
C     CK        : Coef. of time deriv. in induction equation.        C
C                                                                    C
C     CD        : Coef. of diffusion in heat equation.               C
C     CI        : Coef. of diffusion in vorticity equation.          C
C     CL        : Coef. of diffusion in induction equation.          C
C                                                                    C
C     CFAC      : Determines how explicit/implicit time stepping is. C
C                 (Corresponds to 'c' in above equation).            C
C                                                                    C
C     DELTAT    : Time step.                                         C
C                                                                    C
C     V0        : Vector at time step i, zero^th derivatives.        C
C     V1        : Vector at time step i, first   derivatives.        C
C     V2        : Vector at time step i, second  derivatives.        C
C     V3        : Vector at time step i, third   derivatives.        C
C     V4        : Vector at time step i, fourth  derivatives.        C
C                                                                    C
C  (Pre-calculate these with CASVDR)                                 C
C                                                                    C
C     RHS       : Right hand side vector.                            C
C                                                                    C
C     XARR      : Dim ( NR ). i^{th} element gives x_i.              C
C                                                                    C
C____________________________________________________________________C
C
C*********************************************************************
      SUBROUTINE MCDRHF( NR, NH, MHT, MHL, MHM, MHTR, CA, CE, CK, CD,
     1                   CI, CL, CFAC, DELTAT, V0, V1, V2, V3, V4,
     2                   RHS, XARR )
      IMPLICIT NONE
C____________________________________________________________________C
C Variable declarations - Parameters ................................C
      INTEGER NR, NH, MHT( * ), MHL( * ), MHM( * ), MHTR( * )
      DOUBLE PRECISION CA, CE, CK, CD, CI, CL, CFAC, DELTAT,
     1                 V0( * ), V1( * ), V2( * ), V3( * ), V4( * ),
     2                 RHS( * ), XARR( * )
C____________________________________________________________________C
C Variable declarations - Working variables .........................C
      INTEGER INARR( 3 ), ILNR, IRNR, ICMP, ICMPO
      DOUBLE PRECISION FAC
C____________________________________________________________________C
C START OF PROGRAM **************************************************C
C____________________________________________________________________C
C
      INARR( 1 ) = 4
      INARR( 2 ) = NR
      INARR( 3 ) = NH
C     .
C     . Now we are ready to form the right hand side
C     .
C
      ILNR = 2
      IRNR = NR - 1
C
C Poloidal and toroidal velocity
C
      FAC = CI*DELTAT*CFAC
      CALL SSVLC( NR, V0, V1, V2, V3, V4, INARR, MHT,
     1            MHL, MHM, RHS, INARR, MHTR, MHL, MHM,
     2            FAC, ILNR, IRNR, XARR )
C
      FAC    = CE
      ICMP   = 1
      ICMPO  = 2
      CALL SSVCL( NR, V0, V1, V2, INARR, MHT, MHL, MHM,
     1            ICMP, RHS, INARR, MHTR, MHL, MHM, ICMPO,
     2            FAC, ILNR, IRNR, XARR )
C
      ICMP   = 2
      ICMPO  = 1
      CALL SSVCL( NR, V0, V1, V2, INARR, MHT, MHL, MHM,
     1            ICMP, RHS, INARR, MHTR, MHL, MHM, ICMPO,
     2            FAC, ILNR, IRNR, XARR )
C
C Temperature
C
      ICMP = 3
      FAC  = CD*DELTAT*CFAC
      CALL SSVLP( NR, V0, V1, V2, INARR, MHT, MHL, MHM,
     1            ICMP, RHS, INARR, MHTR, MHL, MHM, FAC,
     2            ILNR, IRNR, XARR )
C
      FAC  = CA
      CALL SSVTA( NR, V0, INARR, MHT, MHL, MHM, ICMP,
     1            RHS, INARR, MHTR, MHL, MHM, ICMP,
     2            FAC, ILNR, IRNR )
C
C Poloidal magnetic field
C
      ICMP = 4
      FAC  = CL*DELTAT*CFAC
      CALL SSVLP( NR, V0, V1, V2, INARR, MHT, MHL, MHM,
     1            ICMP, RHS, INARR, MHTR, MHL, MHM, FAC,
     2            ILNR, IRNR, XARR )
C
      FAC  = CK
      CALL SSVTA( NR, V0, INARR, MHT, MHL, MHM, ICMP,
     1            RHS, INARR, MHTR, MHL, MHM, ICMP,
     2            FAC, ILNR, IRNR )
C
C Toroidal magnetic field
C
      ICMP = 5
      FAC  = CL*DELTAT*CFAC
      CALL SSVLP( NR, V0, V1, V2, INARR, MHT, MHL, MHM,
     1            ICMP, RHS, INARR, MHTR, MHL, MHM, FAC,
     2            ILNR, IRNR, XARR )
C
      FAC  = CK
      CALL SSVTA( NR, V0, INARR, MHT, MHL, MHM, ICMP,
     1            RHS, INARR, MHTR, MHL, MHM, ICMP,
     2            FAC, ILNR, IRNR )
C
      RETURN
      END
C*********************************************************************
