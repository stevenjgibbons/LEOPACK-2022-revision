C*********************************************************************
C subroutine Old Format Vector File Read *****************************
C            -   -      -      -    -    *****************************
C Steve Gibbons Wed Sep 20 14:51:13 BST 2000                         C
C____________________________________________________________________C
C                                                                    C
C Reads in an old Leeds-style .evecs file and converts to the new    C
C format vector. A call to this routine MUST have been preceded by   C
C a call to OFTFR, to obtain the appropriate indices. We must also   C
C have supplied values for the inner and outer radii (as these are   C
C not necessarily stored in the .evecs file). The number of radial   C
C grid nodes are checked against a maximum number. The number of     C
C radial functions must agree exactly with the value given by NH.    C
C The elements of VEC are all divided by the radius, due to the      C
C different definition of the old (Bullard and Gellman) and the      C
C new (Zhang and Busse) radial functions.                            C
C____________________________________________________________________C
C                                                                    C
C Input variables :-                                                 C
C ===============                                                    C
C  Integer                                                           C
C  -------                                                           C
C                                                                    C
C     NH        : Number of vector spherical harmonics.              C
C     NRMAX     : Maximum number of radial grid nodes.               C
C     INARR     : See INDFUN for details.                            C
C                  INARR( 1 ) is automatically set to 4.             C
C                  INARR( 2 ) is set to NR, given in .evecs file     C
C                  INARR( 3 ) = NH.                                  C
C     LH        : Maximum spherical harmonic degree (in file)        C
C     LU        : Logical file unit number.                          C
C                                                                    C
C  Character                                                         C
C  ---------                                                         C
C                                                                    C
C     FNAME     : *(*) File name.                                    C
C                                                                    C
C  Double Precision                                                  C
C  ----------------                                                  C
C                                                                    C
C     VEC       : Vector read from file. Dim ( * )                   C
C     RI        : Radius of inner boundary.                          C
C     RO        : Radius of outer boundary.                          C
C                                                                    C
C____________________________________________________________________C
C
C*********************************************************************
      SUBROUTINE OFVFR( NH, NRMAX, INARR, LH, LU, FNAME, VEC, RI, RO)
      IMPLICIT NONE
C____________________________________________________________________C
C Variable declarations - Parameters ................................C
      INTEGER NH, NRMAX, INARR( * ), LH, LU
      CHARACTER *(*) FNAME
      DOUBLE PRECISION VEC( * ), RI, RO
C____________________________________________________________________C
C Variable declarations - Working variables .........................C
      INTEGER IH, IR, NR, IWR, I, NELTS, IND, INDFUN
      DOUBLE PRECISION DLOW, RAD, H
      CHARACTER *(1) CHJUNK
      PARAMETER ( IWR = 1, DLOW = 1.0d-6 )
C____________________________________________________________________C
C START OF PROGRAM **************************************************C
C____________________________________________________________________C
C
C Open file for reading
C
      CALL FOPEN ( LU, FNAME, IWR )
C
C First and second lines contains junk
C
      READ ( LU, 99 ) CHJUNK
      READ ( LU, 99 ) CHJUNK
C
      INARR( 1 ) = 4
      INARR( 3 ) = NH
C
C Third line contains LH, NH, NR
C
      READ ( LU, * )  LH, NH, NR
C
      IF ( NH.NE.INARR( 3 ) ) THEN
        PRINT *,' Subroutine OFVFR.'
        PRINT *,' NH = ',NH,' INARR(3) = ', INARR( 3 )
        PRINT *,' Program aborted.'
        STOP
      ENDIF
C
      IF ( NR.GT.NRMAX ) THEN
        PRINT *,' Subroutine OFVFR.'
        PRINT *,' NR = ',NR,' NRMAX = ', NRMAX
        PRINT *,' Program aborted.'
        STOP
      ENDIF
C
      INARR( 2 ) = NR
      NELTS = NR*NH
      READ ( LU, * ) ( VEC( I ), I = 1, NELTS )
C
C Close file
C
      CALL FCLOSE ( LU, FNAME, 'Error closing file.' )
C
C Now we need to rescale the vector values.
C
      IF ( NR.LE.2 .OR. (RO-RI).LT.DLOW ) THEN
        PRINT *,' Subroutine OFVFR.'
        PRINT *,' NR = ',NR,' RI = ', RI,' RO = ', RO
        PRINT *,' Program aborted.'
        STOP
      ENDIF
      H = (RO-RI)/DBLE( NR - 1 )
C
      DO IH = 1, NH
        DO IR = 2, NR
          IND = INDFUN( IR, IH, INARR )
          RAD = RI + DBLE( IR - 1 )*H
          VEC( IND ) = VEC( IND )/RAD
        ENDDO
      ENDDO
C
 99   FORMAT(A)
      RETURN
      END
C*********************************************************************
