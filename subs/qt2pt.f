C********************************************************************
C subroutine QT2PT **************************************************
C Steve Gibbons 12.4.97                                             C
C___________________________________________________________________C
C Takes the QST array and for a given harmonic NOHARM, calculates   C
C the Poloidal and Toroidal Coefficients for that harmonic          C
C according to equations 45 and 42b in my notes.                    C
C___________________________________________________________________C
C Input Variables :-                                                C
C ===============                                                   C
C  Integers                                                         C
C  --------                                                         C
C     LHMAX     : Max. number of harmonics (for dimensionning array.C
C     LEVH      : Level of wanted harmonic. If this variable is set C
C                  to -1 on input then LEVH is calculated by the    C
C                  routine LMFIND which calculates the Level, order C
C                  and type (i.e. cos or sin) from the number of    C
C                  the harmonic.                                    C
C     NOHARM    : Number of harmonic to be applied.                 C
C  Double Precision                                                 C
C  ----------------                                                 C
C     QST       : Array containing the QST coefficients as          C
C                  described in the PT2QST notes.                   C
C     RAD       : Radius at particular node.                        C
C___________________________________________________________________C
C Output Variables :-                                               C
C ================                                                  C
C  Double Precision                                                 C
C  ----------------                                                 C
C     POLC      : Poloidal Coefficient                              C
C     TORC      : Toroidal Coefficient                              C
C___________________________________________________________________C
C Working Variables :-                                              C
C -----------------                                                 C
C   Integer                                                         C
C   -------                                                         C
C      M        : Only used by LMFIND routine                       C
C      ITYPE    : ditto                                             C
C   Double Precision                                                C
C   ----------------                                                C
C      RLL1     : DBLE  ( L*(L+1) )                                 C
C      SQRLL1   : SQRT (RLL1)                                       C
C___________________________________________________________________C
C Subroutines used :-                                               C
C ----------------                                                  C
C      LMFIND   : Only if levh=-1                                   C
C___________________________________________________________________C
C
C********************************************************************
      SUBROUTINE QT2PT ( LHMAX, LEVH, NOHARM, QST, RAD, POLC, TORC)
      IMPLICIT NONE
C___________________________________________________________________C
C Variable declarations - Parameters ...............................C
      INTEGER LHMAX, LEVH, NOHARM
      DOUBLE PRECISION RAD, POLC, TORC,
     1                 QST ( LHMAX*(LHMAX + 2) , 3)
C___________________________________________________________________C
C Variable declarations - Working Variables ........................C
      INTEGER M, ITYPE
      DOUBLE PRECISION RLL1,SQRLL1
C___________________________________________________________________C
C START OF PROGRAM *************************************************C
C___________________________________________________________________C
      IF ( LEVH.EQ.-1 ) THEN
         CALL LMFIND( NOHARM, LEVH, M, ITYPE)
      ENDIF    

      RLL1 = DBLE ( LEVH*LEVH + LEVH )
      SQRLL1 = SQRT ( RLL1 )
      POLC = QST( NOHARM, 1) * RAD / RLL1
      TORC = -QST( NOHARM, 3)/SQRLL1

      RETURN
      END
C********************************************************************
