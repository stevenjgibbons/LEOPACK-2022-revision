C*********************************************************************
C subroutine Vector Function Single Component EXtract ****************
C            -      -        -      -         --      ****************
C Steve Gibbons Thu Feb 24 18:01:58 GMT 2000                         C
C____________________________________________________________________C
C                                                                    C
C Given a vector function VF( NPHP, NTHP, 3 ) and a component        C
C number ICMP, VFSCEX will put the                                   C
C                                                                    C
C  SF( IPHI, ITHE ) = VF( IPHI, ITHE, ICMP )                         C
C                        for all values IPHI and ITHE.               C
C                                                                    C
C____________________________________________________________________C
C                                                                    C
C Input variables :-                                                 C
C ===============                                                    C
C  Integer                                                           C
C  -------                                                           C
C                                                                    C
C     NPHP      : Number of points in phi direction.                 C
C     NTHP      : Number of points in theta direction.               C
C     ICMP      : Component to be selected.                          C
C                                                                    C
C  Double Precision                                                  C
C  ----------------                                                  C
C                                                                    C
C     VF        : Vector function dim ( NPHP, NTHP, 3 )              C
C     SF        : Scalar function dim ( NPHP, NTHP )                 C
C____________________________________________________________________C
C
C*********************************************************************
      SUBROUTINE VFSCEX( NPHP, NTHP, ICMP, VF, SF )
      IMPLICIT NONE
C____________________________________________________________________C
C Variable declarations - Parameters ................................C
      INTEGER NPHP, NTHP, ICMP
      DOUBLE PRECISION VF( NPHP, NTHP, 3 ), SF( NPHP, NTHP )
C____________________________________________________________________C
C Variable declarations - Working variables .........................C
      INTEGER IPHI, ITHE
C____________________________________________________________________C
C START OF PROGRAM **************************************************C
C____________________________________________________________________C
C
C Check the validity of arguments .....
C
      IF ( ICMP.NE.1 .AND. ICMP.NE.2 .AND. ICMP.NE.3 ) THEN
        PRINT *,' Subroutine VFSCEX.'
        PRINT *,' ICMP = ', ICMP
        PRINT *,' Program aborted.'
        STOP
      ENDIF
C
      DO ITHE = 1, NTHP
        DO IPHI = 1, NPHP
          SF( IPHI, ITHE ) = VF( IPHI, ITHE, ICMP )
        ENDDO
      ENDDO
C
      RETURN
      END
C*********************************************************************
