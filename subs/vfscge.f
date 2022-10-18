C*********************************************************************
C subroutine Vector Function Single Component Global Extract *********
C            -      -        -      -         -      -       *********
C Steve Gibbons Thu Mar  2 10:52:37 GMT 2000                         C
C____________________________________________________________________C
C                                                                    C
C Given a vector function VF( NPHP, NTHP, 3 ) and a component        C
C number ICMP, VFSCGE will put the                                   C
C                                                                    C
C  SF( ITHE, IPHI ) = VF( IPHI, ITHE, ICMP )                         C
C                        for all values IPHI and ITHE.               C
C  and also SF( ITHE, NPHP + 1) = SF( ITHE, 1 )                      C
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
C     SF        : Scalar function dim ( NTHP, NPHP + 1 )             C
C____________________________________________________________________C
C
C*********************************************************************
      SUBROUTINE VFSCGE( NPHP, NTHP, ICMP, VF, SF )
      IMPLICIT NONE
C____________________________________________________________________C
C Variable declarations - Parameters ................................C
      INTEGER NPHP, NTHP, ICMP
      DOUBLE PRECISION VF( NPHP, NTHP, 3 ), SF( NTHP, NPHP+1 )
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
        PRINT *,' Subroutine VFSCGE.'
        PRINT *,' ICMP = ', ICMP
        PRINT *,' Program aborted.'
        STOP
      ENDIF
C
      DO ITHE = 1, NTHP
        DO IPHI = 1, NPHP
          SF( ITHE, IPHI ) = VF( IPHI, ITHE, ICMP )
        ENDDO
        SF( ITHE, NPHP + 1 ) = SF( ITHE, 1 )
      ENDDO
C
      RETURN
      END
C*********************************************************************
