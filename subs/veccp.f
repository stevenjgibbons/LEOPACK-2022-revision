C*********************************************************************
C subroutine VECtor CoPy *********************************************
C            ---    - -  *********************************************
C Steve Gibbons 30.4.97                                              C
C____________________________________________________________________C
C Simply copies vector V1 into V2 ... where V1, V2 are double        C
C precision arrays of dimension ( NDIM ).                            C
C____________________________________________________________________C
      SUBROUTINE VECCP( V1, V2, NDIM )
      IMPLICIT NONE
      INTEGER NDIM, I
      DOUBLE PRECISION V1( NDIM ), V2( NDIM )

      DO I = 1, NDIM 
         V2( I ) = V1( I )
      ENDDO
      RETURN
      END
C*********************************************************************

