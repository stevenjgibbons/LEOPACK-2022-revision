C*********************************************************************
C function KLM(l,m)=l*(l+1)/2+1+m ... I'm sorry but not even I will  C
C write any more comments for this function!! Steve Gibbons 16.4.97  C
C____________________________________________________________________C
      FUNCTION KLM ( L, M )
      IMPLICIT NONE
      INTEGER KLM,L,M
      KLM = L*(L+1)/2+M+1
      RETURN
      END
C*********************************************************************
