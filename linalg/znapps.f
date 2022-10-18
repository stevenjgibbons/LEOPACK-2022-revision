c\BeginDoc
c
c\Name: znapps
c
c\Description:
c  Given the Arnoldi factorization
c
c     A*V_{k} - V_{k}*H_{k} = r_{k+p}*e_{k+p}^T,
c
c  apply NP implicit shifts resulting in
c
c     A*(V_{k}*Q) - (V_{k}*Q)*(Q^T* H_{k}*Q) = r_{k+p}*e_{k+p}^T * Q
c
c  where Q is an orthogonal matrix which is the product of rotations
c  and reflections resulting from the NP bulge change sweeps.
c  The updated Arnoldi factorization becomes:
c
c     A*VNEW_{k} - VNEW_{k}*HNEW_{k} = rnew_{k}*e_{k}^T.
c
c\Usage:
c  call znapps
c     ( N, KEV, NP, SHIFT, V, LDV, H, LDH, RESID, Q, LDQ, 
c       WORKL, WORKD )
c
c\Arguments
c  N       Integer.  (INPUT)
c          Problem size, i.e. size of matrix A.
c
c  KEV     Integer.  (INPUT/OUTPUT)
c          KEV+NP is the size of the input matrix H.
c          KEV is the size of the updated matrix HNEW.  KEV is only 
c          updated on ouput when fewer than NP shifts are applied in
c          order to keep the conjugate pair together.
c
c  NP      Integer.  (INPUT)
c          Number of implicit shifts to be applied.
c
c  SHIFT   Complex*16 array of length NP.  (INPUT)
c          The shifts to be applied.
c
c  V       Complex*16 N by (KEV+NP) array.  (INPUT/OUTPUT)
c          On INPUT, V contains the current KEV+NP Arnoldi vectors.
c          On OUTPUT, V contains the updated KEV Arnoldi vectors
c          in the first KEV columns of V.
c
c  LDV     Integer.  (INPUT)
c          Leading dimension of V exactly as declared in the calling
c          program.
c
c  H       Complex*16 (KEV+NP) by (KEV+NP) array.  (INPUT/OUTPUT)
c          On INPUT, H contains the current KEV+NP by KEV+NP upper 
c          Hessenberg matrix of the Arnoldi factorization.
c          On OUTPUT, H contains the updated KEV by KEV upper Hessenberg
c          matrix in the KEV leading submatrix.
c
c  LDH     Integer.  (INPUT)
c          Leading dimension of H exactly as declared in the calling
c          program.
c
c  RESID   Complex*16 array of length N.  (INPUT/OUTPUT)
c          On INPUT, RESID contains the the residual vector r_{k+p}.
c          On OUTPUT, RESID is the update residual vector rnew_{k} 
c          in the first KEV locations.
c
c  Q       Complex*16 KEV+NP by KEV+NP work array.  (WORKSPACE)
c          Work array used to accumulate the rotations and reflections
c          during the bulge chase sweep.
c
c  LDQ     Integer.  (INPUT)
c          Leading dimension of Q exactly as declared in the calling
c          program.
c
c  WORKL   Complex*16 work array of length (KEV+NP).  (WORKSPACE)
c          Private (replicated) array on each PE or array allocated on
c          the front end.
c
c  WORKD   Complex*16 work array of length 2*N.  (WORKSPACE)
c          Distributed array used in the application of the accumulated
c          orthogonal matrix Q.
c
c\EndDoc
c
c-----------------------------------------------------------------------
c
c\BeginLib
c
c\Local variables:
c     xxxxxx  Complex*16
c
c\References:
c  1. D.C. Sorensen, "Implicit Application of Polynomial Filters in
c     a k-Step Arnoldi Method", SIAM J. Matr. Anal. Apps., 13 (1992),
c     pp 357-385.
c
c\Routines called:
c     zlartg  LAPACK Givens rotation construction routine.
c     zlacpy  LAPACK matrix copy routine.
c     zlaset  LAPACK matrix initialization routine.
c     dlabad  LAPACK routine for defining the underflow and overflow
c             limits.
c     zaxpy   Level 1 BLAS that computes a vector triad.
c     zscal   Level 1 BLAS that scales a vector.
c     zcopy   Level 1 BLAS that copies one vector to another.
c
c\Author
c     Danny Sorensen               Phuong Vu
c     Richard Lehoucq              CRPC / Rice University
c     Dept. of Computational &     Houston, Texas
c     Applied Mathematics 
c     Rice University           
c     Houston, Texas 
c
c\SCCS Information: @(#)
c FILE: napps.F   SID: 2.1   DATE OF SID: 11/21/95   RELEASE: 2
c
c\Remarks
c  1. In this version, each shift is applied to all the sublocks of
c     the Hessenberg matrix H and not just to the submatrix that it
c     comes from. Deflation as in LAPACK routine clahqr (QR algorithm
c     for upper Hessenberg matrices ) is used.
c     The subdiagonals of H are enforced to be non-negative.
c
c\EndLib
c
c-----------------------------------------------------------------------
c
      subroutine znapps
     &   ( n, kev, np, shift, v, ldv, h, ldh, resid, q, ldq, 
     &     workl, workd )
c
c     %----------------------------------------------------%
c     | Include files for debugging and timing information |
c     %----------------------------------------------------%
c
      include   'debug.h'
      include   'stat.h'
c
c     %------------------%
c     | Scalar Arguments |
c     %------------------%
c
      integer    kev, ldh, ldq, ldv, n, np
c
c     %-----------------%
c     | Array Arguments |
c     %-----------------%
c
      Complex*16
     &           h(ldh,kev+np), resid(n), shift(np), 
     &           v(ldv,kev+np), q(ldq,kev+np), workd(2*n), workl(kev+np)
c
c     %------------%
c     | Parameters |
c     %------------%
c
      Complex*16
     &           one, zero
      parameter (one = (1.0, 0.0), zero = (0.0, 0.0))
c
c     %------------------------%
c     | Local Scalars & Arrays |
c     %------------------------%
c
      integer    i, iend, istart, j, jj, kplusp, msglvl
      logical    first
      Complex*16
     &           f, g, h11, h21, r, s, sigma, t
      Double precision             
     &           c,  ovfl, smlnum, ulp, unfl, tst1
      save       first, ovfl, smlnum, ulp, unfl 
c
c     %----------------------%
c     | External Subroutines |
c     %----------------------%
c
      external   zaxpy, zcopy, zgemv, zscal, zlacpy, zlartg, 
     &           zvout, zlaset, dlabad, zmout, second, ivout
c
c     %--------------------%
c     | External Functions |
c     %--------------------%
c
      Double precision                 
     &           dlamch, zlanhs
      external   dlamch, zlanhs
c
c     %----------------------%
c     | Intrinsics Functions |
c     %----------------------%
c
      intrinsic  abs, max, min
c
c     %----------------%
c     | Data statments |
c     %----------------%
c
      data       first / .true. /
c
c     %-----------------------%
c     | Executable Statements |
c     %-----------------------%
c
      if (first) then
c
c        %-----------------------------------------------%
c        | Set machine-dependent constants for the       |
c        | stopping criterion. If norm(H) <= sqrt(OVFL), |
c        | overflow should not occur.                    |
c        | REFERENCE: LAPACK subroutine clahqr           |
c        %-----------------------------------------------%
c
         unfl = dlamch( 'safe minimum' )
         ovfl = dble(one / unfl)
         call dlabad( unfl, ovfl )
         ulp = dlamch( 'precision' )
         smlnum = unfl*( n / ulp )
         first = .false.
      end if
c
c     %-------------------------------%
c     | Initialize timing statistics  |
c     | & message level for debugging |
c     %-------------------------------%
c
      call second (t0)
      msglvl = mcapps
c 
      kplusp = kev + np 
c 
c     %--------------------------------------------%
c     | Initialize Q to the identity to accumulate |
c     | the rotations and reflections              |
c     %--------------------------------------------%
c
      call zlaset ('All', kplusp, kplusp, zero, one, q, ldq)
c
c     %----------------------------------------------%
c     | Chase the bulge with the application of each |
c     | implicit shift. Each shift is applied to the |
c     | whole matrix including each block.           |
c     %----------------------------------------------%
c
      do 110 jj = 1, np
         sigma = shift(jj)
c
         if (msglvl .gt. 2 ) then
            call ivout (logfil, 1, jj, ndigit, 
     &               '_napps: shift number.')
            call zvout (logfil, 1, sigma, ndigit, 
     &               '_napps: Value of the shift ')
         end if
c
         istart = 1
   20    continue
c
         do 30 i = istart, kplusp-1
c
c           %----------------------------------------%
c           | Check for splitting and deflation. Use |
c           | a standard test as in the QR algorithm |
c           | REFERENCE: LAPACK subroutine clahqr    |
c           %----------------------------------------%
c
            tst1 = abs( h( i, i ) ) + abs( h( i+1, i+1 ) )
            if( tst1.eq.zero )
     &         tst1 = zlanhs( '1', kplusp-jj+1, h, ldh, workl )
            if ( abs( h( i+1,i ) ).le.max( ulp*tst1, smlnum ))  then
               if (msglvl .gt. 0) then
                  call ivout (logfil, 1, i, ndigit, 
     &                 '_napps: matrix splitting at row/column no.')
                  call ivout (logfil, 1, jj, ndigit, 
     &                 '_napps: matrix splitting with shift number.')
                  call zvout (logfil, 1, h(i+1,i), ndigit, 
     &                 '_napps: off diagonal element.')
               end if
               iend = i
               h(i+1,i) = zero
               go to 40
            end if
   30    continue
         iend = kplusp
   40    continue
c
         if (msglvl .gt. 2) then
             call ivout (logfil, 1, istart, ndigit, 
     &                   '_napps: Start of current block ')
             call ivout (logfil, 1, iend, ndigit, 
     &                   '_napps: End of current block ')
         end if
c
c        %------------------------------------------------%
c        | No reason to apply a shift to block of order 1 |
c        | or if the current block starts after the point |
c        | of compression since we'll discard this stuff  |
c        %------------------------------------------------%
c
         if ( istart .eq. iend .or. istart .gt. kev) go to 100
c
         if (msglvl .gt. 3) then
            call zmout (logfil, kplusp, kplusp, h, ldh, ndigit,
     &           '_napps: Hessenberg matrix H before the shift')
         end if
c
         h11 = h(istart,istart)
         h21 = h(istart+1,istart)
            f = h11 - sigma
            g = h21
c 
            do 80 i = istart, iend-1
c
c              %------------------------------------------------------%
c              | Construct the plane rotation G to zero out the bulge |
c              %------------------------------------------------------%
c
               call zlartg (f, g, c, s, r)
               if (i .gt. istart) then
c
c                 %-------------------------------------------%
c                 | The following ensures that h(1:iend-1,1), |
c                 | the first iend-2 off diagonal of elements |
c                 | H, remain non negative.                   |
c                 %-------------------------------------------%
c
c                 if (r .lt. zero) then
c                    r = -r
c                    c = -c
c                    s = -s
c                 end if
                  h(i,i-1) = r
                  h(i+1,i-1) = zero
               end if
c
c              %--------------------------------------------%
c              | Apply rotation to the left of H;  H <- G*H |
c              %--------------------------------------------%
c
               do 50 j = i, kplusp
                  t        =  c*h(i,j) + s*h(i+1,j)
                  h(i+1,j) = -conjg(s)*h(i,j) + c*h(i+1,j)
                  h(i,j)   = t   
   50          continue
c
c              %----------------------------------------------%
c              | Apply rotation to the right of H;  H <- H*G' |
c              %----------------------------------------------%
c
               do 60 j = 1, min(i+2,iend)
                  t        =  c*h(j,i) + conjg(s)*h(j,i+1)
                  h(j,i+1) = -s*h(j,i) + c*h(j,i+1)
                  h(j,i)   = t   
   60          continue
c
c              %-----------------------------------------------------%
c              | Accumulate the rotation in the matrix Q;  Q <- Q*G' |
c              %-----------------------------------------------------%
c
               do 70 j = 1, kplusp
                  t        =   c*q(j,i) + conjg(s)*q(j,i+1)
                  q(j,i+1) = - s*q(j,i) + c*q(j,i+1)
                  q(j,i)   = t   
   70          continue
c
c              %---------------------------%
c              | Prepare for next rotation |
c              %---------------------------%
c
               if (i .lt. iend-1) then
                  f = h(i+1,i)
                  g = h(i+2,i)
               end if
   80       continue
c
c           %-------------------------------%
c           | Finished applying the shift.  |
c           %-------------------------------%
c 
  100    continue
c
c        %---------------------------------------------------------%
c        | Apply the same shift to the next block if there is any. |
c        %---------------------------------------------------------%
c
         istart = iend + 1
         if (iend .lt. kplusp) go to 20
c
         if (msglvl .gt. 3) then
            call zmout (logfil, kplusp, kplusp, h, ldh, ndigit,
     &           '_napps: Hessenberg matrix H after this shift')
         end if
c
c        %---------------------------------------------%
c        | Loop back to the top to get the next shift. |
c        %---------------------------------------------%
c
  110 continue
c
c     %--------------------------------------------------%
c     | Perform a similarity transformation that makes   |
c     | sure that H will have non-negative subdiagonals. |
c     %--------------------------------------------------%
c
c     do 120 j=1,kplusp-1
c        if ( h(j+1,j) .lt. zero ) then
c             call zscal( kplusp-j+1, -one, h(j+1,j), ldh )
c             call zscal( min(j+2, kplusp), -one, h(1,j+1), 1 )
c             call zscal( kplusp, -one, q(1,j+1), 1 )
c        end if
c120  continue
c
      do 130 i = 1, kev-1
c
c        %--------------------------------------------%
c        | Final check for splitting and deflation.   |
c        | Use a standard test as in the QR algorithm |
c        | REFERENCE: LAPACK subroutine clahqr        |
c        %--------------------------------------------%
c
c        tst1 = abs( h( i, i ) ) + abs( h( i+1, i+1 ) )
         tst1 = .00001*(abs( h( i, i ) ) + abs( h( i+1, i+1 ) ))
         if( tst1.eq.zero )
     &       tst1 = zlanhs( '1', kev, h, ldh, workl )
         if( abs( h( i+1,i ) ).le.max( ulp*tst1, smlnum ) ) 
     &       h(i+1,i) = zero
 130  continue
c
      if (msglvl .gt. 4) then
         call zmout (logfil, kplusp, kplusp, q, ldq, ndigit,
     &              '_napps: accumulated Hessenberg matrix Q')
      end if
c
c     %-------------------------------------------------%
c     | Compute the (kev+1)-st column of (V*Q) and      |
c     | temporarily store the result in WORKD(N+1:2*N). |
c     | This is needed in the residual update since we  |
c     | cannot GUARANTEE that the corresponding entry   |
c     | of H would be zero as in exact arithmetic.      |
c     %-------------------------------------------------%
c
      call zgemv ('N', n, kplusp, one, v, ldv,
     &            q(1,kev+1), 1, zero, workd(n+1), 1)
c 
c     %----------------------------------------------------------%
c     | Compute column 1 to kev of (V*Q) in backward order       |
c     | taking advantage of the upper Hessenberg structure of Q. |
c     %----------------------------------------------------------%
c
      do 140 i = 1, kev
         call zgemv ('N', n, kplusp-i+1, one, v, ldv,
     &               q(1,kev-i+1), 1, zero, workd, 1)
         call zcopy (n, workd, 1, v(1,kplusp-i+1), 1)
  140 continue
c
c     %-------------------------------------------------%
c     |  Move v(:,kplusp-kev+1:kplusp) into v(:,1:kev). |
c     %-------------------------------------------------%
c
      call zlacpy ('A', n, kev, v(1,kplusp-kev+1), ldv, v, ldv)
c 
c     %--------------------------------------------------------------%
c     | Copy the (kev+1)-st column of (V*Q) in the appropriate place |
c     %--------------------------------------------------------------%
c
      call zcopy (n, workd(n+1), 1, v(1,kev+1), 1)
c 
c     %-------------------------------------%
c     | Update the residual vector:         |
c     |    r <- sigmak*r + betak*v(:,kev+1) |
c     | where                               |
c     |    sigmak = (e_{kev+p}'*Q)*e_{kev}  |
c     |    betak = e_{kev+1}'*H*e_{kev}     |
c     %-------------------------------------%
c
      call zscal (n, q(kplusp,kev), resid, 1)
      call zaxpy (n, h(kev+1,kev), v(1,kev+1), 1, resid, 1)
c
      if (msglvl .gt. 1) then
         call zvout (logfil, 1, q(kplusp,kev), ndigit,
     &        '_napps: sigmak = (e_{kev+p}^T*Q)*e_{kev}')
         call zvout (logfil, 1, h(kev+1,kev), ndigit,
     &        '_napps: betak = e_{kev+1}^T*H*e_{kev}')
         call ivout (logfil, 1, kev, ndigit, 
     &               '_napps: Order of the final Hessenberg matrix ')
         call zmout (logfil, kev, kev, h, ldh, ndigit,
     &        '_napps: updated Hessenberg matrix H for next iteration')
      end if
      if (msglvl .gt. 3) then
         call zvout (logfil, n, resid, ndigit,
     &        '_napps: updated residual for next iteration')
      end if
      if (msglvl .gt. 4) then
         call zmout (logfil, n, kev+1, v, ldv, ndigit,
     &        '_napps: updated matrix of Arnoldi vectors V')
      end if
c 
      call second (t1)
      tcapps = tcapps + (t1 - t0)
c 
      return
c
c     End of znapps
c
      end
