c
c\SCCS Information: @(#)
c FILE: statn.F   SID: 2.1   DATE OF SID: 11/21/95   RELEASE: 2
c
c     %---------------------------------------------%
c     | Initialize statistic and timing information |
c     | for complex nonsymmetric Arnoldi code.      |
c     %---------------------------------------------%

      subroutine zstatn
c
c     %--------------------------------%
c     | See stat.doc for documentation |
c     %--------------------------------%
c
      include   'stat.h'
 
c     %-----------------------%
c     | Executable Statements |
c     %-----------------------%

      nopx   = 0
      nbx    = 0
      nrorth = 0
      nitref = 0
      nrstrt = 0
 
      tcaupd = 0.0
      tcaup2 = 0.0
      tcaitr = 0.0
      tceigh = 0.0
      tcgets = 0.0
      tcapps = 0.0
      tcconv = 0.0
      titref = 0.0
      tgetv0 = 0.0
      trvec  = 0.0
 
c     %----------------------------------------------------%
c     | User time including reverse communication overhead |
c     %----------------------------------------------------%
      tmvopx = 0.0
      tmvbx  = 0.0
 
      return
c
c     End of zstatn
c
      end
