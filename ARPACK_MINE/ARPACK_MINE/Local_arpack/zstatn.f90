!
!\SCCS Information: @(#)
! FILE: statn.F   SID: 2.2   DATE OF SID: 4/20/96   RELEASE: 2
!
!     %---------------------------------------------%
!     | Initialize statistic and timing information |
!     | for complex nonsymmetric Arnoldi code.      |
!     %---------------------------------------------%

      SUBROUTINE zstatn
!
!     %--------------------------------%
!     | See stat.doc for documentation |
!     %--------------------------------%
!
      USE STAT_MODULE
      IMPLICIT NONE
 
!     %-----------------------%
!     | Executable Statements |
!     %-----------------------%

      nopx   = 0
      nbx    = 0
      nrorth = 0
      nitref = 0
      nrstrt = 0
 
      tcaupd = 0.0D+0
      tcaup2 = 0.0D+0
      tcaitr = 0.0D+0
      tceigh = 0.0D+0
      tcgets = 0.0D+0
      tcapps = 0.0D+0
      tcconv = 0.0D+0
      titref = 0.0D+0
      tgetv0 = 0.0D+0
      trvec  = 0.0D+0
 
!     %----------------------------------------------------%
!     | User time including reverse communication overhead |
!     %----------------------------------------------------%
      tmvopx = 0.0D+0
      tmvbx  = 0.0D+0
 
      return
!
!     %---------------%
!     | End of zstatn |
!     %---------------%
!
      end
