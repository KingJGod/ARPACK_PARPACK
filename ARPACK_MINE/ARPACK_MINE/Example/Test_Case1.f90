SUBROUTINE ZNSIMP_CASE1
!
!     This example program is intended to illustrate the 
!     simplest case of using ARPACK in considerable detail.  
!     This code may be used to understand basic usage of ARPACK
!     and as a template for creating an interface to ARPACK.  
!   
!     This code shows how to use ARPACK to find a few eigenvalues 
!     (lambda) and corresponding eigenvectors (x) for the standard 
!     eigenvalue problem:
!          
!                        A*x = lambda*x
! 
!     where A is a general n by n complex matrix.
!
!     The main points illustrated here are 
!
!        1) How to declare sufficient memory to find NEV 
!           eigenvalues of largest magnitude.  Other options
!           are available.
!
!        2) Illustration of the reverse communication interface 
!           needed to utilize the top level ARPACK routine ZNAUPD  
!           that computes the quantities needed to construct
!           the desired eigenvalues and eigenvectors(if requested).
!
!        3) How to extract the desired eigenvalues and eigenvectors
!           using the ARPACK routine ZNEUPD .
!
!     The only thing that must be supplied in order to use this
!     routine on your problem is to change the array dimensions 
!     appropriately, to specify WHICH eigenvalues you want to compute 
!     and to supply a matrix-vector product
!
!                         w <-  znsav
!
!     in place of the call to znsav( )  below.
!
!
!     Once usage of this routine is understood, you may wish to explore
!     the other available options to improve convergence, to solve generalized
!     problems, etc.  Look at the file ex-complex.doc in DOCUMENTS directory.
!     This codes implements  
!
!
!\Example-1
!     ... Suppose we want to solve A*x = lambda*x in regular mode,
!     ... OP = A  and  B = I.
!     ... Assume "call znsav (nx,x,y)" computes y = A*x
!     ... Use mode 1 of ZNAUPD .
!
!\BeginLib
!
!\Routines called
!     znaupd   ARPACK reverse communication interface routine.
!     zneupd   ARPACK routine that returns Ritz values and (optionally)
!             Ritz vectors.
!     dlapy2   LAPACK routine to compute sqrt(x**2+y**2) carefully.
!     dznrm2   Level 1 BLAS that computes the norm of a complex vector.
!     zaxpy    Level 1 BLAS that computes y <- alpha*x+y.
!     znsav      Matrix vector multiplication routine that computes A*x.
!     znstv      Matrix vector multiplication routine that computes T*x,
!             where T is a tridiagonal matrix.  It is used in routine
!             znsav.
!
!\Author
!     Richard Lehoucq
!     Danny Sorensen 
!     Chao Yang 
!     Dept. of Computational & 
!     Applied Mathematics 
!     Rice University 
!     Houston, Texas
!
!\SCCS Information: @(#)
! FILE: nsimp.F   SID: 2.4   DATE OF SID: 10/20/00   RELEASE: 2
!
!\Remarks
!     1. None
!
!\EndLib
!---------------------------------------------------------------------------
!
!     %------------------------------------------------------%
!     | Storage Declarations:                                |
!     |                                                      |
!     | The maximum dimensions for all arrays are            |
!     | set here to accommodate a problem size of            |
!     | N .le. MAXN                                          |
!     |                                                      |
!     | NEV is the number of eigenvalues requested.          |
!     |     See specifications for ARPACK usage below.       |
!     |                                                      |
!     | NCV is the largest number of basis vectors that will |
!     |     be used in the Implicitly Restarted Arnoldi      |
!     |     Process.  Work per major iteration is            |
!     |     proportional to N*NCV*NCV.                       |
!     |                                                      |
!     | You must set:                                        |
!     |                                                      |
!     | MAXN:   Maximum dimension of the A allowed.          |
!     | MAXNEV: Maximum NEV allowed.                         |
!     | MAXNCV: Maximum NCV allowed.                         |
!     %------------------------------------------------------%
!
  USE DEBUG_MODULE
  USE STAT_MODULE
  IMPLICIT NONE
  
  INTEGER( kind = 4 ), PARAMETER:: maxn=20000
  INTEGER( kind = 4 ), PARAMETER:: maxnev=120
  INTEGER( kind = 4 ), PARAMETER:: maxncv=200
  INTEGER( kind = 4 ), PARAMETER:: ldv=20000
!  
!     %--------------%
!     | Local Arrays |
!     %--------------%
!  
  INTEGER( kind = 4 ):: IPARAM(11),IPNTR(14)
  LOGICAL::SELECT(MAXNCV)
  COMPLEX( kind = 8 ):: ax(maxn), d(maxncv)
  COMPLEX( kind = 8 ):: v(ldv,maxncv), workd(3*maxn)
  COMPLEX( kind = 8 ):: workev(2*maxncv), resid(maxn)
  COMPLEX( kind = 8 ):: workl(3*maxncv*maxncv+5*maxncv)
  REAL( kind = 8 ):: rwork(maxncv), rd(maxncv,3)
  
!  
!     %---------------%
!     | Local Scalars |
!     %---------------%
!  
  CHARACTER(len=1)::bmat
  CHARACTER(len=2)::which
  INTEGER( kind = 4 ):: ido, n, nx, nev, ncv, lworkl, info, ierr
  INTEGER( kind = 4 ):: j, ishfts, maxitr, mode1, nconv
  COMPLEX( kind = 8 ):: sigma
  REAL( kind = 8 )::    tol
  LOGICAL::             rvec
!
!     %-----------------------------%
!     | BLAS & LAPACK routines used |
!     %-----------------------------%
!  
  REAL( kind = 8 )::dznrm2 , dlapy2 
  EXTERNAL dznrm2 , zaxpy , dlapy2
!
!     %-----------------------%
!     | Executable Statements |
!     %-----------------------%
!
!     %-------------------------------------------------%
!     | The following include statement and assignments |
!     | initiate trace output from the internal         |
!     | actions of ARPACK.  See debug.doc in the        |
!     | DOCUMENTS directory for usage.  Initially, the  |
!     | most useful information will be a breakdown of  |
!     | time spent in the various stages of computation |
!     | given by setting mcaupd = 1                     |
!     %-------------------------------------------------%
!
      ndigit = -3
      logfil = 6
      mcaitr = 0 
      mcapps = 0
      mcaupd = 1
      mcaup2 = 0
      mceigh = 0
      mceupd = 0
!
!     %-------------------------------------------------%
!     | The following sets dimensions for this problem. |
!     %-------------------------------------------------%
!
      nx    = 10 
      n     = nx*nx 
!
!     %-----------------------------------------------%
!     |                                               | 
!     | Specifications for ARPACK usage are set       | 
!     | below:                                        |
!     |                                               |
!     |    1) NEV = 4  asks for 4 eigenvalues to be   |  
!     |       computed.                               | 
!     |                                               |
!     |    2) NCV = 20 sets the length of the Arnoldi |
!     |       factorization                           |
!     |                                               |
!     |    3) This is a standard problem              |
!     |         (indicated by bmat  = 'I')            |
!     |                                               |
!     |    4) Ask for the NEV eigenvalues of          |
!     |       largest magnitude                       |
!     |         (indicated by which = 'LM')           |
!     |       See documentation in ZNAUPD  for the     |
!     |       other options SM, LR, SR, LI, SI.       | 
!     |                                               |
!     | Note: NEV and NCV must satisfy the following  |
!     | conditions:                                   |
!     |              NEV <= MAXNEV                    |
!     |          NEV + 2 <= NCV <= MAXNCV             |
!     |                                               |
!     %-----------------------------------------------%
!
      nev   = 4
      ncv   = 20 
      bmat  = 'I'
      which = 'LM'
      
      if ( n .gt. maxn ) then
         print *, ' ERROR with _NSIMP: N is greater than MAXN '
         go to 9000
      else if ( nev .gt. maxnev ) then
         print *, ' ERROR with _NSIMP: NEV is greater than MAXNEV '
         go to 9000
      else if ( ncv .gt. maxncv ) then
         print *, ' ERROR with _NSIMP: NCV is greater than MAXNCV '
         go to 9000
      end if
!
!     %-----------------------------------------------------%
!     |                                                     |
!     | Specification of stopping rules and initial         |
!     | conditions before calling ZNAUPD                     |
!     |                                                     |
!     | TOL  determines the stopping criterion.             |
!     |                                                     |
!     |      Expect                                         |
!     |           abs(lambdaC - lambdaT) < TOL*abs(lambdaC) |
!     |               computed   true                       |
!     |                                                     |
!     |      If TOL .le. 0,  then TOL <- macheps            |
!     |           (machine precision) is used.              |
!     |                                                     |
!     | IDO  is the REVERSE COMMUNICATION parameter         |
!     |      used to specify actions to be taken on return  |
!     |      from ZNAUPD . (see usage below)                 |
!     |                                                     |
!     |      It MUST initially be set to 0 before the first |
!     |      call to ZNAUPD .                                | 
!     |                                                     |
!     | INFO on entry specifies starting vector information |
!     |      and on return indicates error codes            |
!     |                                                     |
!     |      Initially, setting INFO=0 indicates that a     | 
!     |      random starting vector is requested to         |
!     |      start the ARNOLDI iteration.  Setting INFO to  |
!     |      a nonzero value on the initial call is used    |
!     |      if you want to specify your own starting       |
!     |      vector (This vector must be placed in RESID).  | 
!     |                                                     |
!     | The work array WORKL is used in ZNAUPD  as           | 
!     | workspace.  Its dimension LWORKL is set as          |
!     | illustrated below.                                  |
!     |                                                     |
!     %-----------------------------------------------------%
!
      lworkl  = 3*ncv**2+5*ncv 
      tol    = 0.0 
      ido    = 0
      info   = 0
!
!     %---------------------------------------------------%
!     | Specification of Algorithm Mode:                  |
!     |                                                   |
!     | This program uses the exact shift strategy        |
!     | (indicated by setting IPARAM(1) = 1).             |
!     | IPARAM(3) specifies the maximum number of Arnoldi |
!     | iterations allowed.  Mode 1 of ZNAUPD  is used     |
!     | (IPARAM(7) = 1). All these options can be changed |
!     | by the user. For details see the documentation in |
!     | ZNAUPD .                                           |
!     %---------------------------------------------------%
!
      ishfts = 1
      maxitr = 300
      mode1 = 1
!
      iparam(1) = ishfts
!                
      iparam(3) = maxitr
!                  
      iparam(7) = mode1
!
!     %------------------------------------------------%
!     | M A I N   L O O P (Reverse Communication Loop) | 
!     %------------------------------------------------%
!
 10   continue
!   
!        %---------------------------------------------%
!        | Repeatedly call the routine ZNAUPD  and take | 
!        | actions indicated by parameter IDO until    |
!        | either convergence is indicated or maxitr   |
!        | has been exceeded.                          |
!        %---------------------------------------------%

         call znaupd  ( ido, bmat, n, which, nev, tol, resid, ncv,   &
     &                 v, ldv, iparam, ipntr, workd, workl, lworkl,  &
     &                 rwork,info )
!
         if (ido .eq. -1 .or. ido .eq. 1) then
!
!           %-------------------------------------------%
!           | Perform matrix vector multiplication      |
!           |                                           |
!           |                y <--- A*x                 |
!           |                                           |
!           | The user should supply his/her own        |
!           | matrix vector multiplication routine here |
!           | that takes workd(ipntr(1)) as the input   |
!           | vector x , and returns the resulting      |
!           | matrix-vector product y = A*x in the      |
!           | array workd(ipntr(2)).                    | 
!           %-------------------------------------------%
!
            call znsav (nx, workd(ipntr(1)), workd(ipntr(2)))
! 
!           %-----------------------------------------%
!           | L O O P   B A C K to call ZNAUPD  again. |
!           %-----------------------------------------%
!
            go to 10
!
         endif
! 
!     %----------------------------------------%
!     | Either we have convergence or there is |
!     | an error.                              |
!     %----------------------------------------%
!
      if ( info .lt. 0 ) then
!
!        %--------------------------%
!        | Error message, check the |
!        | documentation in ZNAUPD   |
!        %--------------------------%
!
         print *, ' '
         print *, ' Error with _naupd, info = ', info
         print *, ' Check the documentation of _naupd'
         print *, ' '
!
      else 
!
!        %-------------------------------------------%
!        | No fatal errors occurred.                 |
!        | Post-Process using ZNEUPD .                |
!        |                                           |
!        | Computed eigenvalues may be extracted.    |  
!        |                                           |
!        | Eigenvectors may be also computed now if  |
!        | desired.  (indicated by rvec = .true.)    | 
!        |                                           |
!        | The routine ZNEUPD  now called to do this  |
!        | post processing (Other modes may require  |
!        | more complicated post processing than     |
!        | mode1.)                                   |
!        |                                           |
!        %-------------------------------------------%
!           
         rvec = .true.
!
         call zneupd  (rvec, 'A', select, D, V, ldv, sigma,  &
     &        workev, bmat, n, which, nev, tol, resid, ncv,  &
     &        v, ldv, iparam, ipntr, workd, workl, lworkl,   &
     &        rwork, ierr)
!
!        %-----------------------------------------------%
!        | Eigenvalues are returned in the one           |
!        | dimensional array D and the corresponding     |
!        | eigenvectors are returned in the first        |
!        | NCONV (=IPARAM(5)) columns of the two         |
!        | dimensional array V if requested.  Otherwise, |
!        | an orthogonal basis for the invariant         |
!        | subspace corresponding to the eigenvalues in  |
!        | D is returned in V.                           |
!        %-----------------------------------------------%
!
         if ( ierr .ne. 0) then
!
!            %------------------------------------%
!            | Error condition:                   |
!            | Check the documentation of ZNEUPD . |
!            %------------------------------------%
!
             print *, ' '
             print *, ' Error with _neupd, info = ', ierr
             print *, ' Check the documentation of _neupd. '
             print *, ' '
!
         else
!
             nconv =  iparam(5)
             do 20 j=1, nconv
!
!                %---------------------------%
!                | Compute the residual norm |
!                |                           |
!                |   ||  A*x - lambda*x ||   |
!                |                           |
!                | for the NCONV accurately  |
!                | computed eigenvalues and  |
!                | eigenvectors.  (iparam(5) |
!                | indicates how many are    |
!                | accurate to the requested |
!                | tolerance)                |
!                %---------------------------%
!
                 call znsav(nx, v(1,j), ax)
                 call zaxpy (n, -d(j), v(1,j), 1, ax, 1)
                 rd(j,1) = dble (d(j))
                 rd(j,2) = dimag (d(j))
                 rd(j,3) = dznrm2 (n, ax, 1)
                 rd(j,3) = rd(j,3) / dlapy2 (rd(j,1),rd(j,2))
 20          continue
!
!            %-----------------------------%
!            | Display computed residuals. |
!            %-----------------------------%
!
             call dmout (6, nconv, 3, rd, maxncv, -6, 'Ritz values (Real, Imag) and relative residuals')
         end if
!
!        %-------------------------------------------%
!        | Print additional convergence information. |
!        %-------------------------------------------%
!
         if ( info .eq. 1) then
             print *, ' '
             print *, ' Maximum number of iterations reached.'
             print *, ' '
         else if ( info .eq. 3) then
             print *, ' ' 
             print *, ' No shifts could be applied during implicit',' Arnoldi update, try increasing NCV.'
             print *, ' '
         end if      
!
         print *, ' '
         print *, '_NSIMP '
         print *, '====== '
         print *, ' '
         print *, ' Size of the matrix is ', n
         print *, ' The number of Ritz values requested is ', nev
         print *, ' The number of Arnoldi vectors generated (NCV) is ', ncv
         print *, ' What portion of the spectrum: ', which
         print *, ' The number of converged Ritz values is ', nconv 
         print *, ' The number of Implicit Arnoldi update iterations taken is ', iparam(3)
         print *, ' The number of OP*x is ', iparam(9)
         print *, ' The convergence criterion is ', tol
         print *, ' '
!
      end if
!
!     %---------------------------%
!     | Done with program znsimp . |
!     %---------------------------%
!
 9000 continue
!
    END SUBROUTINE ZNSIMP_CASE1
    
! 
!==========================================================================
!
!     matrix vector subroutine
!
!     The matrix used is the convection-diffusion operator
!     discretized using centered difference.
!
      subroutine znsav (nx, v, w)
      integer           nx, j, lo
      Complex*16        v(nx*nx), w(nx*nx), one, h2
      parameter         (one = (1.0D+0, 0.0D+0) )
      external          zaxpy 
!
!     Computes w <--- OP*v, where OP is the nx*nx by nx*nx block 
!     tridiagonal matrix
!
!                  | T -I          | 
!                  |-I  T -I       |
!             OP = |   -I  T       |
!                  |        ...  -I|
!                  |           -I T|
!
!     derived from the standard central difference discretization 
!     of the 2-dimensional convection-diffusion operator 
!                  (Laplacian u) + rho*(du/dx)
!     on the unit squqre with zero boundary condition.
!
!     The subroutine znstv is called to computed y<---T*x.
!
!
      h2 = one / dcmplx ((nx+1)*(nx+1))
!
      call znstv(nx,v(1),w(1))
      call zaxpy (nx, -one/h2, v(nx+1), 1, w(1), 1)
!
      do 10 j = 2, nx-1
         lo = (j-1)*nx
         call znstv(nx, v(lo+1), w(lo+1))
         call zaxpy (nx, -one/h2, v(lo-nx+1), 1, w(lo+1), 1)
         call zaxpy (nx, -one/h2, v(lo+nx+1), 1, w(lo+1), 1)
  10  continue 
!
      lo = (nx-1)*nx
      call znstv(nx, v(lo+1), w(lo+1))
      call zaxpy (nx, -one/h2, v(lo-nx+1), 1, w(lo+1), 1)
!
      return
    end
    
    
    
!=========================================================================
      subroutine znstv (nx, x, y)
!
      integer           nx, j 
      Complex*16        x(nx), y(nx), h, h2, dd, dl, du
      Complex*16        one, rho
      parameter         (one = (1.0D+0, 0.0D+0) , rho = (1.0D+2, 0.0D+0) )
!
!     Compute the matrix vector multiplication y<---T*x
!     where T is a nx by nx tridiagonal matrix with DD on the 
!     diagonal, DL on the subdiagonal, and DU on the superdiagonal
!     
      h   = one / dcmplx (nx+1)
      h2  = h*h
      dd  = (4.0D+0, 0.0D+0)  / h2
      dl  = -one/h2 - (5.0D-1, 0.0D+0) *rho/h
      du  = -one/h2 + (5.0D-1, 0.0D+0) *rho/h
! 
      y(1) =  dd*x(1) + du*x(2)
      do 10 j = 2,nx-1
         y(j) = dl*x(j-1) + dd*x(j) + du*x(j+1) 
 10   continue 
      y(nx) =  dl*x(nx-1) + dd*x(nx) 
      return
      end    
  