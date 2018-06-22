      SUBROUTINE ZNDRV1  
!
!     Example program to illustrate the idea of reverse communication
!     for a standard complex nonsymmetric eigenvalue problem. 
!
!     We implement example one of ex-complex.doc in DOCUMENTS directory
!
!\Example-1
!     ... Suppose we want to solve A*x = lambda*x in regular mode,
!         where A is obtained from the standard central difference
!         discretization of the convection-diffusion operator 
!                 (Laplacian u) + rho*(du / dx)
!         on the unit squre [0,1]x[0,1] with zero Dirichlet boundary
!         condition.
!
!     ... OP = A  and  B = I.
!
!     ... Assume "call av (nx,x,y)" computes y = A*x
!
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
!     av      Matrix vector multiplication routine that computes A*x.
!     tv      Matrix vector multiplication routine that computes T*x,
!             where T is a tridiagonal matrix.  It is used in routine
!             av.
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
! FILE: ndrv1.F   SID: 2.4   DATE OF SID: 10/17/00   RELEASE: 2
!
!\Remarks
!     1. None
!
!\EndLib
!---------------------------------------------------------------------------
!
!     %-----------------------------%
!     | Define maximum dimensions   |
!     | for all arrays.             |
!     | MAXN:   Maximum dimension   |
!     |         of the A allowed.   |
!     | MAXNEV: Maximum NEV allowed |
!     | MAXNCV: Maximum NCV allowed |
!     %-----------------------------%
!
      integer           maxn, maxnev, maxncv, ldv
      parameter         (maxn=256, maxnev=12, maxncv=30, ldv=maxn)
!
!     %--------------%
!     | Local Arrays |
!     %--------------%
!
      integer           iparam(11), ipntr(14)
      logical           select(maxncv)
      Complex*16        ax(maxn), d(maxncv) 
      Complex*16        v(ldv,maxncv), workd(3*maxn) 
      Complex*16         workev(3*maxncv), resid(maxn) 
      Complex*16         workl(3*maxncv*maxncv+5*maxncv)
      Double precision  rwork(maxncv), rd(maxncv,3)
!
!     %---------------%
!     | Local Scalars |
!     %---------------%
!
      character         bmat*1, which*2
      integer           ido, n, nx, nev, ncv, lworkl, info, j, &
     &                  ierr, nconv, maxitr, ishfts, mode
      Complex*16        sigma
      Double precision  tol
      logical           rvec
!
!     %-----------------------------%
!     | BLAS & LAPACK routines used |
!     %-----------------------------%
!
      Double precision  dznrm2 , dlapy2 
      external          dznrm2 , zaxpy , dlapy2  
!
!     %-----------------------%
!     | Executable Statements |
!     %-----------------------%
! 
!     %--------------------------------------------------%
!     | The number NX is the number of interior points   |
!     | in the discretization of the 2-dimensional       |
!     | convection-diffusion operator on the unit        |
!     | square with zero Dirichlet boundary condition.   | 
!     | The number N(=NX*NX) is the dimension of the     |
!     | matrix.  A standard eigenvalue problem is        |
!     | solved (BMAT = 'I').  NEV is the number of       |
!     | eigenvalues to be approximated.  The user can    |
!     | modify NX, NEV, NCV, WHICH to solve problems of  |
!     | different sizes, and to get different parts of   |
!     | the spectrum.  However, The following            |
!     | conditions must be satisfied:                    |
!     |                   N <= MAXN                      |
!     |                 NEV <= MAXNEV                    |
!     |           NEV + 2 <= NCV <= MAXNCV               | 
!     %--------------------------------------------------% 
!
      nx    = 10 
      n     = nx*nx 
      nev   = 4
      ncv   = 20 
      if ( n .gt. maxn ) then
         print *, ' ERROR with _NDRV1: N is greater than MAXN '
         go to 9000
      else if ( nev .gt. maxnev ) then
         print *, ' ERROR with _NDRV1: NEV is greater than MAXNEV '
         go to 9000
      else if ( ncv .gt. maxncv ) then
         print *, ' ERROR with _NDRV1: NCV is greater than MAXNCV '
         go to 9000
      end if
      bmat  = 'I'
      which = 'LM'
!
!     %---------------------------------------------------%
!     | The work array WORKL is used in ZNAUPD  as         | 
!     | workspace.  Its dimension LWORKL is set as        |
!     | illustrated below.  The parameter TOL determines  |
!     | the stopping criterion. If TOL<=0, machine        |
!     | precision is used.  The variable IDO is used for  |
!     | reverse communication, and is initially set to 0. |
!     | Setting INFO=0 indicates that a random vector is  |
!     | generated to start the ARNOLDI iteration.         | 
!     %---------------------------------------------------%
!
      lworkl  = 3*ncv**2+5*ncv 
      tol    = 0.0 
      ido    = 0
      info   = 0
!
!     %---------------------------------------------------%
!     | This program uses exact shift with respect to     |
!     | the current Hessenberg matrix (IPARAM(1) = 1).    |
!     | IPARAM(3) specifies the maximum number of Arnoldi |
!     | iterations allowed.  Mode 1 of ZNAUPD  is used     |
!     | (IPARAM(7) = 1). All these options can be changed |
!     | by the user. For details see the documentation in |
!     | ZNAUPD .                                           |
!     %---------------------------------------------------%
!
      ishfts = 1
      maxitr = 300
      mode   = 1
!
      iparam(1) = ishfts
      iparam(3) = maxitr 
      iparam(7) = mode 
!
!     %-------------------------------------------%
!     | M A I N   L O O P (Reverse communication) | 
!     %-------------------------------------------%
!
 10   continue
!
!        %---------------------------------------------%
!        | Repeatedly call the routine ZNAUPD  and take |
!        | actions indicated by parameter IDO until    |
!        | either convergence is indicated or maxitr   |
!        | has been exceeded.                          |
!        %---------------------------------------------%
!
         call znaupd  ( ido, bmat, n, which, nev, tol, resid, ncv, &
     &        v, ldv, iparam, ipntr, workd, workl, lworkl,         &
     &        rwork,info )
!
         if (ido .eq. -1 .or. ido .eq. 1) then
!
!           %-------------------------------------------%
!           | Perform matrix vector multiplication      |
!           |                y <--- OP*x                |
!           | The user should supply his/her own        |
!           | matrix vector multiplication routine here |
!           | that takes workd(ipntr(1)) as the input   |
!           | vector, and return the matrix vector      |
!           | product to workd(ipntr(2)).               | 
!           %-------------------------------------------%
!
            call av (nx, workd(ipntr(1)), workd(ipntr(2)))
!
!           %-----------------------------------------%
!           | L O O P   B A C K to call ZNAUPD  again. |
!           %-----------------------------------------%
!
            go to 10

         end if
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
!        | Eigenvectors may also be computed now if  |
!        | desired.  (indicated by rvec = .true.)    |
!        %-------------------------------------------%
!
         rvec = .true.
!
         call zneupd  (rvec, 'A', select, d, v, ldv, sigma, &
     &        workev, bmat, n, which, nev, tol, resid, ncv, & 
     &        v, ldv, iparam, ipntr, workd, workl, lworkl,  &
     &        rwork, ierr)
!
!        %----------------------------------------------%
!        | Eigenvalues are returned in the one          |
!        | dimensional array D.  The corresponding      |
!        | eigenvectors are returned in the first NCONV |
!        | (=IPARAM(5)) columns of the two dimensional  | 
!        | array V if requested.  Otherwise, an         |
!        | orthogonal basis for the invariant subspace  |
!        | corresponding to the eigenvalues in D is     |
!        | returned in V.                               |
!        %----------------------------------------------%
!
         if ( ierr .ne. 0) then
! 
!           %------------------------------------%
!           | Error condition:                   |
!           | Check the documentation of ZNEUPD . |
!           %------------------------------------%
!
             print *, ' '
             print *, ' Error with _neupd, info = ', ierr
             print *, ' Check the documentation of _neupd. '
             print *, ' '
!
         else
!
             nconv = iparam(5)
             do 20 j=1, nconv
!
!               %---------------------------%
!               | Compute the residual norm |
!               |                           |
!               |   ||  A*x - lambda*x ||   |
!               |                           |
!               | for the NCONV accurately  |
!               | computed eigenvalues and  |
!               | eigenvectors.  (iparam(5) |
!               | indicates how many are    |
!               | accurate to the requested |
!               | tolerance)                |
!               %---------------------------%
!
                call av(nx, v(1,j), ax)
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
             call dmout (6, nconv, 3, rd, maxncv, -6,'Ritz values (Real, Imag) and relative residuals')
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
             print *, ' No shifts could be applied during implicit Arnoldi update, try increasing NCV.'
             print *, ' '
         end if      
!
         print *, ' '
         print *, '_NDRV1'
         print *, '====== '
         print *, ' '
         print *, ' Size of the matrix is ', n
         print *, ' The number of Ritz values requested is ', nev
         print *, ' The number of Arnoldi vectors generated (NCV) is ', ncv
         print *, ' What portion of the spectrum: ', which
         print *, ' The number of converged Ritz values is ',  nconv 
         print *, ' The number of Implicit Arnoldi update iterations taken is ', iparam(3)
         print *, ' The number of OP*x is ', iparam(9)
         print *, ' The convergence criterion is ', tol
         print *, ' '
!
     end if
!
!     %---------------------------%
!     | Done with program zndrv1 . |
!     %---------------------------%
!
 9000 continue
!
END SUBROUTINE ZNDRV1
! 
!==========================================================================
!
!     matrix vector subroutine
!
!     The matrix used is the convection-diffusion operator
!     discretized using centered difference.
!
      subroutine av (nx, v, w)
      integer           nx, j, lo
      Complex*16        v(nx*nx), w(nx*nx), one, h2
      parameter         (one = (1.0D+0, 0.0D+0) )
      external          zaxpy , tv
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
!     derived from the standard central difference  discretization 
!     of the convection-diffusion operator (Laplacian u) + rho*(du/dx)
!     with zero boundary condition.
!
!     The subroutine TV is called to computed y<---T*x.
!
!
      h2 = one / dcmplx ((nx+1)*(nx+1))
!
      call tv(nx,v(1),w(1))
      call zaxpy (nx, -one/h2, v(nx+1), 1, w(1), 1)
!
      do 10 j = 2, nx-1
         lo = (j-1)*nx
         call tv(nx, v(lo+1), w(lo+1))
         call zaxpy (nx, -one/h2, v(lo-nx+1), 1, w(lo+1), 1)
         call zaxpy (nx, -one/h2, v(lo+nx+1), 1, w(lo+1), 1)
  10  continue 
!
      lo = (nx-1)*nx
      call tv(nx, v(lo+1), w(lo+1))
      call zaxpy (nx, -one/h2, v(lo-nx+1), 1, w(lo+1), 1)
!
      return
      end
!=========================================================================
      subroutine tv (nx, x, y)
!
      integer           nx, j 
      Complex*16        x(nx), y(nx), h, h2, dd, dl, du
!
      Complex*16        one, rho
      parameter         (one = (1.0D+0, 0.0D+0) ,rho = (1.0D+2, 0.0D+0) )
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