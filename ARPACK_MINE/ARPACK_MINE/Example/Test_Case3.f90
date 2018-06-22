      SUBROUTINE ZNDRV2 
!
!     Simple program to illustrate the idea of reverse communication
!     in shift-invert mode for a standard complex nonsymmetric eigenvalue 
!     problem.
!
!     We implement example two of ex-complex.doc in DOCUMENTS directory
!
!\Example-2
!     ... Suppose we want to solve A*x = lambda*x in shift-invert mode,
!         where A is derived from the central difference discretization
!         of the 1-dimensional convection-diffusion operator
!                   (d^2u/dx^2) + rho*(du/dx)
!         on the interval [0,1] with zero Dirichlet boundary condition.
!     ... The shift sigma is a complex number.
!
!     ... OP = inv[A-sigma*I] and  B = I.
!
!     ... Use mode 3 of ZNAUPD .
!
!\BeginLib
!
!\Routines called:
!     znaupd   ARPACK reverse communication interface routine.
!     zneupd   ARPACK routine that returns Ritz values and (optionally)
!             Ritz vectors.
!     zgttrf   LAPACK tridiagonal factorization routine.
!     zgttrs   LAPACK tridiagonal solve routine.
!     dlapy2   LAPACK routine to compute sqrt(x**2+y**2) carefully.
!     zaxpy    Level 1 BLAS that computes y <- alpha*x+y.
!     zcopy    Level 1 BLAS that copies one vector to another.
!     dznrm2   Level 1 BLAS that computes the norm of a vector.
!     zndrv2av Matrix vector multiplication routine that computes A*x.
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
! FILE: ndrv2.F   SID: 2.6   DATE OF SID: 10/18/00   RELEASE: 2
!
!\Remarks
!     1. None
!
!\EndLib
!--------------------------------------------------------------------------
!
!     %-----------------------------%
!     | Define leading dimensions   |
!     | for all arrays.             |
!     | MAXN:   Maximum dimension   |
!     |         of the A allowed.   |
!     | MAXNEV: Maximum NEV allowed |
!     | MAXNCV: Maximum NCV allowed |
!     %-----------------------------%
!
      integer           maxn, maxnev, maxncv, ldv
      parameter         (maxn=256, maxnev=10, maxncv=25, ldv=maxn )
!
!     %--------------%
!     | Local Arrays |
!     %--------------%
!
      integer           iparam(11), ipntr(14), ipiv(maxn)
      logical           select(maxncv)
      Complex*16        ax(maxn), d(maxncv), resid(maxn),  &
     &                  v(ldv, maxncv), workd(3*maxn),     &
     &                  workev(2*maxncv),                  &
     &                  workl(3*maxncv*maxncv+5*maxncv),   &
     &                  dd(maxn), dl(maxn), du(maxn),      &
     &                  du2(maxn)
      Double precision  rwork(maxn), rd(maxncv,3)
!
!     %---------------%
!     | Local Scalars |
!     %---------------%
!
      character         bmat*1, which*2
      integer           ido, n, nev, ncv, lworkl, info, j, ierr, &
     &                  nconv, maxitr, ishfts, mode
      Complex*16        h, h2, s, sigma, s1, s2, s3, rho
!
      Double precision  tol
      logical           rvec
!
!     %------------%
!     | Parameters |
!     %------------%
!
      Complex*16         one, zero, two        
      parameter         (one = (1.0D+0, 0.0D+0) , zero = (0.0D+0, 0.0D+0) ,two = (2.0D+0, 0.0D+0) )
!
!     %-----------------------------%
!     | BLAS & LAPACK routines used |
!     %-----------------------------%
!
      Double precision  dznrm2 , dlapy2 
      external          zgttrf , zgttrs , zaxpy , zcopy , dznrm2 ,dlapy2 
!
!     %-----------------------%
!     | Executable statements |
!     %-----------------------%
!
!     %--------------------------------------------------%
!     | The number N is the dimension of the matrix.  A  |
!     | standard eigenvalue problem is solved (BMAT =    |
!     | 'I').  NEV is the number of eigenvalues (closest |
!     | to the shift SIGMA) to be approximated.  Since   |
!     | the shift-invert mode is used, WHICH is set to   |
!     | 'LM'.  The user can modify NEV, NCV, SIGMA to    |
!     | solve problems of different sizes, and to get    |
!     | different parts of the spectrum.  However, The   |
!     | following conditions must be satisfied:          |
!     |                 N <= MAXN,                       |
!     |               NEV <= MAXNEV,                     |
!     |           NEV + 2 <= NCV <= MAXNCV               |
!     %--------------------------------------------------%

!
      n     = 100
      nev   = 4
      ncv   = 20
      if ( n .gt. maxn ) then
         print *, ' ERROR with _NDRV2: N is greater than MAXN '
         go to 9000
      else if ( nev .gt. maxnev ) then
         print *, ' ERROR with _NDRV2: NEV is greater than MAXNEV '
         go to 9000
      else if ( ncv .gt. maxncv ) then
         print *, ' ERROR with _NDRV2: NCV is greater than MAXNCV '
         go to 9000
      end if
      bmat  = 'I'
      which = 'LM'
      sigma = zero
!
!     %----------------------------------------------------%
!     | Construct C = A - SIGMA*I, factor C in complex     |
!     | arithmetic (using LAPACK subroutine zgttrf ). The   |
!     | matrix A is chosen to be the tridiagonal matrix    |
!     | derived from standard central difference of the    |
!     | 1-d convection diffusion operator - u``+ rho*u` on |
!     | the interval [0, 1] with zero Dirichlet boundary   |
!     | condition.                                         |
!     %----------------------------------------------------%
!
      rho = (1.0D+1, 0.0D+0) 
      h = one / dcmplx (n+1)
      h2 = h*h
      s = rho / two
!
      s1 = -one/h2 - s/h 
      s2 = two/h2 - sigma
      s3 = -one/h2 + s/h  
!
      do 10 j = 1, n-1
         dl(j) = s1 
         dd(j) = s2
         du(j) = s3
  10  continue 
      dd(n) = s2 
! 
      call zgttrf (n, dl, dd, du, du2, ipiv, ierr)
      if ( ierr .ne. 0 ) then
         print*, ' '
         print*, ' ERROR with _gttrf in _NDRV2.'
         print*, ' '
         go to 9000
      end if
!
!     %-----------------------------------------------------%
!     | The work array WORKL is used in ZNAUPD  as           |
!     | workspace.  Its dimension LWORKL is set as          |
!     | illustrated below.  The parameter TOL determines    |
!     | the stopping criterion. If TOL<=0, machine          |
!     | precision is used.  The variable IDO is used for    |
!     | reverse communication, and is initially set to 0.   |
!     | Setting INFO=0 indicates that a random vector is    |
!     | generated in ZNAUPD  to start the Arnoldi iteration. |
!     %-----------------------------------------------------%
!
      lworkl = 3*ncv**2+5*ncv 
      tol    = 0.0 
      ido    = 0
      info   = 0
!
!     %---------------------------------------------------%
!     | This program uses exact shifts with respect to    |
!     | the current Hessenberg matrix (IPARAM(1) = 1).    |
!     | IPARAM(3) specifies the maximum number of Arnoldi |
!     | iterations allowed. Mode 3 of ZNAUPD  is used      |
!     | (IPARAM(7) = 3).  All these options can be        |
!     | changed by the user. For details see the          |
!     | documentation in ZNAUPD .                          |
!     %---------------------------------------------------%
!
      ishfts = 1
      maxitr = 300
      mode   = 3
!
      iparam(1) = ishfts 
      iparam(3) = maxitr 
      iparam(7) = mode 
!
!     %-------------------------------------------%
!     | M A I N   L O O P (Reverse communication) | 
!     %-------------------------------------------%
!
 20   continue
!
!        %---------------------------------------------%
!        | Repeatedly call the routine ZNAUPD  and take | 
!        | actions indicated by parameter IDO until    |
!        | either convergence is indicated or maxitr   |
!        | has been exceeded.                          |
!        %---------------------------------------------%
!
         call znaupd  ( ido, bmat, n, which, nev, tol, resid, ncv,  &
     &        v, ldv, iparam, ipntr, workd, workl, lworkl,          &
     &        rwork,info )
!
         if (ido .eq. -1 .or. ido .eq. 1 ) then
!
!           %-------------------------------------------%
!           | Perform  y <--- OP*x = inv[A-SIGMA*I]*x   |
!           | The user should supply his/her own linear |
!           | system solver here that takes             |
!           | workd(ipntr(1)) as the input, and returns |
!           | the result to workd(ipntr(2)).            |
!           %-------------------------------------------%
!
            call zcopy ( n, workd(ipntr(1)),1, workd(ipntr(2)), 1)
!
            call zgttrs ('N', n, 1, dl, dd, du, du2, ipiv, workd(ipntr(2)), n, ierr) 
            if ( ierr .ne. 0 ) then
               print*, ' '
               print*, ' ERROR with _gttrs in _NDRV2.'
               print*, ' '
               go to 9000
            end if
!
!           %-----------------------------------------%
!           | L O O P   B A C K to call ZNAUPD  again. |
!           %-----------------------------------------%
!
            go to 20
!
         end if
!
!     %-----------------------------------------%
!     | Either we have convergence, or there is |
!     | an error.                               |
!     %-----------------------------------------%
!
      if ( info .lt. 0 ) then
!
!        %--------------------------%
!        | Error message, check the |
!        | documentation in ZNAUPD   |
!        %--------------------------%
!
         print *, ' '
         print *, ' Error with _naupd, info = ',info
         print *, ' Check the documentation in _naupd.'
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
         call zneupd  (rvec, 'A', select, d, v, ldv, sigma,     &
     &                workev, bmat, n, which, nev, tol,         &
     &                resid, ncv, v, ldv, iparam, ipntr, workd, &
     &                workl, lworkl, rwork, ierr)
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
         else
!
            nconv = iparam(5) 
            do 60 j=1, nconv
!
!              %---------------------------%
!              | Compute the residual norm |
!              |                           |
!              |   ||  A*x - lambda*x ||   |
!              |                           |
!              | for the NCONV accurately  |
!              | computed eigenvalues and  |
!              | eigenvectors.  (iparam(5) |
!              | indicates how many are    |
!              | accurate to the requested |
!              | tolerance)                |
!              %---------------------------%
!
               call zndrv2av(n, v(1,j), ax, rho)
               call zaxpy (n, -d(j), v(1,j), 1, ax, 1)
               rd(j,1) = dble (d(j))
               rd(j,2) = dimag (d(j))
               rd(j,3) = dznrm2 (n, ax, 1)
               rd(j,3) = rd(j,3) / dlapy2 (rd(j,1),rd(j,2))
 60         continue
!
!           %-----------------------------%
!           | Display computed residuals. |
!           %-----------------------------%
!
            call dmout (6, nconv, 3, rd, maxncv, -6, 'Ritz values (Real, Imag) and relative residuals')
!
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
             print *, ' No shifts could be applied during implicit'
             print *, ' Arnoldi update, try increasing NCV.'
             print *, ' '
         end if      
!
         print *, ' '
         print *, '_NDRV2 '
         print *, '====== '
         print *, ' '
         print *, ' Size of the matrix is ', n
         print *, ' The number of Ritz values requested is ', nev
         print *, ' The number of Arnoldi vectors generated (NCV) is ', ncv
         print *, ' What portion of the spectrum: ', which
         print *, ' The number of converged Ritz values is ',nconv 
         print *, ' The number of Implicit Arnoldi update iterations taken is ', iparam(3)
         print *, ' The number of OP*x is ', iparam(9)
         print *, ' The convergence criterion is ', tol
         print *, ' '
!
      end if
!
!     %---------------------------%
!     | Done with program zndrv2 . |
!     %---------------------------%
!
 9000 continue
!
      END SUBROUTINE ZNDRV2
!
!-------------------------------------------------------------------
!
!     matrix vector multiplication subroutine
!
      subroutine zndrv2av (n, v, w, rho)
      integer           n, j
      Complex*16        v(n), w(n), rho, two, one, dd, dl, du, s, h,h2
      parameter         (one = (1.0D+0, 0.0D+0) ,two = (2.0D+0, 0.0D+0) )
      
      h = one / dcmplx (n+1)
      h2 = h*h
      s = rho / two
      dd = two / h2
      dl = -one/h2 - s/h
      du = -one/h2 + s/h

      w(1) =  dd*v(1) + du*v(2)
      do 10 j = 2,n-1
         w(j) = dl*v(j-1) + dd*v(j) + du*v(j+1)
 10   continue
      w(n) =  dl*v(n-1) + dd*v(n)
      return
      end
