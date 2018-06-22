      SUBROUTINE ZNDRV3 
!
!     Simple program to illustrate the idea of reverse communication
!     in inverse mode for a generalized complex nonsymmetric eigenvalue 
!     problem.
!
!     We implement example three of ex-complex.doc in DOCUMENTS directory
!
!\Example-3
!     ... Suppose we want to solve A*x = lambda*B*x in regular mode,
!         where A and B are derived from the finite element discretization
!         of the 1-dimensional convection-diffusion operator
!                   (d^2u/dx^2) + rho*(du/dx)
!         on the interval [0,1] with zero boundary condition using 
!         piecewise linear elements. 
!
!     ... OP = inv[M]*A  and  B = M.
!
!     ... Use mode 2 of ZNAUPD .
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
!     dznrm2   Level 1 BLAS that computes the norm of a vector.
!     zndrv3av      Matrix vector multiplication routine that computes A*x.
!     zndrv3mv      Matrix vector multiplication routine that computes M*x.
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
! FILE: ndrv3.F   SID: 2.4   DATE OF SID: 10/18/00   RELEASE: 2
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
      Complex*16        ax(maxn), mx(maxn), d(maxncv), resid(maxn), &
     &                  v(ldv,maxncv), workd(3*maxn),               &
     &                  workev(2*maxncv),                           &
     &                  workl(3*maxncv*maxncv+5*maxncv),            & 
     &                  dd(maxn), dl(maxn), du(maxn), du2(maxn)
      Double precision  rwork(maxn), rd(maxncv,3)
!
!     %---------------%
!     | Local Scalars |
!     %---------------%
!
      character         bmat*1, which*2
      integer           ido, n, nev, ncv, lworkl, info, ierr, j,    &
     &                  nconv, maxitr, ishfts, mode
      Complex*16        sigma, h
      Double precision  tol
      logical           rvec
!
!     %------------%
!     | Parameters |
!     %------------%
!
      Complex*16        zero, one
      parameter         (zero = (0.0D+0, 0.0D+0) ,one = (1.0D+0, 0.0D+0) )
!
!     %-----------------------------%
!     | BLAS & LAPACK routines used |
!     %-----------------------------%
!
      Double precision  dznrm2 , dlapy2 
      external          zaxpy , zcopy , dznrm2 , zgttrf , zgttrs ,dlapy2 
!
!     %-----------------------%
!     | Executable Statements |
!     %-----------------------%
!
!     %----------------------------------------------------%
!     | The number N is the dimension of the matrix.  A    |
!     | generalized eigenvalue problem is solved (BMAT =   |
!     | 'G').  NEV is the number of eigenvalues to be      |
!     | approximated.  The user can modify NEV, NCV, WHICH |
!     | to solve problems of different sizes, and to get   |
!     | different parts of the spectrum.  However, The     |
!     | following conditions must be satisfied:            |
!     |                    N <= MAXN,                      |
!     |                  NEV <= MAXNEV,                    |
!     |              NEV + 2 <= NCV <= MAXNCV              |
!     %----------------------------------------------------%
!
      n     = 100 
      nev   = 4 
      ncv   = 20 
      if ( n .gt. maxn ) then
         print *, ' ERROR with _NDRV3: N is greater than MAXN '
         go to 9000
      else if ( nev .gt. maxnev ) then
         print *, ' ERROR with _NDRV3: NEV is greater than MAXNEV '
         go to 9000
      else if ( ncv .gt. maxncv ) then
         print *, ' ERROR with _NDRV3: NCV is greater than MAXNCV '
         go to 9000
      end if
      bmat  = 'G'
      which = 'LM'
      sigma = zero
!
!     %-----------------------------------------------------%
!     | The matrix M is chosen to be the symmetric tri-     |
!     | diagonal matrix with 4 on the diagonal and 1 on the |
!     | off diagonals. It is factored by LAPACK subroutine  |
!     | zgttrf .                                             |
!     %-----------------------------------------------------%
!
      h = one / dcmplx (n+1)
      do 20 j = 1, n-1
         dl(j) = one*h
         dd(j) = (4.0D+0, 0.0D+0) *h
         du(j) = one*h 
  20  continue 
      dd(n) = (4.0D+0, 0.0D+0) *h 
!
      call zgttrf (n, dl, dd, du, du2, ipiv, ierr) 
      if ( ierr .ne. 0 ) then
         print*, ' '
         print*, ' ERROR with _gttrf. ' 
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
!     | iterations allowed.  Mode 2 of ZNAUPD  is used     |
!     | (IPARAM(7) = 2).  All these options can be        |
!     | changed by the user. For details, see the         |
!     | documentation in ZNAUPD .                          |
!     %---------------------------------------------------%
!
      ishfts = 1
      maxitr = 300
      mode   = 2
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
     &        rwork, info )
!
         if (ido .eq. -1 .or. ido .eq. 1) then
!
!           %----------------------------------------%
!           | Perform  y <--- OP*x = inv[M]*A*x      |
!           | The user should supply his/her own     |
!           | matrix vector routine and a linear     |
!           | system solver.  The matrix-vector      |
!           | subroutine should take workd(ipntr(1)) |
!           | as input, and the final result should  |
!           | be returned to workd(ipntr(2)).        |
!           %----------------------------------------%
!
            call zndrv3av (n, workd(ipntr(1)), workd(ipntr(2)))
            call zgttrs ('N', n, 1, dl, dd, du, du2, ipiv,workd(ipntr(2)), n, ierr)
            if ( ierr .ne. 0 ) then
               print*, ' ' 
               print*, ' ERROR with _gttrs. ' 
               print*, ' '
               go to 9000
            end if
!
!           %-----------------------------------------%
!           | L O O P   B A C K to call ZNAUPD  again. |
!           %-----------------------------------------%
!
            go to 10
!
         else if ( ido .eq. 2) then
!
!           %-------------------------------------%
!           |        Perform  y <--- M*x          |
!           | The matrix vector multiplication    |
!           | routine should take workd(ipntr(1)) |
!           | as input and return the result to   |
!           | workd(ipntr(2)).                    |
!           %-------------------------------------%
!
            call zndrv3mv (n, workd(ipntr(1)), workd(ipntr(2)))
!
!           %-----------------------------------------%
!           | L O O P   B A C K to call ZNAUPD  again. |
!           %-----------------------------------------%
!
            go to 10
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
!        | Error message. Check the |
!        | documentation in ZNAUPD . |
!        %--------------------------%
!
         print *, ' '
         print *, ' Error with _naupd, info = ',info
         print *, ' Check the documentation of _naupd.'
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
         call zneupd  ( rvec, 'A', select, d, v, ldv, sigma,    &
     &        workev, bmat, n, which, nev, tol, resid, ncv, v,  &
     &        ldv, iparam, ipntr, workd, workl, lworkl, rwork,  &
     &        ierr )
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
         if ( ierr .ne. 0 ) then
!
!           %------------------------------------%
!           | Error condition:                   |
!           | Check the documentation of ZNEUPD . |
!           %------------------------------------%
! 
            print *, ' ' 
            print *, ' Error with _neupd, info = ', ierr
            print *, ' Check the documentation of _neupd'
            print *, ' '
!
         else
!
            nconv = iparam(5)
            do 80 j=1, nconv 
!
!              %---------------------------%
!              | Compute the residual norm |
!              |                           |
!              |  ||  A*x - lambda*M*x ||  |
!              |                           |
!              | for the NCONV accurately  |
!              | computed eigenvalues and  |
!              | eigenvectors.  (iparam(5) |
!              | indicates how many are    |
!              | accurate to the requested |
!              | tolerance)                |
!              %---------------------------%
!
               call zndrv3av(n, v(1,j), ax)
               call zndrv3mv(n, v(1,j), mx)
               call zaxpy (n, -d(j), mx, 1, ax, 1)
               rd(j,1) = dble (d(j))
               rd(j,2) = dimag (d(j))
               rd(j,3) = dznrm2 (n, ax, 1)
               rd(j,3) = rd(j,3) / dlapy2 (rd(j,1),rd(j,2))
  80        continue
!
!           %-----------------------------%
!           | Display computed residuals. |
!           %-----------------------------%
!
            call dmout (6, nconv, 3, rd, maxncv, -6, 'Ritz values (Real, Imag) and relative residuals')
!
         end if
!
!        %------------------------------------------%
!        | Print additional convergence information |
!        %------------------------------------------%
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
         print *, '_NDRV3 '
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
 9000 continue
!
      end
! 
!==========================================================================
!
!     matrix vector multiplication subroutine
!
      subroutine zndrv3av (n, v, w)
      integer           n, j
      Complex*16        v(n), w(n), one, two, dd, dl, du, s, h, rho 
      parameter         (one = (1.0D+0, 0.0D+0) ,two = (2.0D+0, 0.0D+0) ,rho = (1.0D+1, 0.0D+0) )
!
!     Compute the matrix vector multiplication y<---A*x
!     where A is the stiffness matrix formed by using piecewise linear 
!     elements on [0,1].
!
      h = one / dcmplx (n+1)
      s = rho / two
      dd = two / h
      dl = -one/h - s
      du = -one/h + s
!
      w(1) =  dd*v(1) + du*v(2)
      do 10 j = 2,n-1
         w(j) = dl*v(j-1) + dd*v(j) + du*v(j+1) 
 10   continue 
      w(n) =  dl*v(n-1) + dd*v(n) 
      return
      end subroutine zndrv3av
!------------------------------------------------------------------------
      subroutine zndrv3mv (n, v, w)
      integer           n, j
      Complex*16        v(n), w(n), one, four, h
      parameter         (one = (1.0D+0, 0.0D+0) , four = (4.0D+0, 0.0D+0) )
!
!     Compute the matrix vector multiplication y<---M*x
!     where M is the mass matrix formed by using piecewise linear elements 
!     on [0,1].
! 
      w(1) =  four*v(1) + one*v(2)
      do 10 j = 2,n-1
         w(j) = one*v(j-1) + four*v(j) + one*v(j+1) 
 10   continue 
      w(n) =  one*v(n-1) + four*v(n) 
!
      h = one / dcmplx (n+1)
      call zscal (n, h, w, 1)
      return
      end subroutine zndrv3mv
