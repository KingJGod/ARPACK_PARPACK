 subroutine pzndrv1
    use DEBUG_MODULE
    use STAT_MODULE
    IMPLICIT NONE
    include 'mpif.h'
    
! 
!     %---------------%
!     | MPI INTERFACE |
!     %---------------%
 
      INTEGER( kind = 4 )::comm, myid, nprocs, rc, nloc
!     %-----------------------------%
!     | Define maximum dimensions   |
!     | for all arrays.             |
!     | MAXN:   Maximum dimension   |
!     |         of the A allowed.   |
!     | MAXNEV: Maximum NEV allowed |
!     | MAXNCV: Maximum NCV allowed |
!     %-----------------------------%
!
      INTEGER( kind = 4 )::maxn, maxnev, maxncv, ldv
      parameter (maxn=2560, maxnev=100, maxncv=300, ldv=maxn)
!
!     %--------------%
!     | Local Arrays |
!     %--------------%
!
      INTEGER( kind = 4 )::iparam(11), ipntr(14)
      logical::select(maxncv)
      Complex( kind = 8 )::ax(maxn), d(maxncv),             &
                        &  v(ldv,maxncv), workd(3*maxn),    &
                        &  workev(3*maxncv), resid(maxn),   &
                        &  workl(3*maxncv*maxncv+5*maxncv)
      REAL( kind = 8 )::rwork(maxncv), rd(maxncv,3)
!
!     %---------------%
!     | Local Scalars |
!     %---------------%
!
      character         bmat*1, which*2
      INTEGER( kind = 4 )::ido, n, nx, nev, ncv, lworkl, info, j
      INTEGER( kind = 4 )::ierr, nconv, maxitr, ishfts, mode
      Complex( kind = 8 )::sigma
      REAL( kind = 8 )::tol
      logical::rvec
!
!     %----------------------------------------------%
!     | Local Buffers needed for MPI communication |
!     %----------------------------------------------%
!
      Complex( kind = 8)::mv_buf(maxn)
!
!     %-----------------------------%
!     | BLAS & LAPACK routines used |
!     %-----------------------------%
!
      REAL( kind = 8 )::pdznorm2
      external          pdznorm2, zaxpy 
!
!     %-----------------------%
!     | Executable Statements |
!     %-----------------------%
! 

      comm = MPI_COMM_WORLD
      call MPI_COMM_RANK( comm, myid, ierr )
      call MPI_COMM_SIZE( comm, nprocs, ierr )

!     These parameter control the output information and do not have
!     influence in main calculations procedual
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
      nev   = 6
      ncv   = 18 
!
!     %--------------------------------------%
!     | Set up distribution of data to nodes |
!     %--------------------------------------%
!
      nloc = (nx / nprocs)*nx
      if ( mod(nx, nprocs) .gt. myid ) nloc = nloc + nx
     
      write(*,*)'myid:',myid,'nloc:',nloc
      
      if ( nloc .gt. maxn ) then
         print *, ' ERROR with _NDRV1: NLOC is greater than MAXN '
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
!     | The work array WORKL is used in ZNAUPD as         | 
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
!     | iterations allowed.  Mode 1 of ZNAUPD is used     |
!     | (IPARAM(7) = 1). All these options can be changed |
!     | by the user. For details see the documentation in |
!     | ZNAUPD.                                           |
!     %---------------------------------------------------%
!
      ishfts = 1
      maxitr = 300
      mode   = 1

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
!        | Repeatedly call the routine ZNAUPD and take |
!        | actions indicated by parameter IDO until    |
!        | either convergence is indicated or maxitr   |
!        | has been exceeded.                          |
!        %---------------------------------------------%
!
         call pznaupd ( comm, ido, bmat, nloc, which,      &
     &        nev, tol, resid, ncv, v, ldv, iparam, ipntr, & 
     &        workd, workl, lworkl, rwork,info )

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
            call pzndrv1av ( comm, nloc, nx, mv_buf,           &
     &                workd(ipntr(1)), workd(ipntr(2)))
!
!           %-----------------------------------------%
!           | L O O P   B A C K to call ZNAUPD again. |
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
!        | documentation in ZNAUPD  |
!        %--------------------------%
!
         if ( myid .eq. 0 ) then
            print *, ' '
            print *, ' Error with _naupd, info = ', info
            print *, ' Check the documentation of _naupd'
            print *, ' '
         endif

      else 
!
!        %-------------------------------------------%
!        | No fatal errors occurred.                 |
!        | Post-Process using ZNEUPD.                |
!        |                                           |
!        | Computed eigenvalues may be extracted.    |
!        |                                           |
!        | Eigenvectors may also be computed now if  |
!        | desired.  (indicated by rvec = .true.)    |
!        %-------------------------------------------%
!
         rvec = .true.

         call pzneupd (comm, rvec, 'A', select, d, v, ldv, sigma, &
     &        workev, bmat, nloc, which, nev, tol, resid, ncv,    &
     &        v, ldv, iparam, ipntr, workd, workl, lworkl,        &
     &        rwork, ierr)
!
!       %----------------------------------------------%
!       | Eigenvalues are returned in the one          |
!       | dimensional array D.  The corresponding      |
!       | eigenvectors are returned in the first NCONV |
!       | (=IPARAM(5)) columns of the two dimensional  | 
!       | array V if requested.  Otherwise, an         |
!       | orthogonal basis for the invariant subspace  |
!       | corresponding to the eigenvalues in D is     |
!       | returned in V.                               |
!       %----------------------------------------------%
!
         if ( ierr .ne. 0) then
!
!          %------------------------------------%
!          | Error condition:                   |
!          | Check the documentation of ZNEUPD. |
!          %------------------------------------%
!
            if ( myid .eq. 0 ) then
                print *, ' '
                print *, ' Error with _neupd, info = ', ierr
                print *, ' Check the documentation of _neupd. '
                print *, ' '
            endif

         else

             nconv = iparam(5)
             do 20 j=1, nconv
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
                call pzndrv1av(comm, nloc, nx, mv_buf, v(1,j), ax)
                call zaxpy(nloc, -d(j), v(1,j), 1, ax, 1)
                rd(j,1) = dble(d(j))
                rd(j,2) = dimag(d(j))
                rd(j,3) = pdznorm2(comm, nloc, ax, 1)

 20          continue
!
!           %-----------------------------%
!           | Display computed residuals. |
!           %-----------------------------%
!
             call pdmout(comm, 6, nconv, 3, rd, maxncv, -6, 'Ritz values (Real, Imag) and direct residuals')
          end if
!
!        %-------------------------------------------%
!        | Print additional convergence information. |
!        %-------------------------------------------%
!
         if (myid .eq. 0)then
         if ( info .eq. 1) then
             print *, ' '
             print *, ' Maximum number of iterations reached.'
             print *, ' '
         else if ( info .eq. 3) then
             print *, ' ' 
             print *, ' No shifts could be applied during implicit Arnoldi update, try increasing NCV.'
             print *, ' '
         end if      

         print *, ' '
         print *, '_NDRV1'
         print *, '====== '
         print *, ' '
         print *, ' Size of the matrix is ', n
         print *, ' The number of processors is ', nprocs
         print *, ' The number of Ritz values requested is ', nev
         print *, ' The number of Arnoldi vectors generated (NCV) is ', ncv
         print *, ' What portion of the spectrum: ', which
         print *, ' The number of converged Ritz values is ', nconv 
         print *, ' The number of Implicit Arnoldi update iterations taken is ', iparam(3)
         print *, ' The number of OP*x is ', iparam(9)
         print *, ' The convergence criterion is ', tol
         print *, ' '

         endif
      end if
!
!     %----------------------------%
!     | Done with program pzndrv1. |
!     %----------------------------%
!
 9000 continue
!
!
!     %-------------------------%
!     | Release resources MPI |
!     %-------------------------%
!
!      call MPI_FINALIZE(rc)
!   
    END SUBROUTINE pzndrv1
    
      subroutine pzndrv1av (comm, nloc, nx, mv_buf, v, w )
!
!     .. MPI Declarations ...
      include           'mpif.h'
      integer           comm, nprocs, myid, ierr,status(MPI_STATUS_SIZE)

      integer           nloc, nx, np, j, lo, next, prev
      Complex*16        v(nloc), w(nloc), mv_buf(nx), one
      parameter         (one = (1.0, 0.0))
      external          zaxpy, tv

      call MPI_COMM_RANK( comm, myid, ierr )
      call MPI_COMM_SIZE( comm, nprocs, ierr )

      np = nloc/nx
      call pzndrv1tv(nx,v(1),w(1))
      call zaxpy(nx, -one, v(nx+1), 1, w(1), 1)

      do 10 j = 2, np-1
         lo = (j-1)*nx
         call pzndrv1tv(nx, v(lo+1), w(lo+1))
         call zaxpy(nx, -one, v(lo-nx+1), 1, w(lo+1), 1)
         call zaxpy(nx, -one, v(lo+nx+1), 1, w(lo+1), 1)
  10  continue 

      lo = (np-1)*nx
      call pzndrv1tv(nx, v(lo+1), w(lo+1))
      call zaxpy(nx, -one, v(lo-nx+1), 1, w(lo+1), 1)

      next = myid + 1
      prev = myid - 1
      if ( myid .lt. nprocs-1 ) then
         call mpi_send( v((np-1)*nx+1), nx, MPI_DOUBLE_COMPLEX, next, myid+1, comm, ierr )
      endif
      if ( myid .gt. 0 ) then
         call mpi_recv( mv_buf, nx, MPI_DOUBLE_COMPLEX, prev, myid, comm, status, ierr )
         call zaxpy( nx, -one, mv_buf, 1, w(1), 1 )
      endif

      if ( myid .gt. 0 ) then
         call mpi_send( v(1), nx, MPI_DOUBLE_COMPLEX, prev, myid-1, comm, ierr )
      endif
      if ( myid .lt. nprocs-1 ) then
         call mpi_recv( mv_buf, nx, MPI_DOUBLE_COMPLEX, next, myid, comm, status, ierr )
         call zaxpy( nx, -one, mv_buf, 1, w(lo+1), 1 )
      endif

      return
    end
    
        subroutine pzndrv1tv (nx, x, y)

      integer           nx, j 
      Complex*16        x(nx), y(nx), h, dd, dl, du

      Complex*16        one, rho
      parameter         (one = (1.0, 0.0), rho = (100.0, 0.0))
!
!     Compute the matrix vector multiplication y<---T*x
!     where T is a nx by nx tridiagonal matrix with DD on the 
!     diagonal, DL on the subdiagonal, and DU on the superdiagonal
!     
      h   = one / dcmplx(nx+1)
      dd  = (4.0, 0.0)
      dl  = -one - (0.5, 0.0)*rho*h
      du  = -one + (0.5, 0.0)*rho*h
 
      y(1) =  dd*x(1) + du*x(2)
      do 10 j = 2,nx-1
         y(j) = dl*x(j-1) + dd*x(j) + du*x(j+1) 
 10   continue 
      y(nx) =  dl*x(nx-1) + dd*x(nx) 
      return
      end  