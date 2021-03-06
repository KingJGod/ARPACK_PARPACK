!\BeginDoc
!
!\Name: pdznorm2
!
! Message Passing Layer: MPI
!
!\Description:
!
!\Usage:
!  call pdznorm2 ( COMM, N, X, INC )
!
!\Arguments
!  COMM    MPI Communicator for the processor grid.  (INPUT)
!
!\SCCS Information:
! FILE: norm2.F   SID: 1.2   DATE OF SID: 3/6/96
!
!-----------------------------------------------------------------------
!
      Double precision function pdznorm2( comm, n, x, inc )
      
      include   'mpif.h'
!
!     %---------------%
!     | MPI Variables |
!     %---------------%
!
      INTEGER( kind = 4 )::comm, ierr
!
!     %------------------%
!     | Scalar Arguments |
!     %------------------%
!
      INTEGER( kind = 4 )::n, inc
!
!     %-----------------%
!     | Array Arguments |
!     %-----------------%
!
      COMPLEX( kind = 8 )::x(n)
!
!     %---------------%
!     | Local Scalars |
!     %---------------%
!
      Double precision max, buf, zero
      parameter    ( zero = 0.0 )
!
!     %---------------------%
!     | Intrinsic Functions |
!     %---------------------%
!
      intrinsic    abs, sqrt
!
!     %--------------------%
!     | External Functions |
!     %--------------------%
!
      Double precision  dznrm2
      External     dznrm2
!
!     %-----------------------%
!     | Executable Statements |
!     %-----------------------%
!
      pdznorm2 = dznrm2( n, x, inc)
!
      buf = pdznorm2
      call MPI_ALLREDUCE( buf, max, 1, MPI_DOUBLE_PRECISION, MPI_MAX, comm, ierr )
      if ( max .eq. zero ) then
         pdznorm2 = zero
      else
         buf = (pdznorm2/max)**2.0
         call MPI_ALLREDUCE( buf, pdznorm2, 1, MPI_DOUBLE_PRECISION, MPI_SUM, comm, ierr )
         pdznorm2 = max * sqrt(abs(pdznorm2))
      endif
!
!     %-----------------%
!     | End of pdznorm2 |
!     %-----------------%
!
      return
      end function 
