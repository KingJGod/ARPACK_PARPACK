    !\BeginDoc
    !
    !\Name: zneigh
    !
    !\Description:
    !  Compute the eigenvalues of the current upper Hessenberg matrix
    !  and the corresponding Ritz estimates given the current residual norm.
    !
    !\Usage:
    !  call zneigh
    !     ( RNORM, N, H, LDH, RITZ, BOUNDS, Q, LDQ, WORKL, RWORK, IERR )
    !
    !\Arguments
    !  RNORM   Double precision scalar.  (INPUT)
    !          Residual norm corresponding to the current upper Hessenberg
    !          matrix H.
    !
    !  N       Integer.  (INPUT)
    !          Size of the matrix H.
    !
    !  H       Complex*16 N by N array.  (INPUT)
    !          H contains the current upper Hessenberg matrix.
    !
    !  LDH     Integer.  (INPUT)
    !          Leading dimension of H exactly as declared in the calling
    !          program.
    !
    !  RITZ    Complex*16 array of length N.  (OUTPUT)
    !          On output, RITZ(1:N) contains the eigenvalues of H.
    !
    !  BOUNDS  Complex*16 array of length N.  (OUTPUT)
    !          On output, BOUNDS contains the Ritz estimates associated with
    !          the eigenvalues held in RITZ.  This is equal to RNORM
    !          times the last components of the eigenvectors corresponding
    !          to the eigenvalues in RITZ.
    !
    !  Q       Complex*16 N by N array.  (WORKSPACE)
    !          Workspace needed to store the eigenvectors of H.
    !
    !  LDQ     Integer.  (INPUT)
    !          Leading dimension of Q exactly as declared in the calling
    !          program.
    !
    !  WORKL   Complex*16 work array of length N**2 + 3*N.  (WORKSPACE)
    !          Private (replicated) array on each PE or array allocated on
    !          the front end.  This is needed to keep the full Schur form
    !          of H and also in the calculation of the eigenvectors of H.
    !
    !  RWORK   Double precision  work array of length N (WORKSPACE)
    !          Private (replicated) array on each PE or array allocated on
    !          the front end.
    !
    !  IERR    Integer.  (OUTPUT)
    !          Error exit flag from zlahqr or ztrevc.
    !
    !\EndDoc
    !
    !-----------------------------------------------------------------------
    !
    !\BeginLib
    !
    !\Local variables:
    !     xxxxxx  Complex*16
    !
    !\Routines called:
    !     ivout   ARPACK utility routine that prints integers.
    !     arscnd  ARPACK utility routine for timing.
    !     zmout   ARPACK utility routine that prints matrices
    !     zvout   ARPACK utility routine that prints vectors.
    !     dvout   ARPACK utility routine that prints vectors.
    !     zlacpy  LAPACK matrix copy routine.
    !     zlahqr  LAPACK routine to compute the Schur form of an
    !             upper Hessenberg matrix.
    !     zlaset  LAPACK matrix initialization routine.
    !     ztrevc  LAPACK routine to compute the eigenvectors of a matrix
    !             in upper triangular form
    !     zcopy   Level 1 BLAS that copies one vector to another.
    !     zdscal  Level 1 BLAS that scales a complex vector by a real number.
    !     dznrm2  Level 1 BLAS that computes the norm of a vector.
    !
    !
    !\Author
    !     Danny Sorensen               Phuong Vu
    !     Richard Lehoucq              CRPC / Rice University
    !     Dept. of Computational &     Houston, Texas
    !     Applied Mathematics
    !     Rice University
    !     Houston, Texas
    !
    !\SCCS Information: @(#)
    ! FILE: neigh.F   SID: 2.2   DATE OF SID: 4/20/96   RELEASE: 2
    !
    !\Remarks
    !     None
    !
    !\EndLib
    !
    !-----------------------------------------------------------------------
    !
    SUBROUTINE zneigh(rnorm, n, h, ldh, ritz, bounds, &
        &            q, ldq, workl, rwork, ierr)

    USE DEBUG_MODULE
    USE STAT_MODULE
    IMPLICIT NONE
    !
    !     %------------------%
    !     | Scalar Arguments |
    !     %------------------%
    !
    REAL( kind = 8 ),INTENT(IN)::RNORM
    INTEGER( kind = 4 ),INTENT(IN)::N
    INTEGER( kind = 4 ),INTENT(IN)::LDH
    INTEGER( kind = 4 ),INTENT(IN)::LDQ
    INTEGER( kind = 4 ),INTENT(OUT)::IERR
    !
    !     %-----------------%
    !     | Array Arguments |
    !     %-----------------%
    !
    COMPLEX( kind = 8 ), INTENT(OUT)::BOUNDS(N)
    COMPLEX( kind = 8 ), INTENT(IN)::H(LDH,N)
    COMPLEX( kind = 8 ), INTENT(INOUT)::Q(LDQ,N)
    COMPLEX( kind = 8 ), INTENT(OUT)::RITZ(N)
    COMPLEX( kind = 8 ), INTENT(INOUT)::WORKL(N*(N+3))
    REAL( kind = 8 ), INTENT(INOUT)::RWORK(N)
    !
    !     %------------%
    !     | Parameters |
    !     %------------%
    !
    COMPLEX( kind = 8 ), PARAMETER::  one = (1.d0,0.d0)
    COMPLEX( kind = 8 ), PARAMETER:: zero = (0.d0,0.d0)
    REAL( kind = 8 ),    PARAMETER:: rone = 1.d0
    !
    !     %------------------------%
    !     | Local Scalars & Arrays |
    !     %------------------------%
    !
    LOGICAL::select(1)
    INTEGER( kind = 4 )::j,  msglvl
    COMPLEX( kind = 8 )::vl(1)
    REAL( kind = 8 )::temp
    !
    !     %----------------------%
    !     | External Subroutines |
    !     %----------------------%
    !      In BLAS & LAPACK
    EXTERNAL zlacpy, zlahqr, ztrevc, zcopy
    EXTERNAL zdscal, zmout, zvout, arscnd
    !
    !     %--------------------%
    !     | External functions |
    !     %--------------------%
    !
    REAL( kind = 8 ) dznrm2
    EXTERNAL   dznrm2
    !
    !
    !     %-----------------------%
    !     | Executable Statements |
    !     %-----------------------%
    !
    !     %-------------------------------%
    !     | Initialize timing statistics  |
    !     | & message level for debugging |
    !     %-------------------------------%
    !
    call arscnd (t0)
    msglvl = mceigh
    !
    if (msglvl .gt. 2) then
        call zmout (logfil, n, n, h, ldh, ndigit, '_neigh: Entering upper Hessenberg matrix H ')
    end if
    !
    !     %----------------------------------------------------------%
    !     | 1. Compute the eigenvalues, the last components of the   |
    !     |    corresponding Schur vectors and the full Schur form T |
    !     |    of the current upper Hessenberg matrix H.             |
    !     |    zlahqr returns the full Schur form of H               |
    !     |    in WORKL(1:N**2), and the Schur vectors in q.         |
    !     %----------------------------------------------------------%
    !
    call zlacpy ('All', n, n, h, ldh, workl, n)
    call zlaset ('All', n, n, zero, one, q, ldq)
    call zlahqr (.true., .true., n, 1, n, workl, ldh, ritz, 1, n, q, ldq, ierr)
    if (ierr .ne. 0) then
        go to 9000
    endif
    call zcopy (n, q(n-1,1), ldq, bounds, 1)
    if (msglvl .gt. 1) then
        call zvout (logfil, n, bounds, ndigit, '_neigh: last row of the Schur matrix for H')
    end if
    !
    !     %----------------------------------------------------------%
    !     | 2. Compute the eigenvectors of the full Schur form T and |
    !     |    apply the Schur vectors to get the corresponding      |
    !     |    eigenvectors.                                         |
    !     %----------------------------------------------------------%
    !
    call ztrevc ('Right', 'Back', select, n, workl, n, vl, n, q, ldq, n, n, workl(n*n+1), rwork, ierr)
    !
    if (ierr .ne. 0) go to 9000
    !
    !     %------------------------------------------------%
    !     | Scale the returning eigenvectors so that their |
    !     | Euclidean norms are all one. LAPACK subroutine |
    !     | ztrevc returns each eigenvector normalized so  |
    !     | that the element of largest magnitude has      |
    !     | magnitude 1; here the magnitude of a complex   |
    !     | number (x,y) is taken to be |x| + |y|.         |
    !     %------------------------------------------------%
    !
    do 10 j=1, n
        temp = dznrm2( n, q(1,j), 1 )
        call zdscal ( n, rone / temp, q(1,j), 1 )
10  continue
    !
    if (msglvl .gt. 1) then
        call zcopy(n, q(n,1), ldq, workl, 1)
        call zvout(logfil, n, workl, ndigit, '_neigh: Last row of the eigenvector matrix for H')
    end if
    !
    !     %----------------------------%
    !     | Compute the Ritz estimates |
    !     %----------------------------%
    !
    call zcopy(n, q(n,1), n, bounds, 1)
    call zdscal(n, rnorm, bounds, 1)
    !
    if (msglvl .gt. 2) then
        call zvout(logfil, n, ritz, ndigit, '_neigh: The eigenvalues of H')
        call zvout(logfil, n, bounds, ndigit, '_neigh: Ritz estimates for the eigenvalues of H')
    end if
    !
    call arscnd(t1)
    tceigh = tceigh + (t1 - t0)
    !
9000 continue
    return
    !
    !     %---------------%
    !     | End of zneigh |
    !     %---------------%
    !
    END SUBROUTINE zneigh