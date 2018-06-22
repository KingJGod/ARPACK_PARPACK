    !\BeginDoc
    !
    !\Name: zgetv0
    !
    !\Description:
    !  Generate a random initial residual vector for the Arnoldi process.
    !  Force the residual vector to be in the range of the operator OP.
    !
    !\Usage:
    !  call zgetv0
    !     ( IDO, BMAT, ITRY, INITV, N, J, V, LDV, RESID, RNORM,
    !       IPNTR, WORKD, IERR )
    !
    !\Arguments
    !  IDO     Integer.  (INPUT/OUTPUT)
    !          Reverse communication flag.  IDO must be zero on the first
    !          call to zgetv0.
    !          -------------------------------------------------------------
    !          IDO =  0: first call to the reverse communication interface
    !          IDO = -1: compute  Y = OP * X  where
    !                    IPNTR(1) is the pointer into WORKD for X,
    !                    IPNTR(2) is the pointer into WORKD for Y.
    !                    This is for the initialization phase to force the
    !                    starting vector into the range of OP.
    !          IDO =  2: compute  Y = B * X  where
    !                    IPNTR(1) is the pointer into WORKD for X,
    !                    IPNTR(2) is the pointer into WORKD for Y.
    !          IDO = 99: done
    !          -------------------------------------------------------------
    !
    !  BMAT    Character*1.  (INPUT)
    !          BMAT specifies the type of the matrix B in the (generalized)
    !          eigenvalue problem A*x = lambda*B*x.
    !          B = 'I' -> standard eigenvalue problem A*x = lambda*x
    !          B = 'G' -> generalized eigenvalue problem A*x = lambda*B*x
    !
    !  ITRY    Integer.  (INPUT)
    !          ITRY counts the number of times that zgetv0 is called.
    !          It should be set to 1 on the initial call to zgetv0.
    !
    !  INITV   Logical variable.  (INPUT)
    !          .TRUE.  => the initial residual vector is given in RESID.
    !          .FALSE. => generate a random initial residual vector.
    !
    !  N       Integer.  (INPUT)
    !          Dimension of the problem.
    !
    !  J       Integer.  (INPUT)
    !          Index of the residual vector to be generated, with respect to
    !          the Arnoldi process.  J > 1 in case of a "restart".
    !
    !  V       Complex*16 N by J array.  (INPUT)
    !          The first J-1 columns of V contain the current Arnoldi basis
    !          if this is a "restart".
    !
    !  LDV     Integer.  (INPUT)
    !          Leading dimension of V exactly as declared in the calling
    !          program.
    !
    !  RESID   Complex*16 array of length N.  (INPUT/OUTPUT)
    !          Initial residual vector to be generated.  If RESID is
    !          provided, force RESID into the range of the operator OP.
    !
    !  RNORM   Double precision scalar.  (OUTPUT)
    !          B-norm of the generated residual.
    !
    !  IPNTR   Integer array of length 3.  (OUTPUT)
    !
    !  WORKD   Complex*16 work array of length 2*N.  (REVERSE COMMUNICATION).
    !          On exit, WORK(1:N) = B*RESID to be used in SSAITR.
    !
    !  IERR    Integer.  (OUTPUT)
    !          =  0: Normal exit.
    !          = -1: Cannot generate a nontrivial restarted residual vector
    !                in the range of the operator OP.
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
    !\References:
    !  1. D.C. Sorensen, "Implicit Application of Polynomial Filters in
    !     a k-Step Arnoldi Method", SIAM J. Matr. Anal. Apps., 13 (1992),
    !     pp 357-385.
    !
    !\Routines called:
    !     arscnd  ARPACK utility routine for timing.
    !     zvout   ARPACK utility routine that prints vectors.
    !     zlarnv  LAPACK routine for generating a random vector.
    !     zgemv   Level 2 BLAS routine for matrix vector multiplication.
    !     zcopy   Level 1 BLAS that copies one vector to another.
    !     zdotc   Level 1 BLAS that computes the scalar product of two vectors.
    !     dznrm2  Level 1 BLAS that computes the norm of a vector.
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
    ! FILE: getv0.F   SID: 2.3   DATE OF SID: 08/27/96   RELEASE: 2
    !
    !\EndLib
    !
    !-----------------------------------------------------------------------
    !

    SUBROUTINE zgetv0(ido ,bmat, itry, initv, n, j, v, ldv, resid, rnorm, &
        & ipntr, workd, ierr)
    USE DEBUG_MODULE
    USE STAT_MODULE
    IMPLICIT NONE
    !
    !     %------------------%
    !     | Scalar Arguments |
    !     %------------------%
    !
    CHARACTER(len = 1 ), INTENT(IN)::BMAT
    INTEGER( kind = 4 ), INTENT(IN)::ITRY
    LOGICAL,             INTENT(IN)::initv
    INTEGER( kind = 4 ), INTENT(IN)::N
    INTEGER( kind = 4 ), INTENT(IN)::J
    INTEGER( kind = 4 ), INTENT(IN)::LDV
    INTEGER( kind = 4 ), INTENT(INOUT)::IDO
    INTEGER( kind = 4 ), INTENT(OUT)::IERR
    REAL( kind = 8 ),    INTENT(OUT)::RNORM
    !
    !     %-----------------%
    !     | Array Arguments |
    !     %-----------------%
    !
    COMPLEX( kind = 8 ), INTENT(IN)::V(LDV,J)
    COMPLEX( kind = 8 ), INTENT(INOUT)::RESID(N)
    COMPLEX( kind = 8 ), INTENT(INOUT)::WORKD(2*N)
    INTEGER( kind = 4 ), INTENT(OUT)::ipntr(3)
    !
    !     %------------%
    !     | Parameters |
    !     %------------%
    !
    COMPLEX( kind = 8 ), PARAMETER:: one = (1.d0,0.d0)
    COMPLEX( kind = 8 ), PARAMETER::zero = (0.d0,0.d0)
    REAL( kind = 8 ),    PARAMETER::rzero= 0.d0
    !
    !     %------------------------%
    !     | Local Scalars & Arrays |
    !     %------------------------%
    !
    LOGICAL::    first, inits, orth
    INTEGER( kind = 4 )::idist, iseed(4), iter, msglvl, jj
    REAL( kind = 8 ):: rnorm0
    COMPLEX( kind = 8):: cnorm
    SAVE first, iseed, inits, iter, msglvl, orth, rnorm0
    !
    !     %----------------------%
    !     | External Subroutines |
    !     %----------------------%
    !      In BLAS & LAPACK
    external zcopy, zgemv, zlarnv, zvout, arscnd
    !
    !     %--------------------%
    !     | External Functions |
    !     %--------------------%
    !
    REAL( kind = 8 ):: dznrm2, dlapy2
    COMPLEX( kind = 8 ):: zdotc
    EXTERNAL zdotc,dznrm2,dlapy2
    !
    !     %-----------------%
    !     | Data Statements |
    !     %-----------------%
    !
    DATA inits /.true./

    !
    !     %-----------------------%
    !     | Executable Statements |
    !     %-----------------------%
    !
    !
    !     %-----------------------------------%
    !     | Initialize the seed of the LAPACK |
    !     | random number generator           |
    !     %-----------------------------------%
    !
    if (inits) then
        iseed(1) = 1
        iseed(2) = 3
        iseed(3) = 5
        iseed(4) = 7
        inits = .false.
    end if
    !
    if (ido .eq.  0) then
        !
        !        %-------------------------------%
        !        | Initialize timing statistics  |
        !        | & message level for debugging |
        !        %-------------------------------%
        !
        call arscnd (t0)
        msglvl = mgetv0
        !
        ierr   = 0
        iter   = 0
        first  = .FALSE.
        orth   = .FALSE.
        !
        !        %-----------------------------------------------------%
        !        | Possibly generate a random starting vector in RESID |
        !        | Use a LAPACK random number generator used by the    |
        !        | matrix generation routines.                         |
        !        |    idist = 1: uniform (0,1)  distribution;          |
        !        |    idist = 2: uniform (-1,1) distribution;          |
        !        |    idist = 3: normal  (0,1)  distribution;          |
        !        %-----------------------------------------------------%
        !
        if (.not.initv) then
            idist = 2
            call zlarnv (idist, iseed, n, resid)
        end if
        !
        !        %----------------------------------------------------------%
        !        | Force the starting vector into the range of OP to handle |
        !        | the generalized problem when B is possibly (singular).   |
        !        %----------------------------------------------------------%
        !
        call arscnd (t2)
        if (bmat .eq. 'G') then
            nopx = nopx + 1
            ipntr(1) = 1
            ipntr(2) = n + 1
            call zcopy (n, resid, 1, workd, 1)
            ido = -1
            go to 9000
        end if
    end if
    !
    !     %----------------------------------------%
    !     | Back from computing B*(initial-vector) |
    !     %----------------------------------------%
    !
    if (first) go to 20
    !
    !     %-----------------------------------------------%
    !     | Back from computing B*(orthogonalized-vector) |
    !     %-----------------------------------------------%
    !
    if (orth)  go to 40
    !
    call arscnd (t3)
    tmvopx = tmvopx + (t3 - t2)
    !
    !     %------------------------------------------------------%
    !     | Starting vector is now in the range of OP; r = OP*r; |
    !     | Compute B-norm of starting vector.                   |
    !     %------------------------------------------------------%
    !
    call arscnd (t2)
    first = .TRUE.
    if (bmat .eq. 'G') then
        nbx = nbx + 1
        call zcopy (n, workd(n+1), 1, resid, 1)
        ipntr(1) = n + 1
        ipntr(2) = 1
        ido = 2
        go to 9000
    else if (bmat .eq. 'I') then
        call zcopy (n, resid, 1, workd, 1)
    end if
    !
20  continue
    !
    if (bmat .eq. 'G') then
        call arscnd (t3)
        tmvbx = tmvbx + (t3 - t2)
    end if
    !
    first = .FALSE.
    if (bmat .eq. 'G') then
        cnorm  = zdotc (n, resid, 1, workd, 1)
        rnorm0 = sqrt(dlapy2(dble(cnorm),dimag(cnorm)))
    else if (bmat .eq. 'I') then
        rnorm0 = dznrm2(n, resid, 1)
    end if
    rnorm  = rnorm0
    !
    !     %---------------------------------------------%
    !     | Exit if this is the very first Arnoldi step |
    !     %---------------------------------------------%
    !
    if (j .eq. 1) go to 50
    !
    !     %----------------------------------------------------------------
    !     | Otherwise need to B-orthogonalize the starting vector against |
    !     | the current Arnoldi basis using Gram-Schmidt with iter. ref.  |
    !     | This is the case where an invariant subspace is encountered   |
    !     | in the middle of the Arnoldi factorization.                   |
    !     |                                                               |
    !     |       s = V^{T}*B*r;   r = r - V*s;                           |
    !     |                                                               |
    !     | Stopping criteria used for iter. ref. is discussed in         |
    !     | Parlett's book, page 107 and in Gragg & Reichel TOMS paper.   |
    !     %---------------------------------------------------------------%
    !
    orth = .TRUE.
30  continue
    !
    call zgemv ('C', n, j-1, one, v, ldv, workd, 1, zero, workd(n+1), 1)
    call zgemv ('N', n, j-1, -one, v, ldv, workd(n+1), 1, one, resid, 1)
    !
    !     %----------------------------------------------------------%
    !     | Compute the B-norm of the orthogonalized starting vector |
    !     %----------------------------------------------------------%
    !
    call arscnd (t2)
    if (bmat .eq. 'G') then
        nbx = nbx + 1
        call zcopy (n, resid, 1, workd(n+1), 1)
        ipntr(1) = n + 1
        ipntr(2) = 1
        ido = 2
        go to 9000
    else if (bmat .eq. 'I') then
        call zcopy (n, resid, 1, workd, 1)
    end if
    !
40  continue
    !
    if (bmat .eq. 'G') then
        call arscnd (t3)
        tmvbx = tmvbx + (t3 - t2)
    end if
    !
    if (bmat .eq. 'G') then
        cnorm = zdotc (n, resid, 1, workd, 1)
        rnorm = sqrt(dlapy2(dble(cnorm),dimag(cnorm)))
    else if (bmat .eq. 'I') then
        rnorm = dznrm2(n, resid, 1)
    end if
    !
    !     %--------------------------------------%
    !     | Check for further orthogonalization. |
    !     %--------------------------------------%
    !
    !if (msglvl .gt. 2) then
    !    call dvout (logfil, 1, rnorm0, ndigit,'_getv0: re-orthonalization ; rnorm0 is')
    !    call dvout (logfil, 1, rnorm, ndigit,'_getv0: re-orthonalization ; rnorm is')
    !end if
    !
    if (rnorm .gt. 0.717*rnorm0) go to 50
    !
    iter = iter + 1
    if (iter .le. 1) then
        !
        !        %-----------------------------------%
        !        | Perform iterative refinement step |
        !        %-----------------------------------%
        !
        rnorm0 = rnorm
        go to 30
    else
        !
        !        %------------------------------------%
        !        | Iterative refinement step "failed" |
        !        %------------------------------------%
        !
        do 45 jj = 1, n
            resid(jj) = zero
45      continue
        rnorm = rzero
        ierr = -1
    end if
    !
50  continue
    !
    !if (msglvl .gt. 0) then
    !    call dvout (logfil, 1, rnorm, ndigit,'_getv0: B-norm of initial / restarted starting vector')
    !end if
    if (msglvl .gt. 2) then
        call zvout (logfil, n, resid, ndigit,'_getv0: initial / restarted starting vector')
    end if
    ido = 99
    !
    call arscnd (t1)
    tgetv0 = tgetv0 + (t1 - t0)
    !
9000 continue
    return
    !
    !     %---------------%
    !     | End of zgetv0 |
    !     %---------------%


    END SUBROUTINE zgetv0