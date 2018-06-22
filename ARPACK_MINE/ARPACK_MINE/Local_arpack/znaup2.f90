    !
    !\Usage:
    !  call znaup2
    !     ( IDO, BMAT, N, WHICH, NEV, NP, TOL, RESID, MODE, IUPD,
    !       ISHIFT, MXITER, V, LDV, H, LDH, RITZ, BOUNDS,
    !       Q, LDQ, WORKL, IPNTR, WORKD, RWORK, INFO )
    !
    !\Arguments
    !
    !  IDO, BMAT, N, WHICH, NEV, TOL, RESID: same as defined in znaupd .
    !  MODE, ISHIFT, MXITER: see the definition of IPARAM in znaupd .
    !
    !  NP      Integer.  (INPUT/OUTPUT)
    !          Contains the number of implicit shifts to apply during
    !          each Arnoldi iteration.
    !          If ISHIFT=1, NP is adjusted dynamically at each iteration
    !          to accelerate convergence and prevent stagnation.
    !          This is also roughly equal to the number of matrix-vector
    !          products (involving the operator OP) per Arnoldi iteration.
    !          The logic for adjusting is contained within the current
    !          subroutine.
    !          If ISHIFT=0, NP is the number of shifts the user needs
    !          to provide via reverse comunication. 0 < NP < NCV-NEV.
    !          NP may be less than NCV-NEV since a leading block of the current
    !          upper Hessenberg matrix has split off and contains "unwanted"
    !          Ritz values.
    !          Upon termination of the IRA iteration, NP contains the number
    !          of "converged" wanted Ritz values.
    !
    !  IUPD    Integer.  (INPUT)
    !          IUPD .EQ. 0: use explicit restart instead implicit update.
    !          IUPD .NE. 0: use implicit update.
    !
    !  V       Complex*16  N by (NEV+NP) array.  (INPUT/OUTPUT)
    !          The Arnoldi basis vectors are returned in the first NEV
    !          columns of V.
    !
    !  LDV     Integer.  (INPUT)
    !          Leading dimension of V exactly as declared in the calling
    !          program.
    !
    !  H       Complex*16  (NEV+NP) by (NEV+NP) array.  (OUTPUT)
    !          H is used to store the generated upper Hessenberg matrix
    !
    !  LDH     Integer.  (INPUT)
    !          Leading dimension of H exactly as declared in the calling
    !          program.
    !
    !  RITZ    Complex*16  array of length NEV+NP.  (OUTPUT)
    !          RITZ(1:NEV)  contains the computed Ritz values of OP.
    !
    !  BOUNDS  Complex*16  array of length NEV+NP.  (OUTPUT)
    !          BOUNDS(1:NEV) contain the error bounds corresponding to
    !          the computed Ritz values.
    !
    !  Q       Complex*16  (NEV+NP) by (NEV+NP) array.  (WORKSPACE)
    !          Private (replicated) work array used to accumulate the
    !          rotation in the shift application step.
    !
    !  LDQ     Integer.  (INPUT)
    !          Leading dimension of Q exactly as declared in the calling
    !          program.
    !
    !  WORKL   Complex*16  work array of length at least
    !          (NEV+NP)**2 + 3*(NEV+NP).  (WORKSPACE)
    !          Private (replicated) array on each PE or array allocated on
    !          the front end.  It is used in shifts calculation, shifts
    !          application and convergence checking.
    !
    !
    !  IPNTR   Integer array of length 3.  (OUTPUT)
    !          Pointer to mark the starting locations in the WORKD for
    !          vectors used by the Arnoldi iteration.
    !          -------------------------------------------------------------
    !          IPNTR(1): pointer to the current operand vector X.
    !          IPNTR(2): pointer to the current result vector Y.
    !          IPNTR(3): pointer to the vector B * X when used in the
    !                    shift-and-invert mode.  X is the current operand.
    !          -------------------------------------------------------------
    !
    !  WORKD   Complex*16  work array of length 3*N.  (WORKSPACE)
    !          Distributed array to be used in the basic Arnoldi iteration
    !          for reverse communication.  The user should not use WORKD
    !          as temporary workspace during the iteration !!!!!!!!!!
    !          See Data Distribution Note in ZNAUPD .
    !
    !  RWORK   Double precision    work array of length  NEV+NP ( WORKSPACE)
    !          Private (replicated) array on each PE or array allocated on
    !          the front end.
    !
    !  INFO    Integer.  (INPUT/OUTPUT)
    !          If INFO .EQ. 0, a randomly initial residual vector is used.
    !          If INFO .NE. 0, RESID contains the initial residual vector,
    !                          possibly from a previous run.
    !          Error flag on output.
    !          =     0: Normal return.
    !          =     1: Maximum number of iterations taken.
    !                   All possible eigenvalues of OP has been found.
    !                   NP returns the number of converged Ritz values.
    !          =     2: No shifts could be applied.
    !          =    -8: Error return from LAPACK eigenvalue calculation;
    !                   This should never happen.
    !          =    -9: Starting vector is zero.
    !          = -9999: Could not build an Arnoldi factorization.
    !                   Size that was built in returned in NP.
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
    !  2. R.B. Lehoucq, "Analysis and Implementation of an Implicitly
    !     Restarted Arnoldi Iteration", Rice University Technical Report
    !     TR95-13, Department of Computational and Applied Mathematics.
    !
    !\Routines called:
    !     zgetv0   ARPACK initial vector generation routine.
    !     znaitr   ARPACK Arnoldi factorization routine.
    !     znapps   ARPACK application of implicit shifts routine.
    !     zneigh   ARPACK compute Ritz values and error bounds routine.
    !     zngets   ARPACK reorder Ritz values and error bounds routine.
    !     zsortc   ARPACK sorting routine.
    !     ivout   ARPACK utility routine that prints integers.
    !     arscnd  ARPACK utility routine for timing.
    !     zmout    ARPACK utility routine that prints matrices
    !     zvout    ARPACK utility routine that prints vectors.
    !     dvout    ARPACK utility routine that prints vectors.
    !     dlamch   LAPACK routine that determines machine constants.
    !     dlapy2   LAPACK routine to compute sqrt(x**2+y**2) carefully.
    !     zcopy    Level 1 BLAS that copies one vector to another .
    !     zdotc    Level 1 BLAS that computes the scalar product of two vectors.
    !     zswap    Level 1 BLAS that swaps two vectors.
    !     dznrm2   Level 1 BLAS that computes the norm of a vector.
    !
    !\Author
    !     Danny Sorensen               Phuong Vu
    !     Richard Lehoucq              CRPC / Rice Universitya
    !     Chao Yang                    Houston, Texas
    !     Dept. of Computational &
    !     Applied Mathematics
    !     Rice University
    !     Houston, Texas
    !
    !\SCCS Information: @(#)
    ! FILE: naup2.F   SID: 2.6   DATE OF SID: 06/01/00   RELEASE: 2
    !
    !\Remarks
    !     1. None
    !
    !\EndLib
    !
    !-----------------------------------------------------------------------
    !
    SUBROUTINE znaup2( ido, bmat, n, which, nev, np, tol, resid, mode, iupd, &
        &     ishift, mxiter, v, ldv, h, ldh, ritz,                 &
        &     bounds, q, ldq, workl,                                &
        &     ipntr, workd, rwork, info )
    USE DEBUG_MODULE
    USE STAT_MODULE
    IMPLICIT NONE
    !
    !     %------------------%
    !     | Scalar Arguments |
    !     %------------------%
    !
    INTEGER( kind = 4 ), INTENT(INOUT)::IDO
    CHARACTER( len = 1 ),INTENT(IN)::BMAT
    INTEGER( kind = 4 ), INTENT(IN)::N
    CHARACTER( len = 2 ),INTENT(IN)::WHICH
    INTEGER( kind = 4 ), INTENT(INOUT)::NEV
    REAL( kind = 8 ), INTENT(INOUT)::TOL
    INTEGER( kind = 4 ), INTENT(INOUT)::INFO
    INTEGER( kind = 4 ), INTENT(INOUT)::MODE
    INTEGER( kind = 4 ), INTENT(INOUT)::ISHIFT
    INTEGER( kind = 4 ), INTENT(INOUT)::MXITER
    INTEGER( kind = 4 ), INTENT(INOUT)::NP
    INTEGER( kind = 4 ), INTENT(IN)::IUPD
    INTEGER( kind = 4 ), INTENT(IN)::LDV
    INTEGER( kind = 4 ), INTENT(IN)::LDH
    INTEGER( kind = 4 ), INTENT(IN)::LDQ
    !
    !     %-----------------%
    !     | Array Arguments |
    !     %-----------------%
    !
    COMPLEX( kind = 8 ), INTENT(INOUT)::RESID(N)
    COMPLEX( kind = 8 ), INTENT(INOUT)::V(LDV,NEV+NP)
    COMPLEX( kind = 8 ), INTENT(OUT)::H(LDH,NEV+NP)
    COMPLEX( kind = 8 ), INTENT(OUT)::RITZ(NEV+NP)
    COMPLEX( kind = 8 ), INTENT(OUT)::BOUNDS(NEV+NP)
    COMPLEX( kind = 8 ), INTENT(INOUT)::Q(LDQ,NEV+NP)
    COMPLEX( kind = 8 ), INTENT(INOUT)::WORKL((NEV+NP)*(NEV+NP+3))
    INTEGER( kind = 4 ), INTENT(OUT)::IPNTR(13)
    COMPLEX( kind = 8 ), INTENT(INOUT)::WORKD(3*N)
    REAL( kind = 8 ),    INTENT(INOUT)::RWORK(NEV+NP)
    !
    !     %------------%
    !     | Parameters |
    !     %------------%
    !
    COMPLEX( kind = 8 ), PARAMETER::  one = (1.d0,0.d0)
    COMPLEX( kind = 8 ), PARAMETER:: zero = (0.d0,0.d0)
    REAL( kind = 8 ),    PARAMETER::rzero = 0.d0
    !
    !     %---------------%
    !     | Local Scalars |
    !     %---------------%
    !
    LOGICAL    cnorm , getv0, initv , update, ushift
    INTEGER( kind = 4 )::ierr  , iter , kplusp, msglvl, nconv
    INTEGER( kind = 4 )::nevbef, nev0 , np0 , nptemp, i ,j
    COMPLEX( kind = 8 )::cmpnorm
    REAL( kind = 8 )::rnorm , eps23, rtemp
    CHARACTER( len = 2)::wprime
    !
    SAVE  cnorm,  getv0, initv , update, ushift
    SAVE  rnorm,  iter , kplusp, msglvl, nconv
    SAVE  nevbef, nev0 , np0   , eps23
    !
    !     %-----------------------%
    !     | Local array arguments |
    !     %-----------------------%
    !
    INTEGER( kind = 4 )::kp(3)
    !
    !     %----------------------%
    !     | External Subroutines |
    !     %----------------------%
    !      In BLAS & LAPACK
    EXTERNAL zcopy , zgetv0 , znaitr , zneigh , zngets , znapps
    EXTERNAL zsortc , zswap , zmout , zvout , ivout, arscnd
    !
    !     %--------------------%
    !     | External functions |
    !     %--------------------%
    !
    COMPLEX( kind = 8 )::zdotc
    REAL( kind = 8 )::dznrm2 , dlamch , dlapy2
    EXTERNAL zdotc , dznrm2 , dlamch , dlapy2
    !
    !     %---------------------%
    !     | Intrinsic Functions |
    !     %---------------------%
    !
    INTRINSIC  dimag , dble , min, max
    !
    !     %-----------------------%
    !     | Executable Statements |
    !     %-----------------------%
    !
    if (ido .eq. 0) then
        !
        call arscnd (t0)
        !
        msglvl = mcaup2
        !
        nev0   = nev
        np0    = np
        !
        !        %-------------------------------------%
        !        | kplusp is the bound on the largest  |
        !        |        Lanczos factorization built. |
        !        | nconv is the current number of      |
        !        |        "converged" eigenvalues.     |
        !        | iter is the counter on the current  |
        !        |      iteration step.                |
        !        %-------------------------------------%
        !
        kplusp = nev + np
        nconv  = 0
        iter   = 0
        !
        !        %---------------------------------%
        !        | Get machine dependent constant. |
        !        %---------------------------------%
        !
        eps23 = dlamch ('Epsilon-Machine')
        eps23 = eps23**(2.0D+0  / 3.0D+0 )
        !
        !        %---------------------------------------%
        !        | Set flags for computing the first NEV |
        !        | steps of the Arnoldi factorization.   |
        !        %---------------------------------------%
        !
        getv0    = .true.
        update   = .false.
        ushift   = .false.
        cnorm    = .false.
        !
        if (info .ne. 0) then
            !
            !           %--------------------------------------------%
            !           | User provides the initial residual vector. |
            !           %--------------------------------------------%
            !
            initv = .true.
            info  = 0
        else
            initv = .false.
        end if
    end if
    !
    !     %---------------------------------------------%
    !     | Get a possibly random starting vector and   |
    !     | force it into the range of the operator OP. |
    !     %---------------------------------------------%
    !
10  continue
    !
    if (getv0) then
        call zgetv0(ido, bmat, 1, initv, n, 1, v, ldv, resid, rnorm, ipntr, workd, info)
        !
        if (ido .ne. 99) go to 9000
        !
        if (rnorm .eq. rzero) then
            !
            !           %-----------------------------------------%
            !           | The initial vector is zero. Error exit. |
            !           %-----------------------------------------%
            !
            info = -9
            go to 1100
        end if
        getv0 = .false.
        ido  = 0
    end if
    !
    !     %-----------------------------------%
    !     | Back from reverse communication : |
    !     | continue with update step         |
    !     %-----------------------------------%
    !
    if (update) go to 20
    !
    !     %-------------------------------------------%
    !     | Back from computing user specified shifts |
    !     %-------------------------------------------%
    !
    if (ushift) go to 50
    !
    !     %-------------------------------------%
    !     | Back from computing residual norm   |
    !     | at the end of the current iteration |
    !     %-------------------------------------%
    !
    if (cnorm)  go to 100
    !
    !     %----------------------------------------------------------%
    !     | Compute the first NEV steps of the Arnoldi factorization |
    !     %----------------------------------------------------------%
    !
    call znaitr(ido, bmat, n, 0, nev, mode, resid, rnorm, v, ldv, &
        &             h, ldh, ipntr, workd, info)
    !
    if (ido .ne. 99) go to 9000
    !
    if (info .gt. 0) then
        np   = info
        mxiter = iter
        info = -9999
        go to 1200
    end if
    !
    !     %--------------------------------------------------------------%
    !     |                                                              |
    !     |           M A I N  ARNOLDI  I T E R A T I O N  L O O P       |
    !     |           Each iteration implicitly restarts the Arnoldi     |
    !     |           factorization in place.                            |
    !     |                                                              |
    !     %--------------------------------------------------------------%
    !
1000 continue
    !
    iter = iter + 1
    !
    !if (msglvl .gt. 0) then
    !    call ivout(logfil, 1, iter, ndigit, '_naup2: **** Start of major iteration number ****')
    !end if
    !
    !        %-----------------------------------------------------------%
    !        | Compute NP additional steps of the Arnoldi factorization. |
    !        | Adjust NP since NEV might have been updated by last call  |
    !        | to the shift application routine znapps .                  |
    !        %-----------------------------------------------------------%
    !
    np  = kplusp - nev
    !
    !if (msglvl .gt. 1) then
    !    call ivout(logfil, 1, nev, ndigit, '_naup2: The length of the current Arnoldi factorization')
    !    call ivout(logfil, 1, np, ndigit, '_naup2: Extend the Arnoldi factorization by')
    !end if
    !
    !        %-----------------------------------------------------------%
    !        | Compute NP additional steps of the Arnoldi factorization. |
    !        %-----------------------------------------------------------%
    !
    ido = 0
20  continue
    update = .true.
    !
    call znaitr(ido, bmat, n, nev, np,    mode,  resid, rnorm, &
        &               v  , ldv , h, ldh, ipntr, workd, info)
    !
    if (ido .ne. 99) go to 9000
    !
    if (info .gt. 0) then
        np = info
        mxiter = iter
        info = -9999
        go to 1200
    end if
    update = .false.
    !
    !if (msglvl .gt. 1) then
    !    call dvout(logfil, 1, rnorm, ndigit, '_naup2: Corresponding B-norm of the residual')
    !end if
    !
    !        %--------------------------------------------------------%
    !        | Compute the eigenvalues and corresponding error bounds |
    !        | of the current upper Hessenberg matrix.                |
    !        %--------------------------------------------------------%
    !
    call zneigh(rnorm, kplusp, h, ldh, ritz, bounds, &
        &                q, ldq, workl, rwork,  ierr)
    !
    if (ierr .ne. 0) then
        info = -8
        go to 1200
    end if
    !
    !        %---------------------------------------------------%
    !        | Select the wanted Ritz values and their bounds    |
    !        | to be used in the convergence test.               |
    !        | The wanted part of the spectrum and corresponding |
    !        | error bounds are in the last NEV loc. of RITZ,    |
    !        | and BOUNDS respectively.                          |
    !        %---------------------------------------------------%
    !
    nev = nev0
    np = np0
    !
    !        %--------------------------------------------------%
    !        | Make a copy of Ritz values and the corresponding |
    !        | Ritz estimates obtained from zneigh .             |
    !        %--------------------------------------------------%
    !
    call zcopy(kplusp,ritz,1,workl(kplusp**2+1),1)
    call zcopy(kplusp,bounds,1,workl(kplusp**2+kplusp+1),1)
    !
    !        %---------------------------------------------------%
    !        | Select the wanted Ritz values and their bounds    |
    !        | to be used in the convergence test.               |
    !        | The wanted part of the spectrum and corresponding |
    !        | bounds are in the last NEV loc. of RITZ           |
    !        | BOUNDS respectively.                              |
    !        %---------------------------------------------------%
    !
    call zngets(ishift, which, nev, np, ritz, bounds)
    !
    !        %------------------------------------------------------------%
    !        | Convergence test: currently we use the following criteria. |
    !        | The relative accuracy of a Ritz value is considered        |
    !        | acceptable if:                                             |
    !        |                                                            |
    !        | error_bounds(i) .le. tol*max(eps23, magnitude_of_ritz(i)). |
    !        |                                                            |
    !        %------------------------------------------------------------%
    !
    nconv  = 0
    !
    do 25 i = 1, nev
        rtemp = max( eps23, dlapy2 ( dble (ritz(np+i)),  &
            &                                  dimag (ritz(np+i)) ) )
        if ( dlapy2 (dble (bounds(np+i)),dimag (bounds(np+i))) &
            &                 .le. tol*rtemp ) then
        nconv = nconv + 1
        end if
25  continue
    !
    if (msglvl .gt. 2) then
        kp(1) = nev
        kp(2) = np
        kp(3) = nconv
        call ivout(logfil, 3, kp, ndigit, '_naup2: NEV, NP, NCONV are')
        call zvout(logfil, kplusp, ritz, ndigit,'_naup2: The eigenvalues of H')
        call zvout(logfil, kplusp, bounds, ndigit, '_naup2: Ritz estimates of the current NCV Ritz values')
    end if
    !
    !        %---------------------------------------------------------%
    !        | Count the number of unwanted Ritz values that have zero |
    !        | Ritz estimates. If any Ritz estimates are equal to zero |
    !        | then a leading block of H of order equal to at least    |
    !        | the number of Ritz values with zero Ritz estimates has  |
    !        | split off. None of these Ritz values may be removed by  |
    !        | shifting. Decrease NP the number of shifts to apply. If |
    !        | no shifts may be applied, then prepare to exit          |
    !        %---------------------------------------------------------%
    !
    nptemp = np
    do 30 j=1, nptemp
        if (bounds(j) .eq. zero) then
            np = np - 1
            nev = nev + 1
        end if
30  continue
    !
    if ( (nconv .ge. nev0) .or. (iter .gt. mxiter) .or. (np .eq. 0) ) then
        !
        if (msglvl .gt. 4) then
            call zvout(logfil, kplusp, workl(kplusp**2+1), ndigit, '_naup2: Eigenvalues computed by _neigh:')
            call zvout(logfil, kplusp, workl(kplusp**2+kplusp+1),ndigit,'_naup2: Ritz estimates computed by _neigh:')
        end if
        !
        !           %------------------------------------------------%
        !           | Prepare to exit. Put the converged Ritz values |
        !           | and corresponding bounds in RITZ(1:NCONV) and  |
        !           | BOUNDS(1:NCONV) respectively. Then sort. Be    |
        !           | careful when NCONV > NP                        |
        !           %------------------------------------------------%
        !
        !           %------------------------------------------%
        !           |  Use h( 3,1 ) as storage to communicate  |
        !           |  rnorm to zneupd  if needed               |
        !           %------------------------------------------%

        h(3,1) = dcmplx (rnorm,rzero)
        !
        !           %----------------------------------------------%
        !           | Sort Ritz values so that converged Ritz      |
        !           | values appear within the first NEV locations |
        !           | of ritz and bounds, and the most desired one |
        !           | appears at the front.                        |
        !           %----------------------------------------------%
        !
        if (which .eq. 'LM') wprime = 'SM'
        if (which .eq. 'SM') wprime = 'LM'
        if (which .eq. 'LR') wprime = 'SR'
        if (which .eq. 'SR') wprime = 'LR'
        if (which .eq. 'LI') wprime = 'SI'
        if (which .eq. 'SI') wprime = 'LI'
        !
        call zsortc (wprime, .true., kplusp, ritz, bounds)
        !
        !           %--------------------------------------------------%
        !           | Scale the Ritz estimate of each Ritz value       |
        !           | by 1 / max(eps23, magnitude of the Ritz value).  |
        !           %--------------------------------------------------%
        !
        do 35 j = 1, nev0
            rtemp = max( eps23, dlapy2 ( dble (ritz(j)),     &
                &                                       dimag (ritz(j)) ) )
            bounds(j) = bounds(j)/rtemp
35      continue
        !
        !           %---------------------------------------------------%
        !           | Sort the Ritz values according to the scaled Ritz |
        !           | estimates.  This will push all the converged ones |
        !           | towards the front of ritz, bounds (in the case    |
        !           | when NCONV < NEV.)                                |
        !           %---------------------------------------------------%
        !
        wprime = 'LM'
        call zsortc (wprime, .true., nev0, bounds, ritz)
        !
        !           %----------------------------------------------%
        !           | Scale the Ritz estimate back to its original |
        !           | value.                                       |
        !           %----------------------------------------------%
        !
        do 40 j = 1, nev0
            rtemp = max( eps23, dlapy2 ( dble (ritz(j)),    &
                &                                       dimag (ritz(j)) ) )
            bounds(j) = bounds(j)*rtemp
40      continue
        !
        !           %-----------------------------------------------%
        !           | Sort the converged Ritz values again so that  |
        !           | the "threshold" value appears at the front of |
        !           | ritz and bound.                               |
        !           %-----------------------------------------------%
        !
        call zsortc (which, .true., nconv, ritz, bounds)
        !
        if (msglvl .gt. 1) then
            call zvout  (logfil, kplusp, ritz, ndigit,'_naup2: Sorted eigenvalues')
            call zvout  (logfil, kplusp, bounds, ndigit,'_naup2: Sorted ritz estimates.')
        end if
        !
        !           %------------------------------------%
        !           | Max iterations have been exceeded. |
        !           %------------------------------------%
        !
        if (iter .gt. mxiter .and. nconv .lt. nev0) info = 1
        !
        !           %---------------------%
        !           | No shifts to apply. |
        !           %---------------------%
        !
        if (np .eq. 0 .and. nconv .lt. nev0)  info = 2
        !
        np = nconv
        go to 1100
        !
    else if ( (nconv .lt. nev0) .and. (ishift .eq. 1) ) then
        !
        !           %-------------------------------------------------%
        !           | Do not have all the requested eigenvalues yet.  |
        !           | To prevent possible stagnation, adjust the size |
        !           | of NEV.                                         |
        !           %-------------------------------------------------%
        !
        nevbef = nev
        nev = nev + min(nconv, np/2)
        if (nev .eq. 1 .and. kplusp .ge. 6) then
            nev = kplusp / 2
        else if (nev .eq. 1 .and. kplusp .gt. 3) then
            nev = 2
        end if
        np = kplusp - nev
        !
        !           %---------------------------------------%
        !           | If the size of NEV was just increased |
        !           | resort the eigenvalues.               |
        !           %---------------------------------------%
        !
        if (nevbef .lt. nev) then
            call zngets  (ishift, which, nev, np, ritz, bounds)
        endif
    endif
    !
    if (msglvl .gt. 0) then
    !    call ivout (logfil, 1, nconv, ndigit,'_naup2: no. of "converged" Ritz values at this iter.')
        if (msglvl .gt. 1) then
            kp(1) = nev
            kp(2) = np
            call ivout (logfil, 2, kp, ndigit,'_naup2: NEV and NP are')
            call zvout  (logfil, nev, ritz(np+1), ndigit,'_naup2: "wanted" Ritz values ')
            call zvout  (logfil, nev, bounds(np+1), ndigit,'_naup2: Ritz estimates of the "wanted" values ')
        end if
    end if
    !
    if (ishift .eq. 0) then
        !
        !           %-------------------------------------------------------%
        !           | User specified shifts: pop back out to get the shifts |
        !           | and return them in the first 2*NP locations of WORKL. |
        !           %-------------------------------------------------------%
        !
        ushift = .true.
        ido = 3
        go to 9000
    end if
50  continue
    ushift = .false.
    !
    if ( ishift .ne. 1 ) then
        !
        !            %----------------------------------%
        !            | Move the NP shifts from WORKL to |
        !            | RITZ, to free up WORKL           |
        !            | for non-exact shift case.        |
        !            %----------------------------------%
        !
        call zcopy  (np, workl, 1, ritz, 1)
    end if
    !
    if (msglvl .gt. 2) then
    !    call ivout(logfil, 1, np, ndigit, '_naup2: The number of shifts to apply ')
        call zvout(logfil, np, ritz, ndigit, '_naup2: values of the shifts')
        if ( ishift .eq. 1 ) then
            call zvout(logfil, np, bounds, ndigit, '_naup2: Ritz estimates of the shifts')
        endif
    end if
    !
    !        %---------------------------------------------------------%
    !        | Apply the NP implicit shifts by QR bulge chasing.       |
    !        | Each shift is applied to the whole upper Hessenberg     |
    !        | matrix H.                                               |
    !        | The first 2*N locations of WORKD are used as workspace. |
    !        %---------------------------------------------------------%
    !
    call znapps  (n, nev, np, ritz, v, ldv, h, ldh, resid, q, ldq, workl, workd)
    !
    !        %---------------------------------------------%
    !        | Compute the B-norm of the updated residual. |
    !        | Keep B*RESID in WORKD(1:N) to be used in    |
    !        | the first step of the next call to znaitr .  |
    !        %---------------------------------------------%
    !
    cnorm = .true.
    call arscnd (t2)
    if (bmat .eq. 'G') then
        nbx = nbx + 1
        call zcopy  (n, resid, 1, workd(n+1), 1)
        ipntr(1) = n + 1
        ipntr(2) = 1
        ido = 2
        !
        !           %----------------------------------%
        !           | Exit in order to compute B*RESID |
        !           %----------------------------------%
        !
        go to 9000
    else if (bmat .eq. 'I') then
        call zcopy  (n, resid, 1, workd, 1)
    end if
    !
100 continue
    !
    !        %----------------------------------%
    !        | Back from reverse communication; |
    !        | WORKD(1:N) := B*RESID            |
    !        %----------------------------------%
    !
    if (bmat .eq. 'G') then
        call arscnd (t3)
        tmvbx = tmvbx + (t3 - t2)
    end if
    !
    if (bmat .eq. 'G') then
        cmpnorm = zdotc  (n, resid, 1, workd, 1)
        rnorm = sqrt(dlapy2 (dble (cmpnorm),dimag (cmpnorm)))
    else if (bmat .eq. 'I') then
        rnorm = dznrm2 (n, resid, 1)
    end if
    cnorm = .false.
    !
    if (msglvl .gt. 2) then
        !call dvout  (logfil, 1, rnorm, ndigit, '_naup2: B-norm of residual for compressed factorization')
        call zmout  (logfil, nev, nev, h, ldh, ndigit, '_naup2: Compressed upper Hessenberg matrix H')
    end if
    !
    go to 1000
    !
    !     %---------------------------------------------------------------%
    !     |                                                               |
    !     |  E N D     O F     M A I N     I T E R A T I O N     L O O P  |
    !     |                                                               |
    !     %---------------------------------------------------------------%
    !
1100 continue
    !
    mxiter = iter
    nev = nconv
    !
1200 continue
    ido = 99
    !
    !     %------------%
    !     | Error Exit |
    !     %------------%
    !
    call arscnd (t1)
    tcaup2 = t1 - t0
    !
9000 continue
    !
    !     %---------------%
    !     | End of znaup2  |
    !     %---------------%
    !
    return
    END SUBROUTINE znaup2