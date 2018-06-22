        !COMPILER-GENERATED INTERFACE MODULE: Thu Jun 21 16:36:40 2018
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE PZNEIGH__genmod
          INTERFACE 
            SUBROUTINE PZNEIGH(COMM,RNORM,N,H,LDH,RITZ,BOUNDS,Q,LDQ,    &
     &WORKL,RWORK,IERR)
              INTEGER(KIND=4), INTENT(IN) :: LDQ
              INTEGER(KIND=4), INTENT(IN) :: LDH
              INTEGER(KIND=4), INTENT(IN) :: N
              INTEGER(KIND=4), INTENT(INOUT) :: COMM
              REAL(KIND=8), INTENT(IN) :: RNORM
              COMPLEX(KIND=8), INTENT(IN) :: H(LDH,N)
              COMPLEX(KIND=8), INTENT(OUT) :: RITZ(N)
              COMPLEX(KIND=8), INTENT(OUT) :: BOUNDS(N)
              COMPLEX(KIND=8), INTENT(INOUT) :: Q(LDQ,N)
              COMPLEX(KIND=8), INTENT(INOUT) :: WORKL(N*(N+3))
              REAL(KIND=8), INTENT(INOUT) :: RWORK(N)
              INTEGER(KIND=4), INTENT(OUT) :: IERR
            END SUBROUTINE PZNEIGH
          END INTERFACE 
        END MODULE PZNEIGH__genmod
