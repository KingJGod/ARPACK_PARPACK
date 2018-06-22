        !COMPILER-GENERATED INTERFACE MODULE: Thu Jun 21 16:36:42 2018
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE PZNAUP2__genmod
          INTERFACE 
            SUBROUTINE PZNAUP2(COMM,IDO,BMAT,N,WHICH,NEV,NP,TOL,RESID,  &
     &MODE,IUPD,ISHIFT,MXITER,V,LDV,H,LDH,RITZ,BOUNDS,Q,LDQ,WORKL,IPNTR,&
     &WORKD,RWORK,INFO)
              INTEGER(KIND=4), INTENT(IN) :: LDQ
              INTEGER(KIND=4), INTENT(IN) :: LDH
              INTEGER(KIND=4), INTENT(IN) :: LDV
              INTEGER(KIND=4), INTENT(INOUT) :: NP
              INTEGER(KIND=4), INTENT(INOUT) :: NEV
              INTEGER(KIND=4), INTENT(IN) :: N
              INTEGER(KIND=4), INTENT(INOUT) :: COMM
              INTEGER(KIND=4), INTENT(INOUT) :: IDO
              CHARACTER(LEN=1), INTENT(IN) :: BMAT
              CHARACTER(LEN=2), INTENT(IN) :: WHICH
              REAL(KIND=8), INTENT(INOUT) :: TOL
              COMPLEX(KIND=8), INTENT(INOUT) :: RESID(N)
              INTEGER(KIND=4), INTENT(INOUT) :: MODE
              INTEGER(KIND=4), INTENT(IN) :: IUPD
              INTEGER(KIND=4), INTENT(INOUT) :: ISHIFT
              INTEGER(KIND=4), INTENT(INOUT) :: MXITER
              COMPLEX(KIND=8), INTENT(INOUT) :: V(LDV,NEV+NP)
              COMPLEX(KIND=8), INTENT(OUT) :: H(LDH,NEV+NP)
              COMPLEX(KIND=8), INTENT(OUT) :: RITZ(NEV+NP)
              COMPLEX(KIND=8), INTENT(OUT) :: BOUNDS(NEV+NP)
              COMPLEX(KIND=8), INTENT(INOUT) :: Q(LDQ,NEV+NP)
              COMPLEX(KIND=8), INTENT(INOUT) :: WORKL((NEV+NP)*(NEV+NP+3&
     &))
              INTEGER(KIND=4), INTENT(OUT) :: IPNTR(13)
              COMPLEX(KIND=8), INTENT(INOUT) :: WORKD(3*N)
              REAL(KIND=8), INTENT(INOUT) :: RWORK(NEV+NP)
              INTEGER(KIND=4), INTENT(INOUT) :: INFO
            END SUBROUTINE PZNAUP2
          END INTERFACE 
        END MODULE PZNAUP2__genmod
