        !COMPILER-GENERATED INTERFACE MODULE: Thu Jun 21 16:36:41 2018
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE ZNAUPD__genmod
          INTERFACE 
            SUBROUTINE ZNAUPD(IDO,BMAT,N,WHICH,NEV,TOL,RESID,NCV,V,LDV, &
     &IPARAM,IPNTR,WORKD,WORKL,LWORKL,RWORK,INFO)
              INTEGER(KIND=4), INTENT(IN) :: LWORKL
              INTEGER(KIND=4), INTENT(IN) :: LDV
              INTEGER(KIND=4), INTENT(IN) :: NCV
              INTEGER(KIND=4), INTENT(IN) :: N
              INTEGER(KIND=4), INTENT(INOUT) :: IDO
              CHARACTER(LEN=1), INTENT(IN) :: BMAT
              CHARACTER(LEN=2), INTENT(IN) :: WHICH
              INTEGER(KIND=4), INTENT(IN) :: NEV
              REAL(KIND=8), INTENT(INOUT) :: TOL
              COMPLEX(KIND=8), INTENT(INOUT) :: RESID(N)
              COMPLEX(KIND=8), INTENT(OUT) :: V(LDV,NCV)
              INTEGER(KIND=4), INTENT(INOUT) :: IPARAM(11)
              INTEGER(KIND=4), INTENT(OUT) :: IPNTR(14)
              COMPLEX(KIND=8), INTENT(INOUT) :: WORKD(3*N)
              COMPLEX(KIND=8), INTENT(INOUT) :: WORKL(LWORKL)
              REAL(KIND=8), INTENT(INOUT) :: RWORK(NCV)
              INTEGER(KIND=4), INTENT(INOUT) :: INFO
            END SUBROUTINE ZNAUPD
          END INTERFACE 
        END MODULE ZNAUPD__genmod
