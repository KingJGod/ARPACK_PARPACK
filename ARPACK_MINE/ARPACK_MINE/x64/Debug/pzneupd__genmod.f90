        !COMPILER-GENERATED INTERFACE MODULE: Thu Jun 21 16:36:40 2018
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE PZNEUPD__genmod
          INTERFACE 
            SUBROUTINE PZNEUPD(COMM,RVEC,HOWMNY,SELECT,D,Z,LDZ,SIGMA,   &
     &WORKEV,BMAT,N,WHICH,NEV,TOL,RESID,NCV,V,LDV,IPARAM,IPNTR,WORKD,   &
     &WORKL,LWORKL,RWORK,INFO)
              INTEGER(KIND=4), INTENT(IN) :: LWORKL
              INTEGER(KIND=4), INTENT(IN) :: LDV
              INTEGER(KIND=4), INTENT(IN) :: NCV
              INTEGER(KIND=4), INTENT(IN) :: NEV
              INTEGER(KIND=4), INTENT(IN) :: N
              INTEGER(KIND=4), INTENT(IN) :: LDZ
              INTEGER(KIND=4), INTENT(INOUT) :: COMM
              LOGICAL(KIND=4), INTENT(IN) :: RVEC
              CHARACTER(LEN=1), INTENT(IN) :: HOWMNY
              LOGICAL(KIND=4), INTENT(INOUT) :: SELECT(NCV)
              COMPLEX(KIND=8), INTENT(OUT) :: D(NEV)
              COMPLEX(KIND=8), INTENT(OUT) :: Z(LDZ,NEV)
              COMPLEX(KIND=8), INTENT(IN) :: SIGMA
              COMPLEX(KIND=8), INTENT(INOUT) :: WORKEV(2*NCV)
              CHARACTER(LEN=1), INTENT(IN) :: BMAT
              CHARACTER(LEN=2), INTENT(IN) :: WHICH
              REAL(KIND=8), INTENT(INOUT) :: TOL
              COMPLEX(KIND=8), INTENT(INOUT) :: RESID(N)
              COMPLEX(KIND=8), INTENT(INOUT) :: V(LDV,NCV)
              INTEGER(KIND=4), INTENT(INOUT) :: IPARAM(11)
              INTEGER(KIND=4), INTENT(OUT) :: IPNTR(14)
              COMPLEX(KIND=8), INTENT(INOUT) :: WORKD(3*N)
              COMPLEX(KIND=8), INTENT(INOUT) :: WORKL(LWORKL)
              REAL(KIND=8), INTENT(INOUT) :: RWORK(NCV)
              INTEGER(KIND=4), INTENT(INOUT) :: INFO
            END SUBROUTINE PZNEUPD
          END INTERFACE 
        END MODULE PZNEUPD__genmod
