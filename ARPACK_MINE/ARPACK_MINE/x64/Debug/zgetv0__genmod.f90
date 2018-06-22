        !COMPILER-GENERATED INTERFACE MODULE: Thu Jun 21 16:36:41 2018
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE ZGETV0__genmod
          INTERFACE 
            SUBROUTINE ZGETV0(IDO,BMAT,ITRY,INITV,N,J,V,LDV,RESID,RNORM,&
     &IPNTR,WORKD,IERR)
              INTEGER(KIND=4), INTENT(IN) :: LDV
              INTEGER(KIND=4), INTENT(IN) :: J
              INTEGER(KIND=4), INTENT(IN) :: N
              INTEGER(KIND=4), INTENT(INOUT) :: IDO
              CHARACTER(LEN=1), INTENT(IN) :: BMAT
              INTEGER(KIND=4), INTENT(IN) :: ITRY
              LOGICAL(KIND=4), INTENT(IN) :: INITV
              COMPLEX(KIND=8), INTENT(IN) :: V(LDV,J)
              COMPLEX(KIND=8), INTENT(INOUT) :: RESID(N)
              REAL(KIND=8), INTENT(OUT) :: RNORM
              INTEGER(KIND=4), INTENT(OUT) :: IPNTR(3)
              COMPLEX(KIND=8), INTENT(INOUT) :: WORKD(2*N)
              INTEGER(KIND=4), INTENT(OUT) :: IERR
            END SUBROUTINE ZGETV0
          END INTERFACE 
        END MODULE ZGETV0__genmod
