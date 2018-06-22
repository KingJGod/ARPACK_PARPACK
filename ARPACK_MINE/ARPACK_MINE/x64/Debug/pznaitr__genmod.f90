        !COMPILER-GENERATED INTERFACE MODULE: Thu Jun 21 16:38:14 2018
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE PZNAITR__genmod
          INTERFACE 
            SUBROUTINE PZNAITR(COMM,IDO,BMAT,N,K,NP,NB,RESID,RNORM,V,LDV&
     &,H,LDH,IPNTR,WORKD,WORKL,INFO)
              INTEGER(KIND=4), INTENT(IN) :: LDH
              INTEGER(KIND=4), INTENT(IN) :: LDV
              INTEGER(KIND=4), INTENT(IN) :: NP
              INTEGER(KIND=4), INTENT(IN) :: K
              INTEGER(KIND=4), INTENT(IN) :: N
              INTEGER(KIND=4), INTENT(INOUT) :: COMM
              INTEGER(KIND=4), INTENT(INOUT) :: IDO
              CHARACTER(LEN=1), INTENT(IN) :: BMAT
              INTEGER(KIND=4), INTENT(IN) :: NB
              COMPLEX(KIND=8), INTENT(INOUT) :: RESID(N)
              REAL(KIND=8), INTENT(INOUT) :: RNORM
              COMPLEX(KIND=8), INTENT(INOUT) :: V(LDV,K+NP)
              COMPLEX(KIND=8), INTENT(INOUT) :: H(LDH,K+NP)
              INTEGER(KIND=4), INTENT(OUT) :: IPNTR(3)
              COMPLEX(KIND=8), INTENT(INOUT) :: WORKD(3*N)
              COMPLEX(KIND=8), INTENT(INOUT) :: WORKL(2*LDH)
              INTEGER(KIND=4), INTENT(OUT) :: INFO
            END SUBROUTINE PZNAITR
          END INTERFACE 
        END MODULE PZNAITR__genmod
