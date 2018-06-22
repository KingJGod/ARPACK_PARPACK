        !COMPILER-GENERATED INTERFACE MODULE: Thu Jun 21 16:36:40 2018
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE ZNAPPS__genmod
          INTERFACE 
            SUBROUTINE ZNAPPS(N,KEV,NP,SHIFT,V,LDV,H,LDH,RESID,Q,LDQ,   &
     &WORKL,WORKD)
              INTEGER(KIND=4), INTENT(IN) :: LDQ
              INTEGER(KIND=4), INTENT(IN) :: LDH
              INTEGER(KIND=4), INTENT(IN) :: LDV
              INTEGER(KIND=4), INTENT(IN) :: NP
              INTEGER(KIND=4), INTENT(INOUT) :: KEV
              INTEGER(KIND=4), INTENT(IN) :: N
              COMPLEX(KIND=8), INTENT(IN) :: SHIFT(NP)
              COMPLEX(KIND=8), INTENT(INOUT) :: V(LDV,KEV+NP)
              COMPLEX(KIND=8), INTENT(INOUT) :: H(LDH,KEV+NP)
              COMPLEX(KIND=8), INTENT(INOUT) :: RESID(N)
              COMPLEX(KIND=8), INTENT(INOUT) :: Q(LDQ,KEV+NP)
              COMPLEX(KIND=8), INTENT(INOUT) :: WORKL(KEV+NP)
              COMPLEX(KIND=8), INTENT(INOUT) :: WORKD(2*N)
            END SUBROUTINE ZNAPPS
          END INTERFACE 
        END MODULE ZNAPPS__genmod
