        !COMPILER-GENERATED INTERFACE MODULE: Thu Jun 21 16:36:41 2018
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE PZNGETS__genmod
          INTERFACE 
            SUBROUTINE PZNGETS(COMM,ISHIFT,WHICH,KEV,NP,RITZ,BOUNDS)
              INTEGER(KIND=4), INTENT(IN) :: NP
              INTEGER(KIND=4), INTENT(IN) :: KEV
              INTEGER(KIND=4), INTENT(INOUT) :: COMM
              INTEGER(KIND=4), INTENT(IN) :: ISHIFT
              CHARACTER(LEN=2), INTENT(IN) :: WHICH
              COMPLEX(KIND=8), INTENT(INOUT) :: RITZ(KEV+NP)
              COMPLEX(KIND=8), INTENT(INOUT) :: BOUNDS(KEV+NP)
            END SUBROUTINE PZNGETS
          END INTERFACE 
        END MODULE PZNGETS__genmod
