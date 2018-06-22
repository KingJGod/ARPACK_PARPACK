        !COMPILER-GENERATED INTERFACE MODULE: Thu Jun 21 16:36:36 2018
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE ZSORTC__genmod
          INTERFACE 
            SUBROUTINE ZSORTC(WHICH,APPLY,N,X,Y)
              INTEGER(KIND=4), INTENT(IN) :: N
              CHARACTER(LEN=2), INTENT(IN) :: WHICH
              LOGICAL(KIND=4), INTENT(IN) :: APPLY
              COMPLEX(KIND=8), INTENT(INOUT) :: X(0:N-1)
              COMPLEX(KIND=8), INTENT(INOUT) :: Y(0:N-1)
            END SUBROUTINE ZSORTC
          END INTERFACE 
        END MODULE ZSORTC__genmod
