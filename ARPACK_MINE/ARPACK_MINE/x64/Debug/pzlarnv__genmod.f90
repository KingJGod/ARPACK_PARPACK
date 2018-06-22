        !COMPILER-GENERATED INTERFACE MODULE: Thu Jun 21 16:43:00 2018
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE PZLARNV__genmod
          INTERFACE 
            SUBROUTINE PZLARNV(COMM,IDIST,ISEED,N,X)
              INTEGER(KIND=4) :: COMM
              INTEGER(KIND=4) :: IDIST
              INTEGER(KIND=4) :: ISEED(4)
              INTEGER(KIND=4) :: N
              COMPLEX(KIND=8) :: X(*)
            END SUBROUTINE PZLARNV
          END INTERFACE 
        END MODULE PZLARNV__genmod
