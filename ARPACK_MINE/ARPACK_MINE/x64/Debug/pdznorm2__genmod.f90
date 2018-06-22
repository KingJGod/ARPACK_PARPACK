        !COMPILER-GENERATED INTERFACE MODULE: Thu Jun 21 16:38:13 2018
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE PDZNORM2__genmod
          INTERFACE 
            FUNCTION PDZNORM2(COMM,N,X,INC)
              INTEGER(KIND=4) :: N
              INTEGER(KIND=4) :: COMM
              COMPLEX(KIND=8) :: X(N)
              INTEGER(KIND=4) :: INC
              REAL(KIND=8) :: PDZNORM2
            END FUNCTION PDZNORM2
          END INTERFACE 
        END MODULE PDZNORM2__genmod
