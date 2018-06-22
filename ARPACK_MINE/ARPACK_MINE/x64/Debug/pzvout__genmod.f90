        !COMPILER-GENERATED INTERFACE MODULE: Thu Jun 21 16:36:36 2018
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE PZVOUT__genmod
          INTERFACE 
            SUBROUTINE PZVOUT(COMM,LOUT,N,CX,IDIGIT,IFMT)
              INTEGER(KIND=4) :: COMM
              INTEGER(KIND=4) :: LOUT
              INTEGER(KIND=4) :: N
              COMPLEX(KIND=8) :: CX(*)
              INTEGER(KIND=4) :: IDIGIT
              CHARACTER(*) :: IFMT
            END SUBROUTINE PZVOUT
          END INTERFACE 
        END MODULE PZVOUT__genmod
