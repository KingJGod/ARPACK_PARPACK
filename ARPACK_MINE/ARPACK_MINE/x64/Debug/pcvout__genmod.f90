        !COMPILER-GENERATED INTERFACE MODULE: Thu Jun 21 16:36:37 2018
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE PCVOUT__genmod
          INTERFACE 
            SUBROUTINE PCVOUT(COMM,LOUT,N,CX,IDIGIT,IFMT)
              INTEGER(KIND=4) :: COMM
              INTEGER(KIND=4) :: LOUT
              INTEGER(KIND=4) :: N
              COMPLEX(KIND=4) :: CX(*)
              INTEGER(KIND=4) :: IDIGIT
              CHARACTER(*) :: IFMT
            END SUBROUTINE PCVOUT
          END INTERFACE 
        END MODULE PCVOUT__genmod
