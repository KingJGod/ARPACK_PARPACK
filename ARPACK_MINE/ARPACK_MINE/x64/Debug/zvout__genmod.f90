        !COMPILER-GENERATED INTERFACE MODULE: Thu Jun 21 16:36:37 2018
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE ZVOUT__genmod
          INTERFACE 
            SUBROUTINE ZVOUT(LOUT,N,CX,IDIGIT,IFMT)
              INTEGER(KIND=4) :: LOUT
              INTEGER(KIND=4) :: N
              COMPLEX(KIND=8) :: CX(*)
              INTEGER(KIND=4) :: IDIGIT
              CHARACTER(*) :: IFMT
            END SUBROUTINE ZVOUT
          END INTERFACE 
        END MODULE ZVOUT__genmod
