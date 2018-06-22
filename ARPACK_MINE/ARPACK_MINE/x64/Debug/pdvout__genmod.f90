        !COMPILER-GENERATED INTERFACE MODULE: Thu Jun 21 16:36:39 2018
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE PDVOUT__genmod
          INTERFACE 
            SUBROUTINE PDVOUT(COMM,LOUT,N,SX,IDIGIT,IFMT)
              INTEGER(KIND=4) :: COMM
              INTEGER(KIND=4) :: LOUT
              INTEGER(KIND=4) :: N
              REAL(KIND=8) :: SX(*)
              INTEGER(KIND=4) :: IDIGIT
              CHARACTER(*) :: IFMT
            END SUBROUTINE PDVOUT
          END INTERFACE 
        END MODULE PDVOUT__genmod
