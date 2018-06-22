        !COMPILER-GENERATED INTERFACE MODULE: Thu Jun 21 16:36:37 2018
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE PSVOUT__genmod
          INTERFACE 
            SUBROUTINE PSVOUT(COMM,LOUT,N,SX,IDIGIT,IFMT)
              INTEGER(KIND=4) :: COMM
              INTEGER(KIND=4) :: LOUT
              INTEGER(KIND=4) :: N
              REAL(KIND=4) :: SX(*)
              INTEGER(KIND=4) :: IDIGIT
              CHARACTER(*) :: IFMT
            END SUBROUTINE PSVOUT
          END INTERFACE 
        END MODULE PSVOUT__genmod
