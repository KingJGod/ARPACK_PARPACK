        !COMPILER-GENERATED INTERFACE MODULE: Thu Jun 21 16:36:36 2018
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE DMOUT__genmod
          INTERFACE 
            SUBROUTINE DMOUT(LOUT,M,N,A,LDA,IDIGIT,IFMT)
              INTEGER(KIND=4) :: LDA
              INTEGER(KIND=4) :: LOUT
              INTEGER(KIND=4) :: M
              INTEGER(KIND=4) :: N
              REAL(KIND=8) :: A(LDA,*)
              INTEGER(KIND=4) :: IDIGIT
              CHARACTER(*) :: IFMT
            END SUBROUTINE DMOUT
          END INTERFACE 
        END MODULE DMOUT__genmod
