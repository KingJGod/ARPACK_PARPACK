        !COMPILER-GENERATED INTERFACE MODULE: Thu Jun 21 16:36:39 2018
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE SMOUT__genmod
          INTERFACE 
            SUBROUTINE SMOUT(LOUT,M,N,A,LDA,IDIGIT,IFMT)
              INTEGER(KIND=4) :: LDA
              INTEGER(KIND=4) :: LOUT
              INTEGER(KIND=4) :: M
              INTEGER(KIND=4) :: N
              REAL(KIND=4) :: A(LDA,*)
              INTEGER(KIND=4) :: IDIGIT
              CHARACTER(*) :: IFMT
            END SUBROUTINE SMOUT
          END INTERFACE 
        END MODULE SMOUT__genmod
