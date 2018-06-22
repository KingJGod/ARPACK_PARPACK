        !COMPILER-GENERATED INTERFACE MODULE: Thu Jun 21 16:36:36 2018
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE PZMOUT__genmod
          INTERFACE 
            SUBROUTINE PZMOUT(COMM,LOUT,M,N,A,LDA,IDIGIT,IFMT)
              INTEGER(KIND=4) :: LDA
              INTEGER(KIND=4) :: COMM
              INTEGER(KIND=4) :: LOUT
              INTEGER(KIND=4) :: M
              INTEGER(KIND=4) :: N
              COMPLEX(KIND=8) :: A(LDA,*)
              INTEGER(KIND=4) :: IDIGIT
              CHARACTER(*) :: IFMT
            END SUBROUTINE PZMOUT
          END INTERFACE 
        END MODULE PZMOUT__genmod
