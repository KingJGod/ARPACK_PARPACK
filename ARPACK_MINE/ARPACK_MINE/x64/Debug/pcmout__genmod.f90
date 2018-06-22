        !COMPILER-GENERATED INTERFACE MODULE: Thu Jun 21 16:36:38 2018
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE PCMOUT__genmod
          INTERFACE 
            SUBROUTINE PCMOUT(COMM,LOUT,M,N,A,LDA,IDIGIT,IFMT)
              INTEGER(KIND=4) :: LDA
              INTEGER(KIND=4) :: COMM
              INTEGER(KIND=4) :: LOUT
              INTEGER(KIND=4) :: M
              INTEGER(KIND=4) :: N
              COMPLEX(KIND=4) :: A(LDA,*)
              INTEGER(KIND=4) :: IDIGIT
              CHARACTER(*) :: IFMT
            END SUBROUTINE PCMOUT
          END INTERFACE 
        END MODULE PCMOUT__genmod
