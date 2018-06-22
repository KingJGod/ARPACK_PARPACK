        !COMPILER-GENERATED INTERFACE MODULE: Thu Jun 21 16:45:06 2018
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE PARALLEL_AV__genmod
          INTERFACE 
            SUBROUTINE PARALLEL_AV(COMM,NLOC,NX,MV_BUF,V,W)
              INTEGER(KIND=4) :: NX
              INTEGER(KIND=4) :: NLOC
              INTEGER(KIND=4) :: COMM
              REAL(KIND=8) :: MV_BUF(NX)
              REAL(KIND=8) :: V(NLOC)
              REAL(KIND=8) :: W(NLOC)
            END SUBROUTINE PARALLEL_AV
          END INTERFACE 
        END MODULE PARALLEL_AV__genmod
