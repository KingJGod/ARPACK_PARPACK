        !COMPILER-GENERATED INTERFACE MODULE: Fri Jun 22 12:16:29 2018
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE PZNDRV1AV__genmod
          INTERFACE 
            SUBROUTINE PZNDRV1AV(COMM,NLOC,NX,MV_BUF,V,W)
              INTEGER(KIND=4) :: NX
              INTEGER(KIND=4) :: NLOC
              INTEGER(KIND=4) :: COMM
              COMPLEX(KIND=8) :: MV_BUF(NX)
              COMPLEX(KIND=8) :: V(NLOC)
              COMPLEX(KIND=8) :: W(NLOC)
            END SUBROUTINE PZNDRV1AV
          END INTERFACE 
        END MODULE PZNDRV1AV__genmod
