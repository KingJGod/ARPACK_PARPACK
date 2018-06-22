        !COMPILER-GENERATED INTERFACE MODULE: Fri Jun 22 12:16:28 2018
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE AV__genmod
          INTERFACE 
            SUBROUTINE AV(NX,V,W)
              INTEGER(KIND=4) :: NX
              COMPLEX(KIND=8) :: V(NX*NX)
              COMPLEX(KIND=8) :: W(NX*NX)
            END SUBROUTINE AV
          END INTERFACE 
        END MODULE AV__genmod
