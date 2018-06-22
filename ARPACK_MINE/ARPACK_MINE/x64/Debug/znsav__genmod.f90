        !COMPILER-GENERATED INTERFACE MODULE: Fri Jun 22 12:23:16 2018
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE ZNSAV__genmod
          INTERFACE 
            SUBROUTINE ZNSAV(NX,V,W)
              INTEGER(KIND=4) :: NX
              COMPLEX(KIND=8) :: V(NX*NX)
              COMPLEX(KIND=8) :: W(NX*NX)
            END SUBROUTINE ZNSAV
          END INTERFACE 
        END MODULE ZNSAV__genmod
