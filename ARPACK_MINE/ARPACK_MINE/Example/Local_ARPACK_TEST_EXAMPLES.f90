SUBROUTINE LOCAL_ARPACK_TEST_EXAMPLES
    IMPLICIT NONE
    INTEGER( kind = 4 )::Index_Test_Case, IF_CONTINUE
    WRITE(*,*)"This is the test for F90 version COMPLEX ARPACK Modified By XYC"
    
1010 CONTINUE
    WRITE(*,*)"Input the test case number:"
    WRITE(*,*)"1. For Simplest Case of Using ARPACK"
    WRITE(*,*)"2. For standard eigenvalue case of using ARPACK"
    WRITE(*,*)"3. For standard eigenvalue case of Shift-invert mode Using ARPACK"
    WRITE(*,*)"4. For generalized eigenvalue case of using ARPACK"
    WRITE(*,*)"5. For generalized eigenvalue case of Shift-invert mode Using ARPACK"
    WRITE(*,*)"6. For Parallel test"    
    READ(*,*)Index_Test_Case
    IF((Index_Test_Case .NE. 1).AND.(Index_Test_Case .NE. 2).AND. & 
       & (Index_Test_Case .NE. 3).AND.(Index_Test_Case .NE. 4).AND. &  
       &   (Index_Test_Case .NE. 5).AND.(Index_Test_Case .NE. 6)) THEN
    WRITE(*,*)
    WRITE(*,*)"YOU HAVE INOUT WRONG EXAMPLE CASE NUMBER !!!!!!!!!!!!!!!!!!!"
    GOTO 1010
    END IF
    
    SELECT CASE(Index_Test_Case)

    CASE(1)
        !============================EXAMPLE_CASE1==================================
        !
        !     This example program is intended to illustrate the
        !     simplest case of using ARPACK in considerable detail.
        !     This code may be used to understand basic usage of ARPACK
        !     and as a template for creating an interface to ARPACK.
        !
        !     This code shows how to use ARPACK to find a few eigenvalues
        !     (lambda) and corresponding eigenvectors (x) for the standard
        !     eigenvalue problem:
        !
        !                        A*x = lambda*x
        !
        !     where A is a general n by n complex matrix.
        !
        WRITE(*,*)
        WRITE(*,*)"1. For Simplest Case of Using ARPACK"
        CALL ZNSIMP_CASE1
        !
        !===========================================================================
        
    CASE(2)
        !============================EXAMPLE_CASE2==================================
        !     Example program to illustrate the idea of reverse communication
        !     for a standard complex nonsymmetric eigenvalue problem.
        WRITE(*,*)
        WRITE(*,*)"2. For standard eigenvalue case of using ARPACK"
        CALL ZNDRV1
        !
        !===========================================================================
        
    CASE(3)
        !============================EXAMPLE_CASE3==================================
        !     Simple program to illustrate the idea of reverse communication
        !     in shift-invert mode for a standard complex nonsymmetric eigenvalue
        !     problem.
        WRITE(*,*)
        WRITE(*,*)"3. For standard eigenvalue case of Shift-invert mode Using ARPACK"
        CALL ZNDRV2
        !
        !===========================================================================
    CASE(4)
        !============================EXAMPLE_CASE4==================================
        WRITE(*,*)
        WRITE(*,*)"4. For generalized eigenvalue case of using ARPACK"
        CALL ZNDRV3
        !
        !===========================================================================
    CASE(5)
        !============================EXAMPLE_CASE5==================================
        WRITE(*,*)
        WRITE(*,*)"5. For generalized eigenvalue case of Shift-invert mode Using ARPACK"
        CALL ZNDRV4
        !
        !===========================================================================
    CASE(6)
        !============================EXAMPLE_CASE6==================================
        WRITE(*,*)
        WRITE(*,*)"6. For Parallel test" 
        GOTO 1020
        !
        !===========================================================================
    END SELECT
    WRITE(*,*)
    WRITE(*,*)"CONTINUE TEST or NOT:1 for continue;2 not"
    READ(*,*)IF_CONTINUE
    IF(IF_CONTINUE.EQ.1) GOTO 1010
1020 CONTINUE
    END SUBROUTINE LOCAL_ARPACK_TEST_EXAMPLES