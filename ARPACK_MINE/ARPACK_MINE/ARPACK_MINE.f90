    !| This is the modification version of the parallel arpack for F2003  |
    !| all the modifications are individualy developed by YC Xi           |
    !| at Tsinghua University For the purpose of using, only the          |
    !| complex version of double precision is considered                  |
    !| The whole structure of the lib files are divided into three main   |
    !| parts                                                              |
    !|    COMPLEX PARALLEL ARPACK:                                        |
    !|     TWO MAIN MODULE FOR DEBUGING AND TIME STAT                     |
    !|     DEBUG_MODULE.f90   STAT_MODULE.f90                             |
    !|                                                                    | 
    !|     1.  LOCAL ARPACK + LOCAL ARPACK UTIL                           |
    !|         The source files are in folder Local_arpack                |
    !|         The user define function are in folder Local_arpack_util   |
    !|                                                                    |
    !|     2.  PARALLEL ARPACK + PARAPACK UTIL                            |
    !|         The source files are in folder Paral_arpack                |
    !|         The user define function are in folder Paral_arpack_util   |
    !|                                                                    |
    !|     3.  EXAMPLES and USING                                         |
    !|         Based on previous version, six examples are designed for   |
    !|         illustrating                                               |
    !|         "1. For Simplest Case of Using ARPACK"                     |
    !|                                                                    | 
    !|         "2. For standard eigenvalue case "                         |
    !|                                                                    | 
    !|         "3. For standard eigenvalue case of Shift-invert mode "    |
    !|                                                                    |
    !|         "4. For generalized eigenvalue case "                      |
    !|                                                                    |
    !|         "5. For generalized eigenvalue case of Shift-invert mode " |
    !|                                                                    |
    !|         "6. For Parallel test "                                    |
    !|                                                                    |
    !=======================================================================

    PROGRAM MAIN_ARPACK
    IMPLICIT NONE
    include 'mpif.h'
    INTEGER( kind = 4 )::ierr,rc,myid
    
    CALL MPI_INIT( ierr )
    CALL MPI_COMM_RANK( MPI_COMM_WORLD, myid, ierr )
    
    if( myid == 0) then
    CALL LOCAL_ARPACK_TEST_EXAMPLES
    endif
    
    CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
    
    if( myid == 0) then
        WRITE(*,*)        
        WRITE(*,*)"Parallel test"
    end if
    
    CALL pzndrv1
    
    CALL MPI_FINALIZE(rc)
    
    END PROGRAM MAIN_ARPACK