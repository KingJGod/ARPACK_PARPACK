# ARPACK_PARPACK
For the purpose of calculation very large size eigenvalue problems. 
The complex arpack_parpack packages for double precision has been rewritten in fortran 90 standard.
Various test files are being pasted in system in Windows\Linux\Mac OS

    !| This is the modification version of the parallel arpack for F2003  |
    !| all the modifications are individualy developed by YC Xi           |
    !| at Tsinghua University For the purpose of using, only the          |
    !| complex version of double precision is considered                  |
    !| The whole structure of the lib files are divided into three main   |
    !| parts                                                              |
    !|    COMPLEX PARALLEL ARPACK:                                        |
    !|     TWO MAIN MODULE FOR DEBUGING AND TIME STAT                     |
    !|     DEBUG_MODULE.f90   STAT_MODULE.f90   PCONTEXT_MODULE.f90       |
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

BLAS LAPACK and MPI Libraries are needed.
In windows environment, VS2012 + Intel Visual Fortran 2017 with Intel MPI have been used.
And in Mac OS and Linux System. Both Intel MPI and MPICH have been tested. 
Intel MKL Library has been used to provide the basic environment for using package. You can modified the linking package by yourself.

COMPILING THE CODE IN MAC OS/Linux Systems:
"cd ARPACK_MINE/ARPACK_MINE"
modified the Makefile
" make "
And then you could get the excutable file "ARPACK_TEST" and lib "myparpack.a" in the Lib file folder

COMPILING THE CODE IN WINDOWS SYSTEM:
You can open the ARPACK_MINE.vfproj project file and run the project after linking the MPI and MKL library.


