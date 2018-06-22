#COMPILE THE FILE
FC := mpiifort
AR := ar
FCFLAGS = -c -xHost -O2 -fPIC
FLFLAGS =
#MKLROOT = /opt/intel/mkl 
INC = -I$(MKLROOT)/include/intel64/lp64 -I$(MKLROOT)/include
LIBS = $(MKLROOT)/lib/intel64/libmkl_blas95_lp64.a \
       $(MKLROOT)/lib/intel64/libmkl_lapack95_lp64.a -L$(MKLROOT)/lib/intel64 \
		-lmkl_intel_lp64 -lmkl_intel_thread \
		-lmkl_core -liomp5 -lpthread -lm -ldl -lstdc++

#DIRs
#DIRs = Local_arpack Local_arpack_util Paral_arpack Paral_arpack_util
ARPACK_SRC = $(wildcard Local_arpack/*.f90)\
	$(wildcard Paral_arpack/*.f90)
ARPACK_OBJ = $(patsubst %.f90, %.o, $(ARPACK_SRC))

ARPACKF_SRC = $(wildcard Local_arpack_util/*.f)\
	$(wildcard Paral_arpack_util/*.f)
ARPACKF_OBJ = $(patsubst %.f, %.o, $(ARPACKF_SRC))

ARPACKEX_SRC = $(wildcard Example/*.f90)
ARPACKEX_OBJ = $(patsubst %.f90, %.o, $(ARPACKEX_SRC))

#PROG NAME
PROGRAM = ARPACK_TEST
TARGET = ./Lib/myparpack.a

all: $(PROGRAM) $(TARGET)

$(PROGRAM): $(ARPACK_OBJ) $(ARPACKF_OBJ) $(ARPACKEX_OBJ) ARPACK_MINE.o DEBUG_MODULE.o PCONTEXT_MODULE.o STAT_MODULE.o
	$(FC) $(INC) $(FLFLAGS) -o $@ $^ $(LIBS) 

$(TARGET): $(ARPACK_OBJ) $(ARPACKF_OBJ) DEBUG_MODULE.o PCONTEXT_MODULE.o STAT_MODULE.o
	$(AR)  cr $@ $^

$(ARPACK_OBJ): DEBUG_MODULE.o PCONTEXT_MODULE.o STAT_MODULE.o

$(ARPACKEX_OBJ): DEBUG_MODULE.o PCONTEXT_MODULE.o STAT_MODULE.o

debug_module.mod: DEBUG_MODULE.o

pcontext_module.mod: PCONTEXT_MODULE.o

stat_module.mod: STAT_MODULE.o

%.o: %.f90
	$(FC) $(FCFLAGS) -o $@ $<

%.o: %.f
	$(FC) $(FCFLAGS) -o $@ $<

clean:
	rm -rf *.o */*.o *.mod $(PROGRAM)
