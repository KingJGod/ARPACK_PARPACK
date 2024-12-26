#COMPILE THE FILE
FC := mpifort
AR := ar
FCFLAGS = -c -O2 -fPIC
FLFLAGS =
#MKLROOT = /opt/intel/mkl 
#INC = 
#LIBS = -L/Users/xiyoucheng/Documents/My_Works/StabilityWork/Global_Linear_Stability_Tools/CLSTOF/biglobal_tools/blaslapack -lfblas -lflapack

MKLROOT = /opt/intel/mkl 
INC = -I/opt/intel/mkl/include
LIBS = -L/opt/intel/mkl/lib/intel64 \
		-lmkl_intel_lp64 -lmkl_intel_thread \
		-lmkl_core -liomp5 -lpthread -lm -ldl -lc++

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
