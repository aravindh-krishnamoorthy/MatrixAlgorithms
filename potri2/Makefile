.PHONY: all

SHLIB := potri2.so
SRCS  := dpotri2.f90 dpotri2s.f90 dpotri2b.f90
FOPTS += -fexternal-blas -Ofast -march=native -masm=intel -fverbose-asm -fPIC -shared -save-temps
FLIBS := -lgfortran -lm
FLIBS += -L/opt/intel/mkl/lib/intel64/ -lmkl_rt

all:
	gfortran $(FOPTS) -o $(SHLIB) $(SRCS) -L/usr/local/lib $(FLIBS)

