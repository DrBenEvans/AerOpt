#Compilers and flags

MF	= ifort
MC	= cc
ML	= ar
RM	= rm
MKL_INC = ${MKLROOT}/include
MKL_LIB = ${MKLROOT}/lib/intel64

FFLAGS	=  -r8 -O3 -mkl=sequential -I$(MKL_INC)/intel64/lp64 
CFLAGS	=  -r8 -O3 
#FFLAGS	= -g 
#CFLAGS	= -g 
LFLAGS	= crv

#C and Fortran Object Files
FOBJS	=  Toolbox.o  InputData.o  ReadData.o  CreateSnapshots.o FDGD.o  Smoothing.o  GenerateMesh.o  CFD.o  Optimization.o  Main.o

OBJS	=  $(FOBJS) 

#Compile and link

default:
	make -s comment

all:  $(OBJS)
	$(MF) $(FFLAGS) -o AerOpt_v3.5 $(OBJS) -L$(MKL_LIB) \
		$(MKL_LIB)/libmkl_blas95_lp64.a -Wl,--start-group $(MKL_LIB)/libmkl_intel_lp64.a $(MKL_LIB)/libmkl_core.a $(MKL_LIB)/libmkl_sequential.a -Wl,--end-group -lpthread -lm

comment:
		echo " "
		echo "   Makefile for AerOpt       "
		echo "   ------------------------- "
		echo "    clean  to remove object files "              
		echo "    all    to compile and link "
		echo " "
		
clean:
	find .. -name '*.o' -delete

Toolbox.o:
	$(MF) $(FFLAGS) -c Toolbox.f90
InputData.o:
	$(MF) $(FFLAGS) -c InputData.f90
ReadData.o:
	$(MF) $(FFLAGS) -c ReadData.f90
CreateSnapshots.o:
	$(MF) $(FFLAGS) -c CreateSnapshots.f90
FDGD.o:
	$(MF) $(FFLAGS) -c FDGD.f90
Smoothing.o:
	$(MF) $(FFLAGS) -c Smoothing.f90
GenerateMesh.o:
	$(MF) $(FFLAGS) -c GenerateMesh.f90
CFD.o:
	$(MF) $(FFLAGS) -c CFD.f90
Optimization.o:
	$(MF) $(FFLAGS) -c Optimization.f90
Main.o:
	$(MF) $(FFLAGS) -c Main.f90
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::#



