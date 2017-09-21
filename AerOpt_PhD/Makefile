# Start of the makefile

# Defining variables
#ifneq (,$(filter Windows%,$(OS)))
ifeq ($(OS),Windows_NT)
EXT = .exe
else
EXT =
endif
TARGET = AerOpt$(EXT)
OBJECTS = Toolbox.o  InputData.o  ReadData.o  CreateSnapshots.o  FDGD.o Smoothing.o  GenerateMesh.o  CFD.o  Optimization.o  Main.o

# Find all sources
#SOURCES:=$(wildcard *.f90)
# Make targets out of them
#OBJECTS:=$(SOURCES:%.f90=%.o) #NOTE Order is important sometimes

# Set default output
default :
	make -s comment
comment :
	echo " "
	echo "   Makefile for mgns2d        "
	echo "   -------------------------- "
	echo "   clean  to remove object files "
	echo "   ifort  to compile and link with intel compiler"
	echo "   gfort  to compile and link with gfortran compiler"
	echo " "

ifort : COMP = ifort
ifort : FLAGS += -r8 -O3 -static -static-intel -mkl -I$(MKL_INC)/em64t/lp64 -I$(MKL_INC)/include
ifort : $(OBJECTS)
	$(COMP) $(FLAGS) $(OBJECTS) -o $(TARGET) -L$(MKL_LIB) -Wl, -mkl -lpthread -lm

gfort : COMP = gfortran
gfort : FLAGS = -fdefault-real-8 -O3 -static -static-libgfortran
gfort : $(OBJECTS)
	$(COMP) $(FLAGS) $(OBJECTS) -o $(TARGET)

# List of object files
%.o : %.f90
	$(COMP) $(FLAGS) -c $<

# Cleaning everything
clean :
	find .. -name '*.o' -delete
	find .. -name '*.mod' -delete
# End of the makefile
