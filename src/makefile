# Makefile for xCFD1D code

# Platform compile options 
include makefile.def

# Define the executable
EXEC = xCFD1D

# Create the modules directory
$(shell [ -d ${MODD} ] || mkdir -p ${MODD})

# List of objects
OBJS =	mathFuncs.o \
	cfdParams.o \
	gasConstants.o \
	stateVector.o \
	inputParams.o \
	gridBlock.o \
	fluxVector.o \
	exactSolns.o \
	solnBlock.o \
	gradRecon.o \
	explSolver.o \
	epsPlots.o \
	plotData.o

.SUFFIXES : .f90 .o

$(EXEC) :
	@echo Building $(EXEC) for a $(PLATFORM) platform.
	make $(OBJS)
	$(FC) $(LDFLAGS) $(OBJS) $(EXEC).f90 -o $(EXEC)
	@mv $(EXEC) ${BIND}
	@echo ' '

$(OBJS) :

.f90.o :
	$(FC) $(FFLAGS) -c $*.f90

clean:
	@echo Cleaning directory.
	rm -f *.o *.mod *.bak *~
	@echo ' '

clobber:
	@echo Clobbering directory.
	rm -f *.o *.mod *.bak *~
	rm -rf ${MODD} $(BIND)/${EXEC}
	@echo ' '
