# FORTRAN compiler
FC = gfortran
FFLAGS = -O3 -std=f2008 -fno-range-check

#FC = ifort
#FFLAGS = -O3 -fast -parallel -ipo

TARGET = libraryTest
DATADIR = data
SRCDIR = src
OBJDIR = obj
TESTDIR = test

VPATH = $(SRCDIR):$(OBJDIR):$(TESTDIR)

# add any extra objects here
OBJFILES = libraryTest.o nfConstants.o nfIntegration.o
# OBJFILES = libraryTest.o constants.o cubicSpline.o odeSolver.o \
#            randomNumbers.o statTests.o integration.o specialFunctions.o \
#            sorting.o matrixSolver.o testRandNum.o \
#            nfConstants.o nfIntegration.o

FULLTARGET = $(TARGET)

# Rules to build the fortran files

.SUFFIXES: .f90 .o
.f90.o: ; @mkdir -p $(OBJDIR) $(DATADIR)
	$(FC) -c $(FFLAGS) $(OBJDIR) -o $(OBJDIR)/$@ $<

.SUFFIXES: .F90 .o
.F90.o: ; @mkdir -p $(OBJDIR) $(DATADIR)
	$(FC) -c $(FFLAGS) $(DFLAGS) $(OBJDIR) -o $(OBJDIR)/$@ $<

%.o: %.mod

$(FULLTARGET): $(OBJFILES)
	$(FC) $(FFLAGS) $(OBJDIR) -o $@ $(addprefix $(OBJDIR)/, $(OBJFILES))

.PHONEY: clean
clean:
	@rm -rf *~ $(BINDIR) $(OBJDIR)  $(SRCDIR)/*~ *.out *.err *.log *.ipo

.PHONEY: tidy
tidy:
	@rm -rf $(OBJDIR) $(SRCDIR)/*~ *.out *.err *.log *.ipo

# Add Dependencies here
libraryTest.o : nfConstants.o nfIntegration.o
# libraryTest.o : constants.o cubicSpline.o odeSolver.o randomNumbers.o \
# 								statTests.o integration.o specialFunctions.o sorting.o \
# 								matrixSolver.o testRandNum.o

nfIntegration.o : nfConstants.o

cubicSpline.o : matrixSolver.o

specialFunctions.o: constants.o integration.o odeSolver.o

statTests.o : randomNumbers.o sorting.o
