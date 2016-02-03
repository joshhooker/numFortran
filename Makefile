SHELL = /bin/sh

# FORTRAN compiler
FC = gfortran-mp-5

FFLAGS = -O2 -std=f2008 -fno-range-check -ffree-form

TARGET = libraryTest
DATADIR = data
SRCDIR = src
OBJDIR = obj

MODULEFLAG = -J
VPATH = $(SRCDIR):$(OBJDIR)

# add any extra objects here
OBJFILES = libraryTest.o constants.o cubicSpline.o odeSolver.o \
           randomNumbers.o statTests.o integration.o specialFunctions.o \

FULLTARGET = $(TARGET)

# Rules to build the fortran files

.SUFFIXES: .f08 .o
.f08.o: ; @mkdir -p $(OBJDIR) $(DATADIR)
	$(FC) -c $(FFLAGS)  $(MODULEFLAG) $(OBJDIR) -o $(OBJDIR)/$@ $<

.SUFFIXES: .F08 .o
.F08.o: ; @mkdir -p $(OBJDIR) $(DATADIR)
	$(FC) -c $(FFLAGS) $(DFLAGS) $(MODULEFLAG) $(OBJDIR) -o $(OBJDIR)/$@ $<

%.o: %.mod

$(FULLTARGET): $(OBJFILES)
	$(FC) $(FFLAGS) $(MODULEFLAG) $(OBJDIR) -o $@ $(addprefix $(OBJDIR)/, $(OBJFILES))


.PHONEY: deepclean
deepclean:
	@rm -rf *~ $(BINDIR) $(OBJDIR) $(DATADIR) $(SRCDIR)/*~ *.out *.err *.log *.ipo

.PHONEY: clean
clean:
	@rm -rf *~ $(BINDIR) $(OBJDIR)  $(SRCDIR)/*~ *.out *.err *.log *.ipo

.PHONEY: tidy
tidy:
	@rm -rf $(OBJDIR) $(SRCDIR)/*~ *.out *.err *.log *.ipo

# Add Dependencies here
libraryTest.o : constants.o cubicSpline.o odeSolver.o randomNumbers.o \
								statTests.o integration.o specialFunctions.o
