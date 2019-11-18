MPI=true

ifeq ($(MPI),true)
   FC=mpif90
   MPIFILE=mympi
   EXEC=dmcbcsmpi
else
   FC=ifort
   MPIFILE=nompi
   EXEC=dmcbcs
endif

LDFLAGS   = -L$(MKLROOT)/lib/intel64 -I$(MKLROOT)/include
LIBS=-lmkl_intel_lp64 -lmkl_core -lmkl_intel_thread -openmp -check all -warn all

FFLAGS=-O3 #-pg #-prof-use
#FFLAGS=-check all -O0 -warn all

#LDFLAGS=-L$(HOME)/home/lib
#LIBS=-llapackmin -lblasmin -check all -warn all #-pg #-prof-use

POTEN=coshpot
#POTEN=doublegauss
#POTEN=coshpotbox
#POTEN=aziz
#POTEN=toypot
#POTEN=tabpot
#POTEN=v1pot
#POTEN=vbarrier

BACKFLOW=backflow
#BACKFLOW=backflowhe

MATRIX=matrixlapack

.SUFFIXES: .f90

.f90.o:
	$(FC) $(FFLAGS) -c $<

OBJECTS=\
   dmcbcs.o\
   estimator.o\
   lattice.o\
   ran.o\
   stack.o\
   step.o\
   wavefunction.o\
   $(MATRIX).o\
   $(MPIFILE).o\
   kshell.o\
   gofr.o\
   betafun.o\
   psipfaffian.o\
   optimizer.o\
   psibose.o\
   psislater.o\
   psislaterbcs.o\
   psislaterbcsloff.o\
   psibcstrapiso_xyz.o\
   pfaffianmod.o\
   sofq.o\
   psiatomtab.o\
   density.o\
   radii.o\
   operators.o\
   momentum3d.o\
   $(BACKFLOW).o\
   ylm.o\
   vext.o\
   euclidean.o\
   debug.o\
   jastrowlongrange.o

dmcbcs: $(OBJECTS)
	$(FC) $(LDFLAGS) -o ./$(EXEC) $(OBJECTS) $(LIBS) $(OBJPOT)

include makefile.$(POTEN)

clean:
	rm -f *\.o *\.mod *~ $(EXEC)

dmcbcs.o: stack.o lattice.o ran.o wavefunction.o step.o estimator.o\
   $(MPIFILE).o $(POTEN).o dmcbcs.f90 optimizer.o operators.o
	$(FC) $(FFLAGS) -c dmcbcs.f90

estimator.o: $(MPIFILE).o estimator.f90
	$(FC) $(FFLAGS) -c estimator.f90

lattice.o: lattice.f90
	$(FC) $(FFLAGS) -c lattice.f90

ran.o: ran.f90
	$(FC) $(FFLAGS) -c ran.f90

stack.o: stack.f90
	$(FC) $(FFLAGS) -c stack.f90

step.o: stack.o ran.o wavefunction.o estimator.o debug.o step.f90
	$(FC) $(FFLAGS) -c step.f90

wavefunction.o: $(POTEN).o $(MATRIX).o $(BACKFLOW).o stack.o psibose.o jastrowlongrange.o psipfaffian.o\
   psislater.o psislaterbcs.o psislaterbcsloff.o psiatomtab.o psibcstrapiso_xyz.o wavefunction.f90
	$(FC) $(FFLAGS) -c wavefunction.f90

psibose.o: $(MPIFILE).o psibose.f90
	$(FC) $(FFLAGS) -c psibose.f90

psislater.o: $(MPIFILE).o $(MATRIX).o $(BACKFLOW).o psislater.f90
	$(FC) $(FFLAGS) -c psislater.f90

psislaterbcs.o: $(MPIFILE).o $(MATRIX).o kshell.o betafun.o psislaterbcs.f90
	$(FC) $(FFLAGS) -c psislaterbcs.f90

psislaterbcsloff.o: $(MPIFILE).o $(MATRIX).o kshell.o betafun.o psislaterbcsloff.f90
	$(FC) $(FFLAGS) -c psislaterbcsloff.f90

psipfaffian.o: $(MPIFILE).o $(MATRIX).o betafun.o kshell.o pfaffianmod.o psipfaffian.f90
	$(FC) $(FFLAGS) -c psipfaffian.f90

pfaffianmod.o: pfaffianmod.f90
	$(FC) $(FFLAGS) -c pfaffianmod.f90

psiatomtab.o: $(MPIFILE).o $(MATRIX).o psiatomtab.f90
	$(FC) $(FFLAGS) -c psiatomtab.f90

psibcstrapiso_xyz.o: $(MPIFILE).o $(MATRIX).o ylm.o psibcstrapiso_xyz.f90
	$(FC) $(FFLAGS) -c psibcstrapiso_xyz.f90

$(MATRIX).o: $(MATRIX).f90
	$(FC) $(FFLAGS) -c $(MATRIX).f90

$(MPIFILE).o: $(MPIFILE).f90
	$(FC) $(FFLAGS) -c $(MPIFILE).f90

kshell.o: kshell.f90
	$(FC) $(FFLAGS) -c kshell.f90

gofr.o: $(MPIFILE).o gofr.f90
	$(FC) $(FFLAGS) -c gofr.f90

betafun.o: betafun.f90
	$(FC) $(FFLAGS) -c betafun.f90

optimizer.o: $(MATRIX).o wavefunction.o optimizer.f90
	$(FC) $(FFLAGS) -c optimizer.f90

operators.o: gofr.o sofq.o density.o radii.o momentum3d.o vext.o euclidean.o operators.f90
	$(FC) $(FFLAGS) -c operators.f90

sofq.o: kshell.o $(MPIFILE).o sofq.f90
	$(FC) $(FFLAGS) -c sofq.f90

radii.o: $(MPIFILE).o radii.f90
	$(FC) $(FFLAGS) -c radii.f90

density.o: $(MPIFILE).o density.f90
	$(FC) $(FFLAGS) -c density.f90

momentum3d.o: $(MPIFILE).o momentum3d.f90 stack.o ran.o wavefunction.o
	$(FC) $(FFLAGS) -c momentum3d.f90

$(BACKFLOW).o: $(BACKFLOW).f90
	$(FC) $(FFLAGS) -c $(BACKFLOW).f90

ylm.o: ylm.f90
	$(FC) $(FFLAGS) -c ylm.f90

vext.o: wavefunction.o vext.f90
	$(FC) $(FFLAGS) -c vext.f90

euclidean.o: kshell.o $(MPIFILE).o euclidean.f90
	$(FC) $(FFLAGS) -c euclidean.f90

debug.o: stack.o debug.f90
	$(FC) $(FFLAGS) -c debug.f90

jastrowlongrange.o: kshell.o jastrowlongrange.f90
	$(FC) $(FFLAGS) -c jastrowlongrange.f90
