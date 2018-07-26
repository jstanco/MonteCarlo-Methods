IDIR = ../include
ODIR = ./obj
LDIR = ../lib
LPLDIR = /Users/johnstanco/PhysicsSoftware/lapack-3.8.0
BLASLDIR = /opt/OpenBlas/lib
ARMALDIR = /usr/local/lib
F90LDIR = /usr/local/gfortran/lib
FFTWLDIR = /usr/local/lib
BINDIR = ../bin
CC = g++
DEBUG = -g
CFLAGS = -c -I$(IDIR) $(DEBUG)

_OBJS = help.o ising_lattice.o
OBJS = $(patsubst %,$(ODIR)/%,$(_OBJS))
_DEPS = markov_chain.hpp markov_chain_mc.hpp ising_state.hpp ising_lattice.hpp pimc.hpp state_base.hpp part_storage.hpp signal.hpp help.hpp
DEPS = $(patsubst %,$(IDIR)/%,$(_DEPS))

LFLAGS = $(DEBUG)
LIBS = $(ARMALDIR)/libarmadillo.dylib $(BLASLDIR)/libopenblas.dylib $(FFTWLDIR)/libfftw3.a

all: ising ising_t_dependence

#PIMC : pimc

Ising : ising ising_t_dependence 

#pimc : $(ODIR)/test_pimc.o $(ODIR)/pimc_state.o $(OBJS)
#	@$(CC) -o $(BINDIR)/$@ $^ $(LFLAGS) $(LIBS)

ising : $(ODIR)/run_ising.o $(OBJS)
	@$(CC) -o $(BINDIR)/$@ $^ $(LFLAGS) $(LIBS)

ising_t_dependence : $(ODIR)/ising_temp_dependence.o $(OBJS)
	@$(CC) -o $(BINDIR)/$@ $^ $(LFLAGS) $(LIBS)

$(ODIR)/%.o : %.cpp $(DEPS)
	@$(CC) -o $@ $< $(CFLAGS)
	
.PHONY: clean

clean : 
	@rm -f $(ODIR)/*.o *~ core $(INCDIR)/*~
