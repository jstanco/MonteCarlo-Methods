IDIR = ../include
ODIR = ./obj
LDIR = ../lib
BINDIR = ../bin
CC = g++
DEBUG = -g
CFLAGS = -c -I$(IDIR) $(DEBUG)

_OBJS = help.o
OBJS = $(patsubst %,$(ODIR)/%,$(_OBJS))
_DEPS = markov_chain.hpp metropolis_hastings_mc.hpp ising_state.hpp ising_lattice.hpp pimc.hpp state_base.hpp help.hpp
DEPS = $(patsubst %,$(IDIR)/%,$(_DEPS))

LFLAGS = $(DEBUG)
LIBS = /usr/local/lib/libarmadillo.dylib

all: main 

Markov : markov

PIMC : pimc

Ising: ising

markov : $(ODIR)/test_markov.o
	@$(CC) -o $(BINDIR)/$@ $^ $(LFLAGS) $(LIBS)

pimc : $(ODIR)/test_pimc.o $(ODIR)/pimc_state.o $(OBJS)
	@$(CC) -o $(BINDIR)/$@ $^ $(LFLAGS) $(LIBS)

ising : $(ODIR)/run_ising.o $(OBJS)
	@$(CC) -o $(BINDIR)/$@ $^ $(LFLAGS) $(LIBS)

$(ODIR)/%.o : %.cpp $(DEPS)
	@$(CC) -o $@ $< $(CFLAGS) 
	
.PHONY: clean

clean : 
	@rm -f $(ODIR)/*.o *~ core $(INCDIR)/*~
