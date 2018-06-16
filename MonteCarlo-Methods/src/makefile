IDIR = ../include
ODIR = ./obj
LDIR = ../lib
BINDIR = ../bin
CC = g++
DEBUG = -g
CFLAGS = -c -I$(IDIR) $(DEBUG)

_OBJS = ising_1d.o ising_2d.o pimc_state.o help.o lattice.o
OBJS = $(patsubst %,$(ODIR)/%,$(_OBJS))
_DEPS = markov_chain.hpp metropolis_hastings_mc.hpp ising_base.hpp ising_1d.hpp ising_2d.hpp pimc.hpp state_base.hpp help.hpp lattice.hpp
DEPS = $(patsubst %,$(IDIR)/%,$(_DEPS))

LFLAGS = $(DEBUG)
LIBS = /usr/local/lib/libarmadillo.dylib

all: main 

Markov : markov

PIMC : pimc

markov : $(ODIR)/test_markov.o
	@$(CC) -o $(BINDIR)/$@ $^ $(LFLAGS) $(LIBS)

pimc : $(ODIR)/test_pimc.o $(OBJS)
	@$(CC) -o $(BINDIR)/$@ $^ $(LFLAGS) $(LIBS)

main : $(ODIR)/test_ising_chain.o $(OBJS)
	@$(CC) -o $(BINDIR)/$@ $^ $(LFLAGS) $(LIBS)

$(ODIR)/%.o : %.cpp $(DEPS)
	@$(CC) -o $@ $< $(CFLAGS) 
	
.PHONY: clean

clean : 
	@rm -f $(ODIR)/*.o *~ core $(INCDIR)/*~
