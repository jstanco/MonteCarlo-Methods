IDIR = ../include
ODIR = ./obj
LDIR = ../lib
BINDIR = ../bin
CC = g++
DEBUG = -g
CFLAGS = -c -I$(IDIR) $(DEBUG)

_OBJS = Ising_1D.o Ising_2D.o PIMC_State.o  help.o lattice.o
OBJS = $(patsubst %,$(ODIR)/%,$(_OBJS))
_DEPS = MarkovChain.hpp Metropolis_Hastings_MC.hpp Ising_Base.hpp Ising_1D.hpp Ising_2D.hpp PIMC_Base.hpp PIMC_State.hpp State_Base_Bones.hpp help.hpp lattice.hpp
DEPS = $(patsubst %,$(IDIR)/%,$(_DEPS))

LFLAGS = $(DEBUG)
LIBS = /usr/local/lib/libarmadillo.dylib

all: main 

Markov : markov

PIMC : pimc

markov : $(ODIR)/testMarkov.o
	@$(CC) -o $(BINDIR)/$@ $^ $(LFLAGS) $(LIBS)

pimc : $(ODIR)/testPIMC.o $(OBJS)
	@$(CC) -o $(BINDIR)/$@ $^ $(LFLAGS) $(LIBS)

main : $(ODIR)/testIsingChain.o $(OBJS)
	@$(CC) -o $(BINDIR)/$@ $^ $(LFLAGS) $(LIBS)

$(ODIR)/%.o : %.cpp $(DEPS)
	@$(CC) -o $@ $< $(CFLAGS) 
	
.PHONY: clean

clean : 
	@rm -f $(ODIR)/*.o *~ core $(INCDIR)/*~
