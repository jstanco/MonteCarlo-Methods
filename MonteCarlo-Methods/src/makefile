IDIR = ../include
ODIR = ./obj
LDIR = ../lib
BINDIR = ../bin
CC = g++
DEBUG = -g
CFLAGS = -c -I$(IDIR) $(DEBUG)

_OBJS = Ising_1D.o Ising_2D.o help.o testIsingChain.o
OBJS = $(patsubst %,$(ODIR)/%,$(_OBJS))
_DEPS = MarkovChain.hpp Metropolis_Hastings_MC.hpp Distribution.hpp Ising_Base.hpp Ising_1D.hpp Ising_2D.hpp help.hpp
DEPS = $(patsubst %,$(IDIR)/%,$(_DEPS))

LFLAGS = $(DEBUG)
LIBS = 

all: main 

Markov : markov

markov : testMarkov.o
	@$(CC) -o $(BINDIR)/$@ $^ $(LFLAGS) $(LIBS)

main : $(OBJS)
	@$(CC) -o $(BINDIR)/$@ $^ $(LFLAGS) $(LIBS)

testIsingChain.o : %.cpp $(DEPS)
	@$(CC) -o $@ $< $(CFLAGS)

$(ODIR)/%.o : %.cpp $(DEPS)
	@$(CC) -o $@ $< $(CFLAGS) 
	
.PHONY: clean

clean : 
	@rm -f $(ODIR)/*.o *~ core $(INCDIR)/*~
