PROGNAME=simulation

IDIR =./src/include
ODIR=./src/obj
SRCDIR=./src

CC=g++
C=gcc

CFLAGS=-I$(IDIR) `pkg-config --cflags gsl` -ggdb -fexceptions -Wall -pg -no-pie # for testing
#CFLAGS=-I$(IDIR) `pkg-config --cflags gsl` -Ofast # for stuff

LIBS= -lboost_system -lm `pkg-config --libs gsl` -fopenmp

_DEPS = randomgen.h
DEPS = $(patsubst %,$(IDIR)/%,$(_DEPS))

_OBJ = diffe.o randomgen.o
OBJ = $(patsubst %,$(ODIR)/%,$(_OBJ))


$(ODIR)/%.o: $(SRCDIR)/%.cpp $(DEPS)
	@mkdir -p ${ODIR}
	$(CC) -c -o $@ $< $(CFLAGS)

$(ODIR)/%.o: $(SRCDIR)/%.c $(DEPS)
	$(C) -c -o $@ $< $(CFLAGS)

$(PROGNAME): $(OBJ)
	$(CC) -o $@ $^ $(CFLAGS) $(LIBS)

.PHONY: clean

clean:
	rm -f $(ODIR)/*.o *~ core $(INCDIR)/*~  

