PROGNAME=simulation

IDIR =./src/include
ODIR=./src/obj
SRCDIR=./src

CC=g++
C=gcc

CFLAGST=-I$(IDIR) -I./src/CLI11/include/ `pkg-config --cflags gsl` -ggdb -fexceptions -Wall -pg -no-pie # for testing
CFLAGS=-I$(IDIR) -I./src/CLI11/include/ `pkg-config --cflags gsl` -Ofast # for stuff

LIBS= -lboost_system -lm `pkg-config --libs gsl` -fopenmp

_DEPS = randomgen.h model.h
DEPS = $(patsubst %,$(IDIR)/%,$(_DEPS))

_OBJ = diffe.o model.o randomgen.o
OBJ = $(patsubst %,$(ODIR)/%,$(_OBJ))


$(ODIR)/%.o: $(SRCDIR)/%.cpp $(DEPS)
	@mkdir -p ${ODIR}
	$(CC) -c -o $@ $< $(CFLAGS)

$(ODIR)/%.o: $(SRCDIR)/%.c $(DEPS)
	@mkdir -p ${ODIR}
	$(C) -c -o $@ $< $(CFLAGS)

$(PROGNAME): $(OBJ)
	$(CC) -o $@ $^ $(CFLAGS) $(LIBS)

.PHONY: gdb
gdb: debug
gdb: CFLAGS=$(CFLAGST)
debug: $(OBJ)
	$(CC) -o $@ $^ $(CFLAGST) $(LIBS)

.PHONY: clean

clean:
	rm -f $(ODIR)/*.o *~ core $(INCDIR)/*~  

