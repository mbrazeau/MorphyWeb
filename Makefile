IDIR =.
CC=gcc
CFLAGS=-I$(IDIR) -Wall

ODIR=obj
LDIR =.

LIBS=-lm
EXE=morphy

_DEPS = morphy.h
DEPS = $(patsubst %,$(IDIR)/%,$(_DEPS))

_OBJ = main.o randtree.o taxpart.o
OBJ = $(patsubst %,$(ODIR)/%,$(_OBJ))

all: $(EXE)

$(ODIR)/%.o: %.c $(DEPS)
	-@mkdir -p obj
	$(CC) -c -o $@ $< $(CFLAGS) -g

$(EXE): $(OBJ)
	gcc -o $@ $^ $(CFLAGS) $(LIBS)

.PHONY: clean

clean:
	rm -f $(ODIR)/*.o *~ core $(INCDIR)/*~ $(EXE)

