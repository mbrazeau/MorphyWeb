TUIEXE=$(ODIR)/morphy
INCLUDEDIRS =. ../mfl
include ../Makefile.inc

_OBJ = main.o randtree.o
OBJ = $(patsubst %,$(ODIR)/%,$(_OBJ))

all: $(TUIEXE)

$(ODIR)/%.o: %.c $(DEPS)
	-@mkdir -p $(ODIR)
	$(CC) -c -o $@ $< $(CFLAGS) -g

$(TUIEXE): $(OBJ)
	gcc -o $@ $^ $(CFLAGS) $(LIBS) ../mfl/$(MFLLIB)

.PHONY: clean

clean:
	rm -f $(ODIR)/*.o *~ core $(INCDIR)/*~ $(TUIEXE)

