NUIEXE=$(ODIR)/nui
INCLUDEDIRS=. 
LIBDIRS=/usr/local/lib/ncl
_LIBS=ncl
LIBS=$(patsubst %,-l%,$(_LIBS)) $(patsubst %,-L%,$(LIBDIRS))

include ../Makefile.inc

_OBJ = nui.o
OBJ = $(patsubst %,$(ODIR)/%,$(_OBJ))

all: $(NUIEXE)

$(ODIR)/%.o: %.cpp $(DEPS)
	-@mkdir -p $(ODIR)
	$(CPP) -c -o $@ $< $(CFLAGS) -g

$(NUIEXE): $(OBJ)
	$(CPP) -o $@ $^ $(CFLAGS) $(LIBS) -g

.PHONY: clean

clean:
	rm -f $(ODIR)/*.o *~ core $(INCDIR)/*~ $(NUIEXE)

