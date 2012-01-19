include ../Makefile.inc

_OBJ = rearrange.o taxpart.o tree.o readnewick.o
OBJ = $(patsubst %,$(ODIR)/%,$(_OBJ))

all: $(MFLLIB)

$(ODIR)/%.o: %.c $(DEPS)
	-@mkdir -p $(ODIR)
	$(CC) -c -o $@ $< $(CFLAGS) -g

$(MFLLIB): $(OBJ)
	$(AR) rcs $@ $^

.PHONY: clean

clean:
	rm -f $(ODIR)/*.o *~ core $(INCDIR)/*~ $(MFLLIB)
