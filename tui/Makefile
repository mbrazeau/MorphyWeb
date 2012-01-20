OUTEXE=$(ODIR)/morphy
INCLUDEDIRS =. ../mfl
OBJS = main.o randtree.o

include ../Makefile.inc

$(OUTEXE): $(OBJECTS)
	$(SILENT)$(CC) -o $@ $^ $(CFLAGS) $(LIBRARIES) ../mfl/$(MFLLIB)
	$(NICEOUTPUT)
    
