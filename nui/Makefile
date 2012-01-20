OUTEXE=$(ODIR)/nui
INCLUDEDIRS=. 
LIBDIRS=/usr/local/lib/ncl
LIBS=ncl
OBJS = nui.o

include ../Makefile.inc

$(OUTEXE): $(OBJECTS)
	$(SILENT)$(CC) -o $@ $^ $(CFLAGS) $(LIBRARIES) 
	$(NICEOUTPUT)

