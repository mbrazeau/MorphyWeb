OUTEXE=$(MFLLIB)
OBJS = rearrange.o taxpart.o tree.o readnewick.o
include ../Makefile.inc

$(MFLLIB): $(OBJECTS)
	$(SILENT)$(AR) rcs $@ $^
	$(NICEOUTPUT)

