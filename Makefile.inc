CC   = gcc
CPP  = g++
AR   = ar
ODIR = obj
OBJECTS     = $(patsubst %,$(ODIR)/%,$(OBJS))
LIBRARIES  += $(patsubst %,-l%,$(LIBS)) $(patsubst %,-L%,$(LIBDIRS))
CFLAGS     += $(patsubst %,-I%,$(INCLUDEDIRS)) -Wall -g
MFLLIB      = $(ODIR)/libmfl.a
ifndef VERBOSE
SILENT=@
NICEOUTPUT=@echo Building $(notdir $(CURDIR))/$@
endif

all: $(OUTEXE)

$(ODIR)/%.o: %.cpp $(DEPS)
	-@mkdir -p $(ODIR)
	$(SILENT)$(CPP) -c -o $@ $< $(CFLAGS)
	$(NICEOUTPUT)

$(ODIR)/%.o: %.c $(DEPS)
	-@mkdir -p $(ODIR)
	$(SILENT)$(CC) -c -o $@ $< $(CFLAGS)
	$(NICEOUTPUT)

.PHONY: clean

clean:
	$(SILENT)rm -f $(ODIR)/*.o *~ core $(INCDIR)/*~ $(OUTEXE)


