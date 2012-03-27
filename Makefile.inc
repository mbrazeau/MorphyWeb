CC = g++
AR          = ar
ODIR        = obj
OBJECTS     = $(patsubst %,$(ODIR)/%,$(OBJS))
LIBRARIES  += $(patsubst %,-l%,$(LIBS)) $(patsubst %,-L%,$(LIBDIRS))
CFLAGS     += $(patsubst %,-I%,$(INCLUDEDIRS)) -Wall -g
MFLLIB      = $(ODIR)/libmfl.a
DEPS        = $(patsubst %,$(ODIR)/%,$(OBJS:.o=.d))
MAKEDEPEND  = $(SILENT)$(CC) -MF"$(@:.o=.d)" -MG -MM -MP -MT"$(@:.d=.o)" "$<" $(CFLAGS) 
MKODIR      = -@mkdir -p $(ODIR)

ifndef VERBOSE
SILENT=@
NICEDEPOUTPUT=@echo "  Dependency $(notdir $(CURDIR))/$@"
NICEBLDOUTPUT=@echo "  Building $(notdir $(CURDIR))/$@"
endif

ifdef DEBUG
CFLAGS += -DMFY_DEBUG
endif

all: $(OUTEXE) $(OUTLIB)

$(ODIR)/%.d: %.cpp 
	$(NICEDEPOUTPUT)
	$(MKODIR)
	$(MAKEDEPEND)

$(ODIR)/%.d: %.c
	$(NICEDEPOUTPUT)
	$(MKODIR)
	$(MAKEDEPEND)

$(ODIR)/%.o: %.cpp 
	$(NICEBLDOUTPUT)
	$(MKODIR)
	$(SILENT)$(CC) -c -o $@ $< $(CFLAGS)

$(ODIR)/%.o: %.c 
	$(NICEBLDOUTPUT)
	$(MKODIR)
	$(SILENT)$(CC) -c -o $@ $< $(CFLAGS)

$(OUTLIB): $(OBJECTS)
	$(NICEBLDOUTPUT)
	$(MKODIR)
	$(SILENT)$(AR) rcs $@ $^

$(OUTEXE): $(OBJECTS)
	$(NICEBLDOUTPUT)
	$(MKODIR)
	$(SILENT)$(CC) -o $@ $^ $(CFLAGS) $(LIBRARIES) -Wl,-rpath -Wl,/usr/local/lib/ncl

.PHONY: clean

clean:
	$(SILENT)rm -f $(ODIR)/*.o $(ODIR)/*.d *~ core $(INCDIR)/*~ $(OUTEXE) $(OUTLIB) 


