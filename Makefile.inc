CC          = g++
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

$(OUTEXE): $(OBJECTS) $(ARCHIVES)
	$(NICEBLDOUTPUT)
	$(MKODIR)
	$(SILENT)$(CC) -o $@ $^ $(CFLAGS) $(LIBRARIES) $(ARCHIVES) -Wl,-rpath -Wl,/usr/local/lib/ncl

.PHONY: clean

clean:
	$(SILENT)rm -f $(ODIR)/*.o $(ODIR)/*.d *~ core $(INCDIR)/*~ $(OUTEXE) $(OUTLIB)

#
# How to add new tests:
#
# Step 1:
# Create a new directory under the tests directory that matches the
# wildcard pattern: test*
#
# Step 2:
# In that new directory create 2 files, one file should match this 
# wildcard pattern: *input.txt
# and the other should match this wildcard pattern: *output.txt
# You can add other files to each directory so long as they do not
# match one of the wildcards, then they will be ignored.
#
# Step 3:
# In the input file add the commands you wish to run during the test
# These can be captured with the LOG command in NUI if you like.
#
# Step 4:
# In the output file put the contents of the output from NUI when the
# input commands are run.
#
# Step 5:
# Run make test from the main Morphy directory
#
# Step 6:
# If the output file does not exactly match the output NUI gave during
# the make test run, then the test failed and you need to fix it!
#
$(TESTDIRS): 
	$(eval TMP_IN := $(wildcard $@/*input.txt))
	$(eval TMP_OUT := $(wildcard $@/*output.txt))
	$(SILENT)./runtest.pl --test=$@ --input="$(TMP_IN)" --expected="$(TMP_OUT)"


