SUBDIRS = mfl tui nui
ALLFILES = $(addsuffix /*.c, $(SUBDIRS)) $(addsuffix /*.cpp, $(SUBDIRS))
CSCOPEFILE = cscope.out
ifndef VERBOSE
SILENT=@
NICECSCOUTPUT=@echo " CScope for $(SUBDIRS)"
endif

all: $(SUBDIRS)

.PHONY: $(SUBDIRS) clean cscope

MFLAGS += --no-print-directory

$(SUBDIRS):
	@echo "Building $@..."
	@cd $@; $(MAKE) $(MFLAGS)

cscope: $(CSCOPEFILE)

$(CSCOPEFILE): $(wildcard $(ALLFILES))
	$(NICECSCOUTPUT)
	$(SILENT)cscope -b -u $(patsubst %,-s %,$(SUBDIRS)) -s ../ncl-2.1.17/ncl
    
clean:
	$(SILENT)rm -f $(CSCOPEFILE)
	@for i in $(SUBDIRS); do \
	echo "Cleaning $$i..."; \
	(cd $$i; $(MAKE) $(MFLAGS) $(MYMAKEFLAGS) clean); done
