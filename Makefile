SUBDIRS = Debug Release
SRCDIRS = mfl nui
CSCOPEFILE = cscope.out
ifndef VERBOSE
SILENT=@
NICECSCOUTPUT=@echo " CScope for $(SRCDIRS)"
endif

all: $(SUBDIRS)
	$(NICECSCOUTPUT)
	$(SILENT)cscope -b -u $(patsubst %,-s %,$(SRCDIRS)) -s ../ncl

.PHONY: $(SUBDIRS) clean realclean test

MFLAGS += --no-print-directory

$(SUBDIRS):
	@echo "Building $@..."
	-@mkdir $@
	@cd $@; cmake -DCMAKE_BUILD_TYPE=$@ ..
	@$(MAKE) $(MFLAGS) -C $@

cleancscope:
	$(SILENT)rm -f $(CSCOPEFILE)

clean: cleancscope
	@for i in $(SUBDIRS); do \
	echo "Cleaning $$i..."; \
	(cd $$i; $(MAKE) $(MFLAGS) clean); done

realclean: cleancscope
	@for i in $(SUBDIRS); do \
	echo "Deleting $$i..."; \
	rm -rf $$i; done

test test%: 
	@for i in $(SRCDIRS); do \
	echo "Testing $$i..."; \
	(cd $$i/tests; $(MAKE) $(MFLAGS) $@); done

