SUBDIRS = mfl tui nui

all: $(SUBDIRS)

.PHONY: $(SUBDIRS) clean
MFLAGS += --no-print-directory
$(SUBDIRS):
	@echo Switching to $@
	@cd $@; $(MAKE) $(MFLAGS) 

clean:
	@for i in $(SUBDIRS); do \
	echo "Cleaning $$i..."; \
	(cd $$i; $(MAKE) $(MFLAGS) $(MYMAKEFLAGS) clean); done
