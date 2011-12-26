SUBDIRS = mfl tui

all: $(SUBDIRS)

.PHONY: $(SUBDIRS) clean

$(SUBDIRS):
	cd $@; $(MAKE) $(MFLAGS)

clean:
	@for i in $(SUBDIRS); do \
	echo "Clearing in $$i..."; \
	(cd $$i; $(MAKE) $(MFLAGS) $(MYMAKEFLAGS) clean); done
