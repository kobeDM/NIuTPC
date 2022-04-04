SUBDIRS += NIuTPCReader
SUBDIRS += NIuTPCMerger
SUBDIRS += NIuTPCMonitor
SUBDIRS += NIuTPCLogger
SUBDIRS += Dispatcher

SUBDIRS += NIuTPCReader0
SUBDIRS += NIuTPCReader1
SUBDIRS += NIuTPCReader2
SUBDIRS += NIuTPCReader3


.PHONY: $(SUBDIRS)

all: $(SUBDIRS)
	@set -e; for dir in $(SUBDIRS); do $(MAKE) -C $${dir} $@; done

clean:
	@set -e; for dir in $(SUBDIRS); do $(MAKE) -C $${dir} $@; done
