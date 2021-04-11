SUBDIRS += NIuTPCReader
SUBDIRS += NIuTPCMerger
SUBDIRS += NIuTPCMonitor
SUBDIRS += NIuTPCLogger
SUBDIRS += Dispatcher

.PHONY: $(SUBDIRS)

all: $(SUBDIRS)
	@set -e; for dir in $(SUBDIRS); do $(MAKE) -C $${dir} $@; done

clean:
	@set -e; for dir in $(SUBDIRS); do $(MAKE) -C $${dir} $@; done
