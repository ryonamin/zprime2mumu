SUBDIRS := $(wildcard src/*)
all:
	@for dir in ${SUBDIRS} ; do \
	(cd $$dir && ${MAKE}) ;\
	(cd ${BASEDIR}) ;\
        done
clean:
	@for dir in ${SUBDIRS} ; do \
	(cd $$dir && ${MAKE} clean) ;\
        done
