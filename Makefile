LIBS = MEF90 m_VarFrac_Struct
PROJECTS = PrepVarFrac
SUBDIRS = ${LIBS} ${PROJECTS}

.PHONY: all ${SUBDIRS}

all: ${SUBDIRS}

${SUBDIRS}:
	cd $@ ; exec ${MAKE}

clean:
	@for dir in ${SUBDIRS} ; do \
		${MAKE} -C $$dir clean ;\
	done


