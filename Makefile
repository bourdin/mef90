LIBS = MEF90 VarStruct
PROJECTS = PrepVarFrac 
DIMPROJECTS = VarFracQS
SUBDIRS = ${LIBS} ${PROJECTS}

.PHONY: all ${SUBDIRS}

all: ${DIMPROJECTS}

${SUBDIRS}:
	cd $@ ; exec ${MAKE}

${DIMPROJECTS}: ${SUBDIRS}
	cd $@ ; exec ${MAKE} DIM=2D
	cd $@ ; exec ${MAKE} DIM=3D

clean:
	@for dir in ${SUBDIRS} ; do \
		${MAKE} -C $$dir clean ;\
	done


