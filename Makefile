LIBS = MEF90 m_VarFrac_Struct
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


