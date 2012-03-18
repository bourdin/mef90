LIBS = MEF90 VarStruct
PROJECTS = PrepVarFrac PrepVarFracNG PrepFilm VarFilmQS 
DIMPROJECTS = PoissonHeat VarFracQS
SUBDIRS = ${LIBS} ${PROJECTS}

.PHONY: all ${SUBDIRS} ${DIMPROJECTS}

help: 
	@echo Makefile system for MEF90 
	@echo usage make proj 
	@echo proj in clean, Frac, Film, Heat 

all: ${DIMPROJECTS}

Frac : MEF90 VarStruct PrepVarFrac PrepVarFracNG VarFracQS 

Film : MEF90 VarStruct PrepFilm VarFilmQS

Heat : MEF90 VarStruct PrepVarFracNG PoissonHeat VarFracQSHeat 

${SUBDIRS}:
	cd $@ ; exec ${MAKE}

${DIMPROJECTS}: 
	cd $@ ; exec ${MAKE} DIM=2D
	cd $@ ; exec ${MAKE} DIM=3D

VarFracQSHeat:
	cd VarFracQS ; exec ${MAKE} DIM=2D OPTS=Heat
	cd VarFracQS ; exec ${MAKE} DIM=3D OPTS=Heat


clean:
	-rm -rf ${MEF90_DIR}/${PETSC_ARCH}/obj
	-rm -rf ${MEF90_DIR}/${PETSC_ARCH}/bin
	-rm -rf ${MEF90_DIR}/${PETSC_ARCH}/lib
#	@for dir in ${SUBDIRS} ; do \
#		${MAKE} -C $$dir clean ;\
#	done


