LIBS = MEF90 VarStruct
PROJECTS = PrepVarFrac PrepVarFracNG PrepFilm VarFilmQS 
DIMPROJECTS = VarFracQS
SUBDIRS = ${LIBS} ${PROJECTS}

.PHONY: all ${SUBDIRS}

help: 
	@echo Makefile system for MEF90 
	@echo usage make proj 
	@echo proj in Frac, FracNG, Film

all: ${DIMPROJECTS}

Frac : MEF90 VarStruct PrepVarFrac VarFracQS 

FracNG :  MEF90 VarStruct PrepVarFracNG VarFracQS 

Film : MEF90 VarStruct PrepFilm VarFilmQS

${SUBDIRS}:
	cd $@ ; exec ${MAKE}

${DIMPROJECTS}: 
	cd $@ ; exec ${MAKE} DIM=2D
	cd $@ ; exec ${MAKE} DIM=3D

clean:
	@for dir in ${SUBDIRS} ; do \
		${MAKE} -C $$dir clean ;\
	done


